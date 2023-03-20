import CustomSolverResources
import pyomo.environ as pyo
from pyomo.opt import SolverFactory

import random
import string
import numpy as np

import nathans_Simplex_solver_py

class NathansSimplexSolver(CustomSolverResources.GenericSolverInterface):
    def __init__(self):
        self.availableOptions = {
            "Basis_IO_Approach": ["MaximizeRC"],
            "maxIter": ["AnyPositiveInt"],
            "maxTime": ["AnyPositiveFloat"],
            "LiveUpdate": [False,"iter","time"],
            "LiveUpdateIter": ["AnyPositiveFloat"],
            "LiveUpdateTime": ["AnyPositiveFloat"]
        }
        self.currentOptions = {
            "Basis_IO_Approach": "MaximizeRC",
            "maxIter": np.iinfo(np.int64).max,
            "maxTime": np.inf,
            "LiveUpdate": False,
            "LiveUpdateIter": np.iinfo(np.int64).max,
            "LiveUpdateTime": np.inf
        }

    def _recurseExprTree(self,expr,varNames,multiplier=1):
        coefs = np.zeros(len(varNames))
        const = 0

        exprType = str(type(expr))

        if isinstance(expr,int) or isinstance(expr,float):
            return [coefs,expr]
        elif "ProductExpression" in exprType or "MonomialTermExpression" in exprType:
            #In order for this to be a linear problem, one of the child nodes here needs to be a number.
            isNumeric = [isinstance(arg,int) or isinstance(arg,float) for arg in expr.args]
            if sum(isNumeric) < len(expr.args)-1: #Only one is allowed to be non-numeric
                raise Exception("Error! Received non-numeric terms on both sides of a multiplication operator: \"{}\". The problem is likely nonlinear and cannot be solved with simplex.".format(expr))

            if len(expr.args) != 2:
                raise Exception("Error! I haven't coded how to handle product expressions with more than two terms.")

            if isNumeric[0]:
                return self._recurseExprTree(expr.args[1],varNames,multiplier*expr.args[0])
            else:
                return self._recurseExprTree(expr.args[0],varNames,multiplier*expr.args[1])
        elif "SumExpression" in exprType:
            for arg in expr.args:
                result = self._recurseExprTree(arg,varNames,1)
                coefs += result[0]
                const += result[1]

            return [coefs * multiplier, const * multiplier]
        elif "GeneralVarData" in exprType:
            varName = str(expr)
            varIndex = varNames.index(varName)
            coefs[varIndex] = multiplier
            
            return [coefs,const]
        else:
            raise Exception("{} is not a recognized expression subtype. It probably just needs to be coded in here.".format(exprType))

    def version(self):# -> Tuple:
        return (0,0,1)

    @property
    def options(self):# -> ConfigDict:
        return self.currentOptions

    @options.setter
    def options(self, val):
        errorMessage = "Options can only be set to a dict object."

        if type(val) != dict:
            raise Exception(errorMessage)

        if "LiveUpdateIter" in val or "LiveUpdateTime" in val:
            self.currentOptions["LiveUpdate"] = True

        for key in val:
            if type(key) != str:
                raise Exception(errorMessage)
            if type(val[key]) not in [int,float,str]:
                raise Exception(errorMessage)

            if key in ["maxIter","LiveUpdateIter"]:
                if val[key] < 0:
                    raise Exception("{} must be a positive integer.".format(key))
                self.currentOptions[key] = int(val[key])
            elif key in ["maxTime","LiveUpdateTime"]:
                if val[key] < 0:
                    raise Exception("{} must be a positive float.".format(key))
                self.currentOptions[key] = float(val[key])
            else:
                if key in self.availableOptions.keys():
                    if val[key] in self.availableOptions[key]:
                        self.currentOptions[key] = val[key]
                    else:
                        raise Exception("{} is not a valid option for {}".format(val[key],key))
                else:
                    raise Exception("{} is not a recognized option".format(key))

    def solve(self,
              model,#: _BlockData,
              tee=False,#: bool = False,
              load_solutions=True,#: bool = True,
              logfile=None,#: Optional[str] = None,
              solnfile=None,#: Optional[str] = None,
              timelimit=None,#: Optional[float] = None,
              report_timing=False,#: bool = False,
              solver_io=None,#: Optional[str] = None,
              suffixes=None,#: Optional[Sequence] = None,
              options=None,#: Optional[Dict] = None,
              keepfiles=False,#: bool = False,
              symbolic_solver_labels=False#: bool = False):
              ):
        modelData = CustomSolverResources.model_interface(model)

        varNames = [var.name for var in modelData._var]
        constrExpressions = [str(constr.expr) for constr in modelData._con]

        #Add Slack Variables
        slackIndex = 0
        slackVarBase = "SLACK_VAR"
        while not all([slackVarBase not in varName for varName in varNames]): #Make sure the slack var name is unique.
            slackVarBase += random.choice(string.ascii_letters)

        slackAdditions = [None for constr in constrExpressions]

        for i in range(len(constrExpressions)):
            if "<=" in constrExpressions[i]:
                newSlackVar = slackVarBase + "[{}]".format(slackIndex)
                varNames.append(newSlackVar)
                slackIndex += 1

                slackAdditions[i] = (1,newSlackVar)

            elif ">=" in constrExpressions:
                newSlackVar = slackVarBase + "[{}]".format(slackIndex)
                varNames.append(newSlackVar)
                slackIndex += 1

                slackAdditions[i] = (-1,newSlackVar)
                

        #Assemble Matrices
        numVar = len(varNames)
        numConstr = len(constrExpressions)

        A = np.zeros((numConstr,numVar))
        b = np.zeros(numConstr)

        for i in range(len(constrExpressions)):
            constrBody = modelData._con[i].expr

            leftData = self._recurseExprTree(constrBody.args[0],varNames)
            rightData = self._recurseExprTree(constrBody.args[1],varNames)

            varData = leftData[0] - rightData[0]
            constData = rightData[1] - leftData[1]

            if slackAdditions[i] != None:
                contribution, varName = slackAdditions[i]
                varIndex = varNames.index(varName)
                varData[varIndex] = contribution

            A[i,:] = varData
            b[i] = constData

        if len(modelData._obj) != 1:
            raise Exception("Error! Haven't programed how to handle multiple objective functions yet.")

        recurseC = self._recurseExprTree(modelData._obj[0].expr,varNames)

        c = np.append(recurseC[0],recurseC[1])
        
        #Now pass these matrices to C++
        basis = []
        for _,varName in slackAdditions:
            basis.append(varNames.index(varName))
        
        if self.currentOptions["Basis_IO_Approach"] == "MaximizeRC":
            solver = nathans_Simplex_solver_py.SimplexSolver_MaximizeRC()
            solver.setMaxIter(self.currentOptions["maxIter"])
            solver.setMaxTime(self.currentOptions["maxTime"])
        else:
            raise Exception("The Basis_IO_Approach \"{}\" went un-initialized. This is a bug.")

        if self.currentOptions["LiveUpdate"]:
            solver.SetLiveUpdateSettings(self.currentOptions["LiveUpdateIter"],self.currentOptions["LiveUpdateTime"])


        solver.EngageModel(numVar,numConstr,varNames,A.flatten(),b,c,basis)
        solver.Solve()

        logs = solver.GetLogs()
        logLevels = ["ERROR"," WARN"," INFO"]
        if len(logs) > 0:
            print("Solver Messages:")
        for level,message in logs:
            print("\t{}: {}".format(logLevels[level],message))
        
        for var in modelData._var:
            var.value = solver.GetVariableValue(var.name)

        
        

            





        


SolverFactory.register('NathansSimplexSolver', doc='Basic Simplex Solver by Nathan Barrett')(NathansSimplexSolver)