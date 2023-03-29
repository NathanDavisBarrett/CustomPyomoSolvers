import CustomSolverResources
import pyomo.environ as pyo
from pyomo.opt import SolverFactory

import random
import string
import numpy as np

import nathans_Simplex_solver_py

class SolverResults:
        def __init__(self,logs,finalStatus,numIter,time,objValue):
            self.logs = logs
            self.finalStatus = finalStatus
            self.numIter = numIter
            self.time = time
            self.objValue = objValue

        def write(self):
            if len(self.logs) > 0:
                print("########## Solver Messages ##########")
                for _,message in self.logs:
                    print("\t{}".format(message))
                print("#####################################")

            print("############## RESULTS ##############")
            print("Final Status         :",self.finalStatus)
            print("Number of Iterations :",self.numIter)
            print("Execution Time       :",self.time)
            print("Final Objective Value:",self.objValue)
            print("#####################################")

class NathansSimplexSolver(CustomSolverResources.GenericSolverInterface):
    def __init__(self):
        self.availableOptions = {
            "maxIter": ["AnyPositiveInt"],
            "maxTime": ["AnyPositiveFloat"],
            "LiveUpdate": [False,"iter","time"],
            "LiveUpdateIter": ["AnyPositiveFloat"],
            "LiveUpdateTime": ["AnyPositiveFloat"],
            "BasisReductionMethod": ["AVERAGE","FIRST"]
        }
        self.currentOptions = {
            "maxIter": np.iinfo(np.int64).max,
            "maxTime": np.inf,
            "LiveUpdate": False,
            "LiveUpdateIter": np.iinfo(np.int64).max,
            "LiveUpdateTime": np.inf,
            "BasisReductionMethod": "FIRST"
        }

        self.codes = {
            "logLevels": {
                0: "ERROR",
                1: "WARN",
                2: "INFO"
            },
            "basisReductionOptions": {
                0: "AVERAGE",
                1: "FIRST",
                100: "INVALID"
            },
            "finalStatus": {
                0: "OPTIMAL_SOLUTION_FOUND",
                1: "STOPPED_AT_MAX_ITER",
                2: "STOPPED_AT_MAX_TIME",
                3: "UNBOUNDED",
                4: "INFEASIBLE",
                5: "UNDETERMINED_STATUS",
                1234567890: "TERMINATE"
            },
            "constraintTypes": {
                0: "EQ",
                1: "LEQ",
                2: "GEQ"
            }
        }
        categories = list(self.codes.keys())
        for category in categories:
            keys = list(self.codes[category].keys())
            for key in keys:
                val = self.codes[category][key]
                self.codes[category][val] = key

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
        elif "GeneralVarData" in exprType or "ScalarVar" in exprType:
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
        varIDs = {varNames[i]: i for i in range(len(varNames))}

        for var in modelData._var:
            if not var.is_continuous():
                raise Exception("Error! The variable {} is not continuous. This Simplex Solver can only handle continuous variables. Did you mean to call NathansMILPSolver?")

        constrExpressions = [str(constr.expr) for constr in modelData._con]

        constrTypes = np.zeros(len(constrExpressions),dtype=int)

        for i in range(len(constrExpressions)):
            if "<=" in constrExpressions[i]:
                constrTypes[i] = self.codes["constraintTypes"]["LEQ"]
            elif ">=" in constrExpressions:
                constrTypes[i] = self.codes["constraintTypes"]["GEQ"]
            else:
                constrTypes[i] = self.codes["constraintTypes"]["EQ"]
                

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

            if constData < 0:
                #Standard form has no negative in the constant column
                constData *= -1
                varData *= -1
                if constrTypes[i] == self.codes["constraintTypes"]["LEQ"]:
                    constrTypes[i] = self.codes["constraintTypes"]["GEQ"]
                elif constrTypes[i] == self.codes["constraintTypes"]["GEQ"]:
                    constrTypes[i] = self.codes["constraintTypes"]["LEQ"]

            A[i,:] = varData
            b[i] = constData

        if len(modelData._obj) != 1:
            raise Exception("Error! Haven't programed how to handle multiple objective functions yet.")

        recurseC = self._recurseExprTree(modelData._obj[0].expr,varNames)

        c = np.append(recurseC[0],recurseC[1])

        #Collect all variable bounds.
        bounds = [tuple([b if b != None else 0 for b in var.bounds]) for var in modelData._var]
        boundActivations = [tuple([b != None for b in var.bounds]) for var in modelData._var]

        maximizationProblem = modelData._obj[0].sense == pyo.maximize

        #Now pass these matrices to C++        
        solver = nathans_Simplex_solver_py.SimplexSolver()
        solver.setMaxIter(self.currentOptions["maxIter"])
        solver.setMaxTime(self.currentOptions["maxTime"])

        if self.currentOptions["LiveUpdate"]:
            solver.SetLiveUpdateSettings(self.currentOptions["LiveUpdateIter"],self.currentOptions["LiveUpdateTime"])

        solver.EngageModel(numVar,numConstr,varNames,A.ravel(),b,c,constrTypes,bounds,boundActivations,maximizationProblem)
        exitCode = solver.Solve()

        logs = solver.GetLogs()
        finalStatus = self.codes["finalStatus"][exitCode]
        numIter = solver.GetNumIterations()
        time = solver.GetSolveTime()
        objValue = solver.GetObjectiveValue()
        
        basisReductionMethodCode = self.codes["basisReductionOptions"][self.currentOptions["BasisReductionMethod"]]
        
        for var in modelData._var:
            var.value = solver.GetVariableValue(var.name,basisReductionMethodCode)

        results = SolverResults(logs,finalStatus,numIter,time,objValue)

        return results

SolverFactory.register('NathansSimplexSolver', doc='Basic Simplex Solver by Nathan Barrett')(NathansSimplexSolver)

class NathansMILPSolver(CustomSolverResources.GenericSolverInterface):
    pass

SolverFactory.register('NathansMILPSolver', doc='Basic MILP Solver by Nathan Barrett')(NathansMILPSolver)