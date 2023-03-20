import pyomo.environ as pyo
import numpy as np
import NathansSimplexSolver

model = pyo.ConcreteModel()

numI = 3
numJ = 4
numK = 5

model.I = pyo.Set(initialize=range(numI))
model.J = pyo.Set(initialize=range(numJ))
model.K = pyo.Set(initialize=range(numK))

model.I_J = pyo.Set(initialize=model.I * model.J)
model.I_J_K = pyo.Set(initialize=model.I * model.J * model.K)

model.x = pyo.Var(model.I_J_K, domain=pyo.NonNegativeReals)

c = np.random.uniform(0,1,(numI,numJ,numK))
b = np.random.uniform(0,1,(numI,numJ))

model.objFunc = pyo.Objective(
    expr = sum([c[i,j,k] * model.x[i,j,k] for i in model.I for j in model.J for k in model.K]),
    sense = pyo.maximize
)

def myConstraint(m,i,j):
    return sum([m.x[i,j,k] for k in m.K]) <= b[i,j]
model.myConstraint = pyo.Constraint(model.I_J, rule=myConstraint)

solver = pyo.SolverFactory('NathansSimplexSolver')
solver.solve(model)#.write()

#model.objFunc.display()
print("ObjFunc = {}".format(model.objFunc()))
model.x.display()

