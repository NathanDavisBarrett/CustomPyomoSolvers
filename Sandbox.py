import pyomo.environ as pyo
import numpy as np
import NathansSimplexSolver

model1 = pyo.ConcreteModel()

numI = 3
numJ = 4
numK = 5

model1.I = pyo.Set(initialize=range(numI))
model1.J = pyo.Set(initialize=range(numJ))
model1.K = pyo.Set(initialize=range(numK))

model1.I_J = pyo.Set(initialize=model1.I * model1.J)
model1.I_J_K = pyo.Set(initialize=model1.I * model1.J * model1.K)

model1.x = pyo.Var(model1.I_J_K, domain=pyo.NonNegativeReals)

c = np.random.uniform(0,1,(numI,numJ,numK))
b = np.random.uniform(0,1,(numI,numJ))

model1.objFunc = pyo.Objective(
    expr = sum([c[i,j,k] * model1.x[i,j,k] for i in model1.I for j in model1.J for k in model1.K]),
    sense = pyo.maximize
)

def myConstraint(m,i,j):
    return sum([m.x[i,j,k] for k in m.K]) <= b[i,j]
model1.myConstraint = pyo.Constraint(model1.I_J, rule=myConstraint)


model2 = pyo.ConcreteModel()

model2.x = pyo.Var(range(3),domain=pyo.NonNegativeReals)
model2.y = pyo.Var(range(2),domain=pyo.NonNegativeReals)#Binary)
model2.c1 = pyo.Constraint(expr=  model2.x[0]                               + 13.2 * model2.y[0] - 0.001 * model2.y[1] <= 1 )
model2.c2 = pyo.Constraint(expr=4*model2.x[0] +   model2.x[1] + 6*model2.x[2] - 43.1 * model2.y[0] + 5.2  * model2.y[1] <= 6 )
model2.c3 = pyo.Constraint(expr=8*model2.x[0] + 4*model2.x[1] +   model2.x[2]                     + 7    * model2.y[1] <= 36)
model2.c4 = pyo.Constraint(expr= model2.x[0] == 1.2 * model2.y[0] + 0.01 * model2.y[1] + 1.0)

model2.objFunc = pyo.Objective(expr=4 * model2.x[0] + 2 * model2.x[1] + model2.x[2], sense = pyo.maximize)


model3 = pyo.ConcreteModel()

model3.x = pyo.Var(range(3),domain=pyo.NonNegativeReals)
model3.c1 = pyo.Constraint(expr=  model3.x[0]                                 >= 1 )
model3.c2 = pyo.Constraint(expr=4*model3.x[0] +   model3.x[1] + 6*model3.x[2] <= 6 )
model3.c3 = pyo.Constraint(expr=8*model3.x[0] + 4*model3.x[1] +   model3.x[2] <= 36)
#model3.c4 = pyo.Constraint(expr=  model3.x[0]                                 == 1 )

model3.objFunc = pyo.Objective(expr=4 * model3.x[0] + 2 * model3.x[1] + model3.x[2], sense = pyo.maximize)


model = model3

solver = pyo.SolverFactory('NathansSimplexSolver')
solver.options = {
    "LiveUpdateTime": 10,
    "maxIter": 25
}

# solver = pyo.SolverFactory('gurobi')

solver.solve(model).write()

#model.objFunc.display()
print("ObjFunc = {}".format(model.objFunc()))
model.x.display()
model.y.display()
model.myConstraint.display()

