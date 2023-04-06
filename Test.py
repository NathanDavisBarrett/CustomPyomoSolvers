import pyomo.environ as pyo
import numpy as np
import NathansSimplexSolver
from copy import deepcopy

percErrorTollerance = 0.01 #Percent

def percError(val1,val2):
    if abs(val1) < 1e-7:
        if abs(val2) < 1e-7:
            return 0
        else:
            return np.inf

    return abs((val1 - val2) / val1) * 100

model = pyo.ConcreteModel()

model.x = pyo.Var(domain=pyo.NonNegativeReals)
model.y = pyo.Var(domain=pyo.NonNegativeReals)

x = model.x
y = model.y

model.c1 = pyo.Constraint(expr= y <= x+1)
model.c2 = pyo.Constraint(expr= y >= -x+1.2)
model.c3 = pyo.Constraint(expr= y >= 1.2*x-0.5)
model.c4 = pyo.Constraint(expr= y <= -0.5*x + 3)
model.c5 = pyo.Constraint(expr= y == 1.1*x+0.5)

model.obj = pyo.Objective(expr= x + y,sense=pyo.maximize) #Fails only on maximize...

solver = pyo.SolverFactory('NathansSimplexSolver')
result = solver.solve(model)#.write()

myX = deepcopy(model.x.value)
myY = deepcopy(model.y.value)
myObj = deepcopy(model.obj())

model.x.display()
model.y.display()

solver = pyo.SolverFactory('gurobi')
result = solver.solve(model)

model.x.display()
model.y.display()

assert percError(model.x.value, myX) < percErrorTollerance
assert percError(model.y.value, myY) < percErrorTollerance
assert percError(model.obj(), myObj) < percErrorTollerance