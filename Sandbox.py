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

################ TEST 1 ################
########################################
#The initial case this solver was debugged with
# From: https://www.youtube.com/watch?v=-RtEzwfMqxk&list=PLg2tfDG3Ww4vyVtIvTUY2JaOZDbQStcsb&index=9
model = pyo.ConcreteModel()

model.x = pyo.Var(domain=pyo.NonNegativeReals)
model.y = pyo.Var(domain=pyo.NonNegativeReals)

model.c1 = pyo.Constraint(expr= model.x + model.y <= 15)
model.c2 = pyo.Constraint(expr= 2*model.x + 3 * model.y <= 40)
model.c3 = pyo.Constraint(expr= model.x + 8 * model.y >= 18)
model.c4 = pyo.Constraint(expr= model.x >= 2)

model.objFunc = pyo.Objective(expr= 5*model.x + 6 * model.y, sense = pyo.maximize)

solver = pyo.SolverFactory('NathansSimplexSolver')
result = solver.solve(model)#.write()

myX = deepcopy(model.x.value)
myY = deepcopy(model.y.value)
myObj = deepcopy(model.objFunc())

solver = pyo.SolverFactory('gurobi')
result = solver.solve(model)

assert percError(model.x.value, myX) < percErrorTollerance
assert percError(model.y.value, myY) < percErrorTollerance
assert percError(model.objFunc(), myObj) < percErrorTollerance
########################################

################ TEST 2 ################
########################################
#Demonstrating basic LEQ constraints
# From: https://math.libretexts.org/Bookshelves/Applied_Mathematics/Applied_Finite_Mathematics_(Sekhon_and_Bloom)/04%3A_Linear_Programming_The_Simplex_Method/4.02%3A_Maximization_By_The_Simplex_Method
model = pyo.ConcreteModel()

model.x1 = pyo.Var(domain=pyo.NonNegativeReals)
model.x2 = pyo.Var(domain=pyo.NonNegativeReals)

model.c1 = pyo.Constraint(expr= model.x1 + model.x2 <= 12)
model.c2 = pyo.Constraint(expr= 2 * model.x1 + model.x2 <= 16)

model.obj = pyo.Objective(expr= 40 * model.x1 + 30 * model.x2, sense= pyo.maximize)

solver = pyo.SolverFactory('NathansSimplexSolver')
result = solver.solve(model)#.write()

myX1 = deepcopy(model.x1.value)
myX2 = deepcopy(model.x2.value)
myObj = deepcopy(model.obj())

solver = pyo.SolverFactory('gurobi')
result = solver.solve(model)

assert percError(model.x1.value, myX1) < percErrorTollerance
assert percError(model.x2.value, myX2) < percErrorTollerance
assert percError(model.obj(), myObj) < percErrorTollerance
########################################

################ TEST 3 ################
########################################
#Demonstrating mixed LEQ and GEQ constraints
# From: https://college.cengage.com/mathematics/larson/elementary_linear/4e/shared/downloads/c09s5.pdf
model = pyo.ConcreteModel()

model.x1 = pyo.Var(domain=pyo.NonNegativeReals)
model.x2 = pyo.Var(domain=pyo.NonNegativeReals)
model.x3 = pyo.Var(domain=pyo.NonNegativeReals)

model.c1 = pyo.Constraint(expr= 3*model.x1 + 2*model.x2 + 5*model.x3 <= 18)
model.c2 = pyo.Constraint(expr= 4*model.x1 + 2*model.x3 + 3*model.x3 <= 16)
model.c3 = pyo.Constraint(expr= 2*model.x1 + model.x2 + model.x3 >= 4)

model.obj = pyo.Objective(expr= 3 * model.x1 + 2 * model.x2 + 4*model.x3, sense= pyo.maximize)

solver = pyo.SolverFactory('NathansSimplexSolver')
result = solver.solve(model)#.write()

myX1 = deepcopy(model.x1.value)
myX2 = deepcopy(model.x2.value)
myX3 = deepcopy(model.x3.value)
myObj = deepcopy(model.obj())

solver = pyo.SolverFactory('gurobi')
result = solver.solve(model)

assert percError(model.x1.value, myX1) < percErrorTollerance
assert percError(model.x2.value, myX2) < percErrorTollerance
assert percError(model.x3.value, myX3) < percErrorTollerance
assert percError(model.obj(), myObj) < percErrorTollerance

########################################

################ TEST 4 ################
########################################
#Demonstrating mixed GEQ, LEQ, and EQ constraints
#Made up by Nathan
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

model.obj = pyo.Objective(expr= x + y)

solver = pyo.SolverFactory('NathansSimplexSolver')
result = solver.solve(model)#.write()

myX = deepcopy(model.x.value)
myY = deepcopy(model.y.value)
myObj = deepcopy(model.obj())

solver = pyo.SolverFactory('gurobi')
result = solver.solve(model)

assert percError(model.x.value, myX) < percErrorTollerance
assert percError(model.y.value, myY) < percErrorTollerance
assert percError(model.obj(), myObj) < percErrorTollerance

########################################

################ TEST 5 ################
########################################
#Demonstrating mixed GEQ, LEQ, and EQ constraints as well as definite upper bounds
#Made up by Nathan

model = pyo.ConcreteModel()

model.x = pyo.Var(domain=pyo.NonNegativeReals,bounds=(None,1.1))
model.y = pyo.Var(domain=pyo.NonNegativeReals)

x = model.x
y = model.y

model.c1 = pyo.Constraint(expr= y <= x+1)
model.c2 = pyo.Constraint(expr= y >= -x+1.2)
model.c3 = pyo.Constraint(expr= y >= 1.2*x-0.5)
model.c4 = pyo.Constraint(expr= y <= -0.5*x + 3)
model.c5 = pyo.Constraint(expr= y == 1.1*x+0.5)

model.obj = pyo.Objective(expr= x + y)

solver = pyo.SolverFactory('NathansSimplexSolver')
result = solver.solve(model)#.write()

myX = deepcopy(model.x.value)
myY = deepcopy(model.y.value)
myObj = deepcopy(model.obj())

solver = pyo.SolverFactory('gurobi')
result = solver.solve(model)

assert percError(model.x.value, myX) < percErrorTollerance
assert percError(model.y.value, myY) < percErrorTollerance
assert percError(model.obj(), myObj) < percErrorTollerance
########################################

################ TEST 6 ################
########################################
#Demonstrating mixed GEQ, LEQ, and EQ constraints as well as definite lower bounds

#Made up by Nathan

model = pyo.ConcreteModel()

model.x = pyo.Var(domain=pyo.NonNegativeReals,bounds=(1,None))
model.y = pyo.Var(domain=pyo.NonNegativeReals)

x = model.x
y = model.y

model.c1 = pyo.Constraint(expr= y <= x+1)
model.c2 = pyo.Constraint(expr= y >= -x+1.2)
model.c3 = pyo.Constraint(expr= y >= 1.2*x-0.5)
model.c4 = pyo.Constraint(expr= y <= -0.5*x + 3)
model.c5 = pyo.Constraint(expr= y == 1.1*x+0.5)

model.obj = pyo.Objective(expr= x + y,sense=pyo.maximize)

solver = pyo.SolverFactory('NathansSimplexSolver')
result = solver.solve(model)#.write()

myX = deepcopy(model.x.value)
myY = deepcopy(model.y.value)
myObj = deepcopy(model.obj())

solver = pyo.SolverFactory('gurobi')
result = solver.solve(model)

assert percError(model.x.value, myX) < percErrorTollerance
assert percError(model.y.value, myY) < percErrorTollerance
assert percError(model.obj(), myObj) < percErrorTollerance

########################################

################ TEST 7 ################
########################################
#Demonstrating mixed GEQ, LEQ, and EQ constraints as well as both definite upper and lower bounds
model = pyo.ConcreteModel()

model.x = pyo.Var(domain=pyo.NonNegativeReals,bounds=(0.5,1.1))
model.y = pyo.Var(domain=pyo.NonNegativeReals)

x = model.x
y = model.y

model.c1 = pyo.Constraint(expr= y <= x+1)
model.c2 = pyo.Constraint(expr= y >= -x+1.2)
model.c3 = pyo.Constraint(expr= y >= 1.2*x-0.5)
model.c4 = pyo.Constraint(expr= y <= -0.5*x + 3)
model.c5 = pyo.Constraint(expr= y == 1.1*x+0.5)

model.obj = pyo.Objective(expr= x + y)

solver = pyo.SolverFactory('NathansSimplexSolver')
result = solver.solve(model)#.write()

myX = deepcopy(model.x.value)
myY = deepcopy(model.y.value)
myObj = deepcopy(model.obj())

solver = pyo.SolverFactory('gurobi')
result = solver.solve(model)

assert percError(model.x.value, myX) < percErrorTollerance
assert percError(model.y.value, myY) < percErrorTollerance
assert percError(model.obj(), myObj) < percErrorTollerance

########################################

################ TEST 8 ################
########################################
#Demonstrating mixed GEQ, LEQ, and EQ constraints as well as unbounded variables
model = pyo.ConcreteModel()

model.x = pyo.Var(domain=pyo.Reals)
model.y = pyo.Var(domain=pyo.NonNegativeReals)

x = model.x
y = model.y

model.c1 = pyo.Constraint(expr= y <= x+1)
model.c2 = pyo.Constraint(expr= y >= -x+1.2)
model.c3 = pyo.Constraint(expr= y >= 1.2*x-0.5)
model.c4 = pyo.Constraint(expr= y <= -0.5*x + 3)
model.c5 = pyo.Constraint(expr= y == 1.1*x+0.5)

model.obj = pyo.Objective(expr= x + y)

solver = pyo.SolverFactory('NathansSimplexSolver')
result = solver.solve(model)#.write()

myX = deepcopy(model.x.value)
myY = deepcopy(model.y.value)
myObj = deepcopy(model.obj())

solver = pyo.SolverFactory('gurobi')
result = solver.solve(model)

assert percError(model.x.value, myX) < percErrorTollerance
assert percError(model.y.value, myY) < percErrorTollerance
assert percError(model.obj(), myObj) < percErrorTollerance

########################################

################ TEST 9 ################
########################################
#Demonstrating mixed GEQ, LEQ, and EQ constraints as well as a combination of def. upper, lower, both, and neither bounds.

########################################

################ TEST 10 ################
#########################################
#Demonstrating the ability to detect and report unbounded problems

#########################################

################ TEST 11 ################
#########################################
#Demonstrating the ability to detect and report infeasible problems

#########################################

################ TEST 12 ################
#########################################
#Demonstrating the ability to handle minimization objectives.

#########################################

print("All Tests Passed!")