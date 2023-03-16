#From https://stackoverflow.com/questions/74179391/custom-python-solver-for-pyomo

import pyomo.environ as pyo
from pyomo.opt import SolverFactory

class GenericSolverInterface(object):
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
        pass

    def available(self, exception_flag=True):#, -> bool:
        pass

    def license_is_valid(self):# -> bool:
        pass

    def version(self):# -> Tuple:
        pass

    @property
    def options(self):# -> ConfigDict:
        pass

    @options.setter
    def options(self, val):
        pass            

# class DemoSolver(GenericSolverInterface):
#     pass

# SolverFactory.register('demo', doc='DEMO Solver Interface')(DemoSolver)

from pyomo.core.base import Objective, Constraint
from pyomo.core.expr.visitor import identify_variables

class model_interface(object):
    def __init__(self, model):
        self.model = model
        self._obj = list(model.component_data_objects(Objective, active=True))
        self._con = list(model.component_data_objects(Constraint, active=True))
        self._var = {}
        for c in self._con:
            self._var.update({id(v): v for v in identify_variables(c.body)})
        self._var = list(self._var.values())

    @property
    def x(self):
        """Return the current list of variable values"""
        return [v.value for v in self._var]

    @x.setter
    def x(self, values):
        """Set the variables to new values"""
        for v, val in zip(self._var, values):
            v.set_value(val)

    @property
    def x_lb(self):
        """Return the list of variable lower bounds (may include None)"""
        return [v.lb for v in self._var]

    @property
    def x_ub(self):
        """Return the list of variable upper bounds (may include None)"""
        return [v.ub for v in self._var]

    @property
    def obj(self):
        """Return the list of objective values (computed using the current
        variable values)"""
        return [value(o) for o in self._obj]

    @property
    def con(self):
        """Return the list of constraint 'body' values (computed using the
        current variable values)
        """
        return [value(c) for c in self._con]

    @property
    def con_lb(self):
        """Return the list of constraint lower bounds (may include None)"""
        return [c.lb for c in self._con]

    @property
    def con_ub(self):
        """Return the list of constraint upper bounds (may include None)"""
        return [c.ub for c in self._con]