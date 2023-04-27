from enum import Enum
import sys

class Method(Enum):
    PRIMAL_SIMPLEX = 0
    DUAL_SIMPLEX = 1
    BARRIER = 2
    CONCURRENT = 3
    PRIMAL_AND_DUAL = 5


class LpSolver(object): 
    def __init__(self, model, out=None, gurobi_out=''):
        if out is None:
            out = sys.stdout
        self._model = model
        self.out = out  
        self.gurobi_out = gurobi_out

    def solve_lp(self, method=Method.BARRIER):
        model = self._model
        model.setParam('Method', method.value) # CONCURRENT
        model.setParam('LogFile', self.gurobi_out)
        try: 
            model.optimize()
            return model.objVal
        except:
            print("Model failed to optimize lp.")
