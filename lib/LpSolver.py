from gurobipy import GurobiError 

class LpSolver(object): 
	def __init__(self, model, out=None, gurobi_out=''):
		if out is None:
			out = sys.stdout
		self._model = model
		self._debug_fn = debug_fn
		self.out = out	
		self_gurobi_out = gurobi_out

	def solve_lp(self, method=Method.CONCURRENT):
		model = self.model
		model.setParam('Method', method.value)
		model.setParam('LogFile', self.gurobi_out)
		try: 
			model.optimize()
			return model.objVal
				
