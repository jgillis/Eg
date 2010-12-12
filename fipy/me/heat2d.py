def hit_continue(Prompt='Hit any key to continue'):
	raw_input(Prompt)

nx = 50
ny = nx
dx = 1.
dy = dx
from fipy.meshes.grid2D import Grid2D
mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)

from fipy.variables.cellVariable import CellVariable
phi = CellVariable(name = "solution variable",mesh = mesh,value = 0)

D = 1.
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
eq = ImplicitDiffusionTerm(coeff=D)

valueLeft = 1
valueRight = 0

from fipy.boundaryConditions.fixedValue import FixedValue
BCs = (FixedValue(faces=mesh.getFacesLeft(),value=0),
       FixedValue(faces=mesh.getFacesRight(),value=1),
	FixedValue(faces=mesh.getFacesTop(),value=0),
	FixedValue(faces=mesh.getFacesBottom(),value=0))


from fipy import viewers
viewer = viewers.make(vars=phi,
                       limits={'datamin': 0., 'datamax': 1.})
viewer.plot()


timeStepDuration = 10 * 0.9 * dx**2 / (2 * D)
steps = 10

eq.solve(var=phi,boundaryConditions=BCs)
viewer.plot()



hit_continue()

