nx = 50
dx = 1.
from fipy.meshes.grid1D import Grid1D
mesh = Grid1D(nx = nx, dx = dx)

from fipy.variables.cellVariable import CellVariable
phi = CellVariable(name="solution variable",
                    mesh=mesh,
                    value=0)
D = 1

valueLeft = 1
valueRight = 0

from fipy.boundaryConditions.fixedValue import FixedValue
BCs = (FixedValue(faces=mesh.getFacesRight(), value=valueRight),
        FixedValue(faces=mesh.getFacesLeft(), value=valueLeft))


from fipy.terms.explicitDiffusionTerm import ExplicitDiffusionTerm
from fipy.terms.transientTerm import TransientTerm
eqX = TransientTerm() == ExplicitDiffusionTerm(coeff=D)


timeStepDuration = 0.9 * dx**2 / (2 * D)
steps=100



from fipy import viewers
viewer = viewers.make(vars=(phi),
                      limits={'datamin': 0., 'datamax': 1.})

def hit_continue(Prompt='Hit any key to continue'):
	raw_input(Prompt)

for step in range(steps):
	eqX.solve(var=phi,
          	boundaryConditions=BCs,
          	dt=timeStepDuration)
	viewer.plot()
	hit_continue()
	

