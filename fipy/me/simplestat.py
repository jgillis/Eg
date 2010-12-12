#Laplacian=0

nx = 50
dx = 1.
from fipy.meshes.grid1D import Grid1D
mesh = Grid1D(nx = nx, dx = dx)

from fipy.variables.cellVariable import CellVariable
phi = CellVariable(name="solution variable",
                    mesh=mesh,
                    value=0)

z = CellVariable(name="dummy",
                    mesh=mesh,
                    value=0)

valueLeft = 1
valueRight = 0

from fipy.boundaryConditions.fixedValue import FixedValue
BCs = (FixedValue(faces=mesh.getFacesRight(), value=valueRight),
        FixedValue(faces=mesh.getFacesLeft(), value=valueLeft))


from fipy.terms.explicitDiffusionTerm import ExplicitDiffusionTerm

eqX = ExplicitDiffusionTerm(coeff=1) 


from fipy import viewers
viewer = viewers.make(vars=(phi),
                      limits={'datamin': 0., 'datamax': 1.})
viewer.plot()
eqX.solve(var=phi,
        	boundaryConditions=BCs)
viewer.plot()


