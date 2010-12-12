#doorbuiging van een plank

nx = 50
dx = 1.
from fipy.meshes.grid1D import Grid1D
mesh = Grid1D(nx = nx, dx = dx)

from fipy.variables.cellVariable import CellVariable
u = CellVariable(name="solution variable",
                    mesh=mesh,
                    value=0,hasOld=1)

from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.transientTerm import TransientTerm
a=1

timeStepDuration = 0.01
steps=1000

eq= ImplicitDiffusionTerm(coeff=(a,a))+(u.getOld().getOld() - 2*u.getOld() + u) / timeStepDuration ** 2




from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition

BCs = ( NthOrderBoundaryCondition(faces=mesh.getFacesRight(), value=0,order=2),
	NthOrderBoundaryCondition(faces=mesh.getFacesRight(), value=1,order=3),
	FixedValue(faces=mesh.getFacesLeft(), value=0),
	FixedFlux(faces=mesh.getFacesLeft(), value=0))

from fipy import viewers
viewer = viewers.make(vars=(u))

def hit_continue(Prompt='Hit any key to continue'):
	raw_input(Prompt)


for step in range(steps):
	u.updateOld()
	eq.solve(var=u,boundaryConditions=BCs,dt=timeStepDuration)
	eq= ImplicitDiffusionTerm(coeff=(a,a))+(u.getOld().getOld() - 2*u.getOld() + u) / timeStepDuration ** 2
	viewer.plot()
	hit_continue()

#getOld seems to be not working
