#!/usr/bin/env python
L = 1.0
N = 40
dL = L / N
viscosity = 1.
pressureRelaxation = 0.2
velocityRelaxation = 0.5

sweeps=10

from fipy.meshes.grid2D import Grid2D
mesh = Grid2D(nx=N, ny=N, dx=dL, dy=dL)

from fipy.variables.cellVariable import CellVariable

psi = CellVariable(mesh=mesh, name='stream function')

Xhat=(1,0)
Yhat=(0,1)

u=psi.getGrad().dot(Yhat)
v=-psi.getGrad().dot(Xhat)

U=psi.getFaceGrad()

from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm

eq=ImplicitDiffusionTerm(0.00001) \
	 + u * u.getGrad().dot(Xhat) \
	 + v * u.getGrad().dot(Yhat) \
	 - 2 * u.getGrad().dot(Yhat).getGrad().dot(Yhat)
from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition

from fipy.boundaryConditions.fixedValue import FixedValue
BC = (NthOrderBoundaryCondition(faces=mesh.getFacesBottom(), value=0,order=1),
      FixedValue(faces=mesh.getFacesBottom(), value=0),
      NthOrderBoundaryCondition(faces=mesh.getFacesLeft(), value=0,order=1),
	
      NthOrderBoundaryCondition(faces=mesh.getFacesTop(), value=1,order=1),
      NthOrderBoundaryCondition(faces=mesh.getFacesTop(), value=0,order=2))

from fipy.boundaryConditions.fixedFlux import FixedFlux


BC = (FixedFlux(faces=mesh.getFacesBottom(), value=0),
      FixedValue(faces=mesh.getFacesBottom(), value=0),
      FixedValue(faces=mesh.getFacesTop(), value=1),
      FixedFlux(faces=mesh.getFacesTop(), value=1),
      FixedFlux(faces=mesh.getFacesLeft(), value=0))



from fipy import viewers
viewer = viewers.make(vars=(psi,u,v,U))

eq.solve(var=psi,boundaryConditions=BC)
viewer.plot()

def hit_continue(Prompt='Hit any key to continue'):
	raw_input(Prompt)

hit_continue()
