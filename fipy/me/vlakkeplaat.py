#!/usr/bin/env python
L = 1.0
N = 10
dL = L / N
viscosity = 1.
pressureRelaxation = 0.2
velocityRelaxation = 0.5

sweeps=10

from fipy.meshes.grid2D import Grid2D
mesh = Grid2D(nx=N, ny=N, dx=dL, dy=dL)

from fipy.variables.cellVariable import CellVariable

u = CellVariable(mesh=mesh, name='X velocity')
v = CellVariable(mesh=mesh, name='Y velocity')

from fipy.variables.vectorFaceVariable import VectorFaceVariable
velocity = VectorFaceVariable(mesh=mesh)

Xhat=(1,0)
Yhat=(0,1)
velocity=u.getFaceGrad().dot(Xhat) * Xhat + v.getFaceGrad().dot(Yhat) * Yhat

from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm

diffTerm=ImplicitDiffusionTerm(coeff=viscosity)

from fipy.terms.exponentialConvectionTerm import ExponentialConvectionTerm

convTerm = ExponentialConvectionTerm(coeff = velocity,diffusionTerm = diffTerm)

eq1= diffTerm - u * u.getGrad().dot(Xhat) - v * u.getGrad().dot(Yhat)
eq2= velocity.getDivergence()
#eq2= u.getGrad() + v.getGrad()
eq2=u.getFaceGrad().dot(Xhat) * Xhat + v.getFaceGrad().dot(Yhat) * Yhat
eq2= u.getGrad().dot(Xhat) + v.getGrad().dot(Yhat)
eq2= ImplicitDiffusionTerm(1) + u.getGrad().dot(Xhat) + v.getGrad().dot(Yhat)

from fipy.boundaryConditions.fixedValue import FixedValue
BCu = (FixedValue(faces=mesh.getFacesBottom(), value=0),
        FixedValue(faces=mesh.getFacesLeft(), value=1))

BCv = (FixedValue(faces=mesh.getFacesBottom(), value=0),
        FixedValue(faces=mesh.getFacesLeft(), value=0))

from fipy import viewers
viewer = viewers.make(vars=(u,v),limits={'datamin': 0., 'datamax': 1.})

def hit_continue(Prompt='Hit any key to continue'):
	raw_input(Prompt)

steps=10

for step in range(steps):
	eq1.solve(var=u,boundaryConditions=BCu)
	eq2.solve(var=v,boundaryConditions=BCv)
	viewer.plot()
	hit_continue()
	
