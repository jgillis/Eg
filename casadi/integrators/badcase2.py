from casadi import *
from integrators import *
from numpy import *
tend = 2*pi
t=SX("t")

n=10

x  = symbolic("x",n,1) 

A = symbolic("A",n,n)
dx = casadi.dot(A,x)

f=SXFunction({'NUM': ODE_NUM_IN, ODE_T: t, ODE_Y: x, ODE_P: A},[dx])
f.init()

ts = linspace(0,tend,100)

F = RK4(f,ts)

F.init()
print countNodes(F.outputSX())

F.evaluate()

J=F.jacobian()
J.init()
print countNodes(J.outputSX())

J.evaluate()
