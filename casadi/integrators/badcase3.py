from casadi import *
from integrators import *
from numpy import *
tend = 2*pi
t=SX("t")

n=2

x  = symbolic("x",n,1) 

A = symbolic("A",n,n)
dx = casadi.dot(A,x)

f=SXFunction({'NUM': ODE_NUM_IN, ODE_T: t, ODE_Y: x, ODE_P: A},[dx])
f.init()

ts = linspace(0,tend,3)

F = Euler(f,ts)

F.init()
print countNodes(F.outputSX())
F.generateCode('F.txt')
print F.outputSX()
F.evaluate()

J=F.jacobian()
J.init()
print countNodes(J.outputSX())
J.generateCode('J.txt')
print J.outputSX()
J.evaluate()
