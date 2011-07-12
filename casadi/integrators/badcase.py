from casadi import *
from integrators import *
from numpy import *
tend = 2*pi
t=SX("t")

x=SX("x") 
dx=SX("dx")

f=SXFunction({'NUM': ODE_NUM_IN, ODE_T: t, ODE_Y: [x,dx]},[[dx,-x]])
f.init()

ts = linspace(0,tend,10000)

F = RK4(f,ts)

F.init()
F.input(0).set([1,0])

F.evaluate()

J=F.jacobian()
J.init()
J.input(0).set([1,0])

J.evaluate()
