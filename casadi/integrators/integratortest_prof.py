from casadi import *
from integrators import *
from numpy import *
from pylab import *
from timeit import timeit
import pickle

tend = 2*pi
t=SX("t")

x=SX("x") 
dx=SX("dx")

f=SXFunction({'NUM': ODE_NUM_IN, ODE_T: t, ODE_Y: [x,dx]},[[dx,-x]])
f.init()

results=[]

for method, name,N in [(Euler,"Euler",logspace(0.5,4,50)),(RK2,"RK2",logspace(0.5,3,50)),(RK4,"RK4",logspace(0.5,2,50))]:
  N = floor(N)
  x=tend/N
  y=zeros(N.shape)
  te=zeros(N.shape)
  tj=zeros(N.shape)
  tc=zeros(N.shape)
    
  for i in range(len(N)):
    Nsteps = N[i]
    print "Nsteps = ", Nsteps

    ts = linspace(0,tend,Nsteps)

    F = method(f,ts)

    F.init()
    F.input(0).set([1,0])
    
    F.evaluate()
    
    J=Jacobian(F)
    J.init()
    J.input(0).set([1,0])

    J.evaluate()
    
    y[i]=abs(F.output()[0]-1)
  
  results.append({'method':method,'N':N,'x':x,'y':y,'name':name})
  
