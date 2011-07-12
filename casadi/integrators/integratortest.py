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
  nce=zeros(N.shape)
  ncj=zeros(N.shape)
  
  for i in range(len(N)):
    Nsteps = N[i]
    print "Nsteps = ", Nsteps

    ts = linspace(0,tend,Nsteps)

    F = method(f,ts)
    tc[i]=timeit(stmt="method(f,ts)", setup="from __main__ import f,ts,method",number=10)/10
    
    F.init()
    F.input(0).set([1,0])
    nce[i] = countNodes(F.outputSX())
    te[i]=timeit(stmt="F.evaluate()", setup="from __main__ import F",number=100)/100
    
    J=F.jacobian() # symbolic jacobian
    J.init()
    ncj[i] = countNodes(J.outputSX())
    J.input(0).set([1,0])

    tj[i]=timeit(stmt="J.evaluate()", setup="from __main__ import J",number=100)/100
    
    y[i]=abs(F.output()[0]-1)
  
  results.append({'method':method,'N':N,'x':x,'y':y,'te':te,'tj':tj,'name':name,'tc':tc,'nce':nce,'ncj':ncj})
  
pickle.dump(results,file('integratortest.dat','w'))
