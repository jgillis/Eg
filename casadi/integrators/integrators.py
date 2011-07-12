from casadi import *
from numpy import *
import casadi
def ExplicitFixedStepIntegrator(f,times=None,a=None,b=None,c=None):
  """ a,b,c are the tableau coefficients
  
  If s is the number of stages, then we have:
  
  a: tril(s-1 x s-1)
  b,c: (sx1)
  
  times may be DVector or SXVector
  
  """
  
  if not(isinstance(times,DMatrix)):
    times = DMatrix(times)
    
  
  def toSX(a):
    return casadi.reshape(SXMatrix(a),a.shape[0],a.shape[1])
    
    
  times = toSX(times)
  a     = toSX(a)
  b     = toSX(b)
  c     = toSX(c)
  
  x_init = f.inputSX(ODE_Y)
  N = x_init.numel()
  p = f.inputSX(ODE_P)
  
  s=b.numel()
  assert(a.size1()==s-1)
  assert(a.size2()==s-1)
  assert(c.numel()==s)
  
  if s>1:
    for lhs,rhs in zip(c[1:,0],casadi.sum(a,1)):
     pass
     #assert(lhs==rhs)
    
  ks = SXMatrix(N,s)
  y = x_init
  
  for k in range(len(times)-1):
    t = times[k]
    h = times[k+1]-times[k]
    for i in range(s):
      if i>0:
        x = y + casadi.dot(ks[:,:i],a[i-1,:i].T)*h
      else:
        x = y
      ks[:,i] = f.eval({ODE_T: t+c[i,0]*h, ODE_Y: x, ODE_P: p})[0]
    y+= casadi.dot(ks,b)*h
    
  return SXFunction([x_init,p],[y])
  
def Euler(f,times=None):
  return ExplicitFixedStepIntegrator(f,times=times,a=DMatrix(),b=DMatrix([1]),c=DMatrix([0]))
  
  
def RK4(f,times=None):
  return ExplicitFixedStepIntegrator(f,times=times,
      a=DMatrix(array([[0.5,0,0],[0,0.5,0],[0,0,1]])),
      b=DMatrix([1.0/6,1.0/3,1.0/3,1.0/6]),
      c=DMatrix([0,0.5,0.5,1]))
      
def RK2(f,times=None):
  return ExplicitFixedStepIntegrator(f,times=times,
      a=DMatrix(array([[2.0/3]])),
      b=DMatrix([1.0/4,3.0/4]),
      c=DMatrix([0,2.0/3]))
      
def debug(f,times=None,s=1):
  return ExplicitFixedStepIntegrator(f,times=times,a=symbolic("a",s-1,s-1),b=symbolic("b",s,1),c=symbolic("c",s,1))
  

if __name__ == "__main__":
  t=SX("t")

  x=SX("x") 
  dx=SX("dx")

  f=SXFunction({'NUM': ODE_NUM_IN, ODE_T: t, ODE_Y: [x,dx]},[[dx,-x]])
  f.init()
  
  if 0:
  
    f=SXFunction({'NUM': ODE_NUM_IN, ODE_T: t, ODE_Y: [x]},[tan(x)+1])
    f.init()
    
    times = DMatrix([1,1.025])
    e = RK2(f,times)
    e.init()
    e.input(0).set([1])
    e.evaluate()
    print e.output()
      
    if 0:
      times = SXMatrix([0,SX("h")])
      e = debug(f,times,3)
    
  
  times = DMatrix(linspace(0,pi/2,5))
  #times = DMatrix([0,0.1,0.2])
  e = Euler(f,times)
  print e.outputSX()
  e.init()
  e.input(0).set([1,0])
  e.evaluate()
  print e.output()
  
  times = DMatrix(linspace(0,pi/2,100))
  e = RK4(f,times)
  e.init()
  e.input(0).set([1,0])
  e.evaluate()
  print e.output()
