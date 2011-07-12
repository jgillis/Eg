from numpy import *
import numpy
from pylab import *
from numpy.linalg import solve, norm
from theano import *
from theano.tensor import *

from time import time

N=100

print "Construction of the expression tree"
t0=time()

phi = vector()

# n2/n3: must be integer and even

a=1.0
b=1.0
m=3
n1=3
n2=14
n3=2

superformula = lambda phi: ((cos(m*phi/4))**n2+(sin(m*phi/4))**n3)**(-1/n1)

a=scalar()

r = superformula(a)
x = r*cos(a)
y = r*sin(a)
dfpx = function(inputs=[a],outputs=[grad(x,a)])
dfpy = function(inputs=[a],outputs=[grad(y,a)])

dfp.init()
n=sqrt(x**2+y**2)


r = superformula(phi)
x = r*cos(phi)
y = r*sin(phi)

n=sqrt(x**2+y**2)
t=horzcat([dfp.call([phi[i]])[0] for i in range(phi.size())]).T

tx=t[:,0]
ty=t[:,1]
n=sqrt(t[:,0]**2+t[:,1]**2)
tx/=n
ty/=n

dx = repmat(x,1,N)-repmat(x.T,N,1)
dy = repmat(y,1,N)-repmat(y.T,N,1)

N_ = dx**2+dy**2
D = N_**(-3/2.0)

diagonal_k = list(getNZDense(sp_diag(N)))

D[diagonal_k]=MX(0)

n = MX.ones(N,1)

Fx_ = c.prod((D*dx),n)
Fy_ = c.prod((D*dy),n)

F = Fx_*tx+Fy_*ty

f = MXFunction([phi],[F,tx])
f.init()
J = Jacobian(f)
J.setOption("ad_mode","forward")
J.init()

print "duration: %f [s]" % (time()-t0)


x_ = array(numpy.linspace(0,2*pi*(1-1.0/N),N),ndmin=2).T
J_ = zeros((N,N))
f_ = zeros((N,1))
fn_ = zeros((N,1))

I=eye(N)
res=Inf

f.input().set(x_)
J.input().set(x_)
J.fwdSeed().set([1]+[0]*(N-1))
f.evaluate(1,0)
J.evaluate()
f.output().get(f_)
J.output().get(J_)
  
c=0
lambd = 0.001

t_solve = 0
t_f =0
t_J =0

while norm(f_.T) > 1e-10:
  JJ=numpy.dot(J_.T,J_)
  JF=numpy.dot(J_.T,f_)
  t0_solve=time()
  dx_ = -solve(JJ+diag(diag(JJ))*lambd,JF)
  t_solve+=time()-t0_solve

  f.input().set(x_+dx_)
  t0_f=time()
  f.evaluate()
  t_f+=(time()-t0_f)
  f.output().get(fn_)
  
  if norm(fn_.T) < norm(f_.T):
    print "%f -> %f " % (norm(f_.T), norm(fn_.T))
    lambd=max(lambd/10,1e-16)
    x_+=dx_
    J.input().set(x_)
    t0_J=time()
    J.evaluate()
    t_J+=(time()-t0_J)
    J.output().get(J_)
    f_[:]=fn_[:]  # f_=fn_ would share the same memory
  else:
    print "%f <-> %f " % (norm(f_.T), norm(fn_.T))
    lambd=min(lambd*10,1e10)
    
  c+=1
  print norm(f_.T), lambd

print "%d steps for convergence" % c

print "duration linear solve: %f [s]" % t_solve
print "duration f: %f [s]" % t_f
print "duration J: %f [s]" % t_J

f = MXFunction([phi],[x,y,tx,ty,Fx_,Fy_])
f.init()
f.input().set(x_)
f.evaluate()

plot(f.output(0),f.output(1),'o')
quiver(f.output(0),f.output(1),f.output(2),f.output(3))
quiver(f.output(0),f.output(1),f.output(4),f.output(5),color='r')
show()



