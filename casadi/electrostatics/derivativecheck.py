from numpy import *
from casadi import *
import casadi as c
import numpy
from pylab import *

from time import time

N=100


print "Construction of the expression tree"
t0=time()

phi = symbolic("phi",N,1)

# n2/n3: must be integer and even

a=1.0
b=1.0
m=3
n1=3
n2=14
n3=2

superformula = lambda phi: ((cos(m*phi/4))**n2+(sin(m*phi/4))**n3)**(-1/n1)

a=SX("a")

r = superformula(a)
x = r*cos(a)
y = r*sin(a)
fp = SXFunction([a],[[x,y]])
fp.init()
dfp = fp.jacobian()
dfp.init()
n=sqrt(x**2+y**2)


r = superformula(phi)
x = r*cos(phi)
y = r*sin(phi)

n=sqrt(x**2+y**2)
t=SXMatrix(matrix([dfp.eval([phi[i]])[0] for i in range(phi.size())]))
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

D[diagonal_k]=SX(0)

Fx_ = c.sum(D*dx,1)
Fy_ = c.sum(D*dy,1)

F = Fx_*tx+Fy_*ty





f = SXFunction([phi],[norm_2(F)**2])
f.init()
f.input().set(numpy.linspace(0,2*pi*(1-1.0/N),N))
f.evaluate()
f_ = f.output()[0]

J_fd=[]
e=eye(N)
for i in range(N):
  step = 1e-6
  f.input().set(numpy.linspace(0,2*pi*(1-1.0/N),N)+e[:,i]*step)
  f.evaluate()
  J_fd.append((f.output()[0]-f_)/step)
J_fd_ = array(J_fd)

J_symbolic = f.jacobian()
J_symbolic.init()
J_symbolic.input().set(numpy.linspace(0,2*pi*(1-1.0/N),N))
J_symbolic.evaluate()

J_symbolic_= J_symbolic.output().toArray()


J_fwd = Jacobian(f)
J_fwd.setOption("ad_mode","forward")
J_fwd.init()
J_fwd.input().set(numpy.linspace(0,2*pi*(1-1.0/N),N))
J_fwd.evaluate()
J_fwd_= J_fwd.output().toArray()

J_adj = Jacobian(f)
J_adj.setOption("ad_mode","adjoint")
J_adj.init()
J_adj.input().set(numpy.linspace(0,2*pi*(1-1.0/N),N))
J_adj.evaluate()
J_adj_= J_adj.output().toArray()

print "symbolic= ", J_symbolic_
print "fwd= ", J_fwd_
print "adj= ", J_adj_


print "res fwd - symbolic:", max(abs((J_fwd_ - J_symbolic_)/J_symbolic_).ravel())
print "res adj - symbolic:", max(abs((J_adj_ - J_symbolic_)/J_symbolic_).ravel())
print "res fd - symbolic:", max(abs((J_fd_  - J_symbolic_)/J_symbolic_).ravel())
print "finite diff= ", (J_fd_  - J_symbolic_)/J_symbolic_
