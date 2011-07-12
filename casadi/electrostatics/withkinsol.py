from numpy import *
from casadi import *
import casadi as c
import numpy
from pylab import *
from numpy.linalg import solve, norm

N=2

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

x_ = array(numpy.linspace(0,2*pi*(1-1.0/N),N),ndmin=2).T

f = SXFunction([phi],[F])
f.init()
f.input().set(x_)
f.evaluate()
  
print f.output()

solver=KinsolSolver(f)
solver.init()
solver.setOption("constraints", [0]*N)
solver.output().set(x_)
print solver.output()
solver.solve()
print solver.output()




