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

s=IpoptSolver(f)
s.setOption('max_iter',1000)
s.setOption('derivative_test','first-order')
s.init()

print " duration: %f [s]" % (t0-time())
print
t=time()


s.input(NLP_X_INIT).set(numpy.linspace(0,2*pi*(1-1.0/N),N))
s.input(NLP_LBX).set([-100]*N)
s.input(NLP_UBX).set([100]*N)
s.solve()

print s.output()

f = SXFunction([phi],[x,y,tx,ty,Fx_,Fy_])
f.init()
f.input().set(s.output())
f.evaluate()
print f.output(),f.output(0).toArray()

plot(f.output(0),f.output(1),'o')
quiver(f.output(0),f.output(1),f.output(2),f.output(3))
quiver(f.output(0),f.output(1),f.output(4),f.output(5),color='r')
show()





