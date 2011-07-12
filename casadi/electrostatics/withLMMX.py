from numpy import *
from casadi import *
import casadi as c
import numpy
from pylab import *
from numpy.linalg import solve, norm

from time import time

# N charged particles are free to find a position along the curve
N=100

print "Construction of the expression tree"
t0=time()


# The angle used in the parametrisation of the curve serves as decision variables
phi = MX("phi",N,1)

m=3
n1=3
n2=14 # n2/n3: must be integer and even
n3=2 # n2/n3: must be integer and even

# Our curve parametrisation
superformula = lambda phi: ((cos(m*phi/4))**n2+(sin(m*phi/4))**n3)**(-1/n1)

# A dummy scalar
a=SX("a")

r = superformula(a)

# A symbolic function that provides the tangent to the curve
dfp = SXFunction([a],[vertcat([jacobian(r*cos(a),a),jacobian(r*sin(a),a)])])
dfp.init()

r = superformula(phi)
x = r*cos(phi)
y = r*sin(phi)

# Evaluating the tangent for every point
t=horzcat([dfp.call([phi[i]])[0] for i in range(phi.size())]).T
tx=t[:,0]
ty=t[:,1]
n=sqrt(t[:,0]**2+t[:,1]**2)
tx/=n
ty/=n

# Evaluating the tangent for every point
dx = repmat(x,1,N)-repmat(x.T,N,1)
dy = repmat(y,1,N)-repmat(y.T,N,1)

# distance^2 matrix
N_ = dx**2+dy**2
# Denominator of Coulomb force
D = N_**(-3/2.0)

# k-indices to get diagonal elements
diagonal_k = list(getNZDense(sp_diag(N)))

# No self-interaction
D[diagonal_k]=MX(0)

# Summing all force contributions
n = MX.ones(N,1)
Fx_ = c.prod((D*dx),n)
Fy_ = c.prod((D*dy),n)

# Projecting the forces on the local tangents
F = Fx_*tx+Fy_*ty

f = MXFunction([phi],[F])
f.init()
J = Jacobian(f)
J.setOption("ad_mode","forward")
J.init()

print "duration: %f [s]" % (time()-t0)

# initial guess: evenly spaced
x_ = array(numpy.linspace(0,2*pi*(1-1.0/N),N),ndmin=2).T  # decision variables
J_ = zeros((N,N))  # Jacobian matrix
f_ = zeros((N,1))  # Function evaluation matrix
fn_ = zeros((N,1)) # Backup function evaluation matrix

f.input().set(x_)
J.input().set(x_)
f.evaluate()
J.evaluate()
f.output().get(f_)
J.output().get(J_)
  
c=0
# The Levenberg-Marquardt parameter, initial value
lambd = 0.001

# variables to hold timings
t_solve = t_f = t_J = 0

# Main loop of LM
while norm(f_.T) > 1e-9:
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
  
  # Did the residual shrink after taking a temptative step?
  if norm(fn_.T) < norm(f_.T):
    # If yes, make the step definitive
    x_+=dx_
    # Decrease the Levenberg-Marquardt parameter
    lambd=max(lambd/10,1e-16)
    # Get the function value and jacobian
    J.input().set(x_)
    t0_J=time()
    J.evaluate()
    t_J+=(time()-t0_J)
    J.output().get(J_)
    f_[:]=fn_[:]  # f_=fn_ would share the same memory
  else:
    # If no, forget the step and increase Levenberg-Marquardt parameter
    lambd=min(lambd*10,1e10)
  c+=1
  print norm(f_.T), lambd

print "%d steps for convergence" % c
print "duration linear solve: %f [s]" % t_solve
print "duration f: %f [s]" % t_f
print "duration J: %f [s]" % t_J

# Post-processing: make fancy plots

f = MXFunction([phi],[x,y,tx,ty,Fx_,Fy_])
f.init()
f.input().set(x_)
f.evaluate()

plot(f.output(0),f.output(1),'o')
quiver(f.output(0),f.output(1),f.output(2),f.output(3))
quiver(f.output(0),f.output(1),f.output(4),f.output(5),color='r')
show()



