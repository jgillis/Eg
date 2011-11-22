from pylab import *
from casadi import *
import numpy as n

q  = ssym("[x,theta]")
x,theta = q
dq = ssym("[dx,dtheta]")
dx,dtheta = dq
v = vertcat([q,dq])
ddq = ssym("[ddx,ddtheta]")

par = ssym("[F,L,I,m,M,L,c,C]")
F,L,I,m,M,L,c,C = par
par_ = [0,1,1,1,1,1,0.1,0.1]

g = 9.81

# The position vector of the pendulum's center of mass
p = SXMatrix([x,0]) + L/2*SXMatrix([sin(theta),cos(theta)])
# The velocity of the pendulum's center of mass
vp = mul(jacobian(p,q),dq)

# Potential energy of the system
E = L/2*cos(theta)*g
# Kinetic energy of the system
T = 1.0/2*m*(dx)**2 + I*dtheta**2/2 + 1.0/2*M*mul(vp.T,vp)

# The lagrangian
L = E + T


W = jacobian(L,dq)

eq = SXMatrix(2,1,0)

forcing = vertcat([F-c*dx,-C*dtheta])

# implicit form of lagrange
for i in range(2):
  eq[i] = mul(jacobian(W[i],q),dq) + mul(jacobian(W[i],dq),ddq) - jacobian(L,q[i]) - forcing[i]

  
A = SXMatrix(2,2)
b = SXMatrix(2,1,0)
  
for i in range(2):
  A[i,:] = jacobian(W[i],dq)
  b[i]   = - mul(jacobian(A[i],q),dq) + jacobian(L,q[i]) + forcing[i]

## Start intermezzo - linearized form
rhs = vertcat([dq,mul(inv(A),b)])

rhs_linear=mtaylor(rhs,[theta,dtheta],[0,0],1)

J = jacobian(rhs_linear,v)
## End intermezzo 

t = SX("t")

# ODE form
f = SXFunction([t,v,par,[]],[rhs])
f.init()
I = CVodesIntegrator(f)
I.init()

I_original = I

# Simulation
tl = n.linspace(0,10,1000)
S = Simulator(I,tl)
S.init()
S.input(0).set([0,0.2,0,0])
S.input(1).set(par_)

S.evaluate()

r = S.output()

period = 2*pi*sqrt(1/g)

print "period = ",  period

figure()
plot(tl,r[:,0],'r',label='x')
plot(tl,r[:,1],'g',label='theta')
legend()
title('Forward in time simulation')

# Time optimal control

T = SX("T") # The end time of intergation

rhs_s = substitute(rhs,par[1:],[1,1,1,1,1,0.01,0.01])
rhs_s = substitute(rhs_s,t,t*T)

# ODE form
f = SXFunction([t,v,F,[]],[rhs_s])
f.init()
I = CVodesIntegrator(f)
I.init()
I_original = I

# ODE form
f = SXFunction([t,vertcat([v,T]),vertcat([F,T]),[]],[vertcat([rhs_s*T,0])])
f.init()
I = CVodesIntegrator(f)
I.init()

m = SXFunction([vertcat([v,T])],[T])
m.init()


#c = SXFunction([t,vertcat([v,d]),F,[]],[v])
#c.init()

ns = 10

ms = MultipleShooting(I,m,FX())
ms.setOption("final_time",1)
ms.setOption("parallelization","expand")
ms.setOption("number_of_grid_points",ns)
ms.setOption("number_of_parameters",1)
ms.setOption("nlp_solver",IpoptSolver)
ms.setOption("nlp_solver_options",{"max_iter": 200})
ms.init()

from numpy import inf

ms.input(OCP_LBX).setAll(-inf)
ms.input(OCP_UBX).setAll(inf)
ms.input(OCP_LBX)[0,:] = -5
ms.input(OCP_LBX)[1,:] = -20
ms.input(OCP_UBX)[0,:] = 5
ms.input(OCP_UBX)[1,:] = 20

ms.input(OCP_LBX)[:-1,0] = ms.input(OCP_UBX)[:-1,0] = DMatrix([0,0,0,0])
ms.input(OCP_LBX)[:-1,-1] = ms.input(OCP_UBX)[:-1,-1] = DMatrix([-1,0,0,0])

ms.input(OCP_X_INIT).setAll(0)

ms.input(OCP_LBU).setAll(-10)
ms.input(OCP_UBU).setAll(10)
ms.input(OCP_U_INIT).setAll(0)

ms.input(OCP_LBP).set([0.99])
ms.input(OCP_UBP).set([1])

print ms.input(OCP_LBX)
print ms.input(OCP_UBX)
ms.solve()

print "t=", ms.output(OCP_P_OPT)

print ms.output(OCP_X_OPT)
print ms.output(OCP_U_OPT)



x = ms.output(OCP_X_OPT)[0,:].toArray().squeeze()
theta = ms.output(OCP_X_OPT)[1,:].toArray().squeeze()
dx = ms.output(OCP_X_OPT)[2,:].toArray().squeeze()
dtheta = ms.output(OCP_X_OPT)[3,:].toArray().squeeze()
t = (ms.input(OCP_T)*ms.output(OCP_P_OPT)).toArray().squeeze()
print "t=", ms.output(OCP_P_OPT)

u = ms.output(OCP_U_OPT)[0,:].toArray().squeeze()
figure()
hold(True)
plot(t,x,'r',label='x')
plot(t,theta,'g',label='theta')
plot(t,dx,'b',label='dx')
plot(t,dtheta,'k',label='dtheta')
plot(t[:-1],u,label='u')
title('Multiple Shooting solution -raw')
legend()
print u
#show()

T = float(ms.output(OCP_P_OPT))
print T

results = []
t = []
for i in range(ns):
  tsim = n.linspace(float(ms.input(OCP_T)[i])*T,float(ms.input(OCP_T)[i+1])*T,100)
  sim = Simulator(I_original,tsim)
  sim.init()
  sim.input(INTEGRATOR_X0).set(ms.output(OCP_X_OPT)[:-1,i])
  sim.input(INTEGRATOR_P).set(ms.output(OCP_U_OPT)[:,i])
  sim.evaluate()
  results.append(array(sim.output()))
  t.append(tsim)


results = n.vstack(results)
t = n.hstack(t)

print t.shape
print results.shape

figure()
plot(t,results[:,0])
title('Multiple Shooting solution - simulated')
show()

#plot(t,r[:,0],'r',t,r[:,1],'g')
