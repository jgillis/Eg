from casadi import *
import numpy as n

(x,dx,theta,dtheta)=ssym("[x dx theta dtheta]")

q  = vertcat([x,theta])
dq = vertcat([dx,dtheta])
v = vertcat([q,dq])
ddq = ssym("[ddx,ddtheta]")

L=SX("L")
I=SX("I")
m=SX("m")
M=SX("M")
F=SX("F")
c=SX("c")
C=SX("C")

par = vertcat([F,L,I,m,M,L,c,C])

p = SXMatrix([x,0]) + L*SXMatrix([sin(theta),cos(theta)])
vp = mul(jacobian(p,q),dq)

g = 9.81

E = L/2*cos(theta)*g
T = 1.0/2*m*(dx)**2 + I*dtheta**2/2 + 1.0/2*M*mul(vp.T,vp)

L = E + T

W = jacobian(L,dq)

eq = SXMatrix(2,1,0)

forcing = vertcat([F-c*dx,-C*dtheta])

print W

# implicit form of lagrange
for i in range(2):
  eq[i] = mul(jacobian(W[i],q),dq) + mul(jacobian(W[i],dq),ddq) - jacobian(L,q[i]) - forcing[i]

  
A = SXMatrix(2,2)
b = SXMatrix(2,1,0)
  
for i in range(2):
  A[i,:] = jacobian(W[i],dq)
  b[i]   = - mul(jacobian(A[i],q),dq) + jacobian(L,q[i]) + forcing[i]

print eq

## Start intermezzo - linearized form
rhs = vertcat([dq,mul(inv(A),b)])

rhs_linear=mtaylor(rhs,[theta,dtheta],[0,0],1)
print rhs_linear

J = jacobian(rhs_linear,v)

print J
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
S.input(1).set([0,1,1,1,1,1,0.1,0.1])

S.evaluate()

r = S.output()

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

m = SXFunction([vertcat([v,T])],[1])
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
ms.setOption("nlp_solver_options",{"max_iter": 20})
ms.init()

from numpy import inf

ms.input(OCP_LBX).setAll(-inf)
ms.input(OCP_UBX).setAll(inf)
ms.input(OCP_LBX)[0,:] = -5
ms.input(OCP_LBX)[1,:] = -20
ms.input(OCP_UBX)[0,:] = 5
ms.input(OCP_UBX)[1,:] = 20

ms.input(OCP_LBX)[:,0] = ms.input(OCP_UBX)[:,0] = DMatrix([0,0,0,0,0])
ms.input(OCP_LBX)[:-1,-1] = ms.input(OCP_UBX)[:-1,-1] = DMatrix([-1,0,0,0])

ms.input(OCP_X_INIT).setAll(0)

ms.input(OCP_LBU).setAll(-1)
ms.input(OCP_UBU).setAll(1)
ms.input(OCP_U_INIT).setAll(0)

ms.input(OCP_LBP).set([0.5])
ms.input(OCP_UBP).set([1])

print ms.input(OCP_LBX)
print ms.input(OCP_UBX)
ms.solve()

print "t=", ms.output(OCP_P_OPT)

print ms.output(OCP_X_OPT)
print ms.output(OCP_U_OPT)

from pylab import *

x = ms.output(OCP_X_OPT)[0,:].toArray().squeeze()
theta = ms.output(OCP_X_OPT)[1,:].toArray().squeeze()
dx = ms.output(OCP_X_OPT)[2,:].toArray().squeeze()
dtheta = ms.output(OCP_X_OPT)[3,:].toArray().squeeze()
t = (ms.input(OCP_T)*ms.output(OCP_P_OPT)).toArray().squeeze()
print "t=", ms.output(OCP_P_OPT)

u = ms.output(OCP_U_OPT)[0,:].toArray().squeeze()
plot(t,x,'r',t,theta,'g',t,dx,'b',t,dtheta,'k')
plot(t[:-1],u)
print u
#show()

T = float(ms.output(OCP_P_OPT))
print T

results = []
for i in range(ns):
  tsim = n.linspace(float(ms.input(OCP_T)[i])*T,float(ms.input(OCP_T)[i+1])*T,100)
  sim = Simulator(I_original,tsim)
  sim.init()
  sim.input(INTEGRATOR_X0).set(ms.output(OCP_X_OPT)[:-1,i])
  sim.input(INTEGRATOR_P).set(ms.output(OCP_U_OPT)[:,i])
  sim.evaluate()
  results.append(array(sim.output()))


results = n.vstack(results)

plot(results[:,0])

show()

#plot(t,r[:,0],'r',t,r[:,1],'g')
