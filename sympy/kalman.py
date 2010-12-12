from sympy import *

N=5

dt=Symbol('dt')
t=Symbol('t')
epsilon=Symbol('epsilon')
q=map(Symbol,["x%d" % i for i in range(N)])

Q=Matrix(q)

M=Matrix(N,N,lambda i,j:1 if j==i+1 else 0)
P=Matrix(N,1,lambda i,j:1 if i==N-1 else 0)

dQ=Q.applyfunc(lambda x: diff(x(t),t))

print "Equations of motion":
e=dQ-(M*Q+P*epsilon)

print "Discretization: s-> (1-1/z)/dT
       q
       q[t-1]
"

qc=map(Symbol,["x%d'" % i for i in range(N)])
Qc=Matrix(qc)

print (Q-Qc)/dt-(M*Q+P*epsilon)
print (eye(N)/dt-M)*Q-P*epsilon-Qc/dt

i=(eye(N)/dt-M).inv()
A=i/dt
B=i*P
print (eye(N)/dt-M)*Q-P*epsilon-Qc/dt


print "Different discretizations:"
(eye(N)/dt+M)
(eye(N)/dt-M).inv()
(eye(N)+M/2*dt)*(eye(N)-M/2*dt).inv()
