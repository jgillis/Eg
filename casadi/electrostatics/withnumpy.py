from numpy import *
from pylab import *

N=5

# n2/n3: must be integer and even

a=1.0
b=1.0
m=3
n1=3
n2=14
n3=2

superformula = lambda phi: ((cos(m*phi/4))**n2+(sin(m*phi/4))**n3)**(-1/n1)
superformula = lambda phi: 1

a=array(linspace(0,2*pi*(1-1.0/N),N),ndmin=2).T

r = superformula(a)
x = r*cos(a)
y = r*sin(a)


dx = (tile(x,(1,N))-tile(x.T,(N,1)))
dy = (tile(y,(1,N))-tile(y.T,(N,1)))

print "dx=",dx
print "dy=",dy

N_ = dx**2+dy**2

print "N=",N_

D = N_**(-3/2.0)

print "D=",D

for i in range(N):
  D[i,i]=0

print "D=",D

print "D*dx=",D*dx
print "D*dy=",D*dy

Fx_ = sum(D*dx,1)
Fy_ = sum(D*dy,1)

print Fx_, Fy_

figure(figsize=(8,8))
plot(x,y,'o')
for i in range(N):
  quiver(x,y,(D*dx)[:,i],(D*dy)[:,i],angles='xy',scale=4)
quiver(x,y,Fx_,Fy_,angles='xy',scale=4,color='r')
axis([-2,4,-2,4])
show()






