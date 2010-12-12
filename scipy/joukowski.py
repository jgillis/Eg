close('all')
Vi=1.0;

alphaDeg=15;
alpha=(alphaDeg+0.0)/180*pi;

a=1.0;

cx=-0.15;
cy=0.1;

c=cx+cy*1j;
s=a;

r=abs(c-s);
beta=-angle(c-s);
C=4*pi*Vi*r*sin(beta+alpha);
omega=-C/2/pi/r**2;


contz=r*exp(1j*mgrid[0:2*pi:100j])+c

n=400;
dx=20.0/n;

x,y=mgrid[-10:10:n*1j,-10:10:n*1j];
z=x+y*1j;


zc=z*exp(-1j*alpha)-c;
F=Vi*(zc+r**2/zc)+1j*r**2*omega*log(zc);

V=mgrid[-4:4:81j]

z[abs(z-c)<r]=NaN;


figure()
contour(real(z),imag(z),imag(F),V)
fill(real(contz),imag(contz),'b')
plot([a],[0],'or')
axis('equal')

J=z+a**2/z;
contJ=contz+a**2/contz;

figure()
contour(real(J),imag(J),imag(F),V)
fill(real(contJ),imag(contJ),'b',zorder=1)
plot([2*a],[0],'or')
axis('equal')


figure()
v=diff(real(F),axis=-1)[1:,0:]/dx;
u=diff(real(F),axis=-2)[0:,1:]/dx;

zs=(z[0:-1,0:-1]+z[1:,1:])/2;
zs[abs(z-c)<r]=NaN;

d=10;

quiver(real(zs[0::d,0::d]),imag(zs[0::d,0::d]),u[0::d,0::d],v[0::d,0::d],scale=50)
axis('equal')

figure()
Vinf=Vi*exp(1j*alpha);
ui=real(Vinf);
vi=imag(Vinf);


v2=v-vi;
u2=u-ui;

quiver(real(zs[0::d,0::d]),imag(zs[0::d,0::d]),u2[0::d,0::d],v2[0::d,0::d],scale=10)
axis('equal')

figure()
from scipy.interpolate import *
#interp2d(real(J),imag(J),real(F))
