from  scipy.interpolate import bisplrep,bisplev
from random import uniform


def r():
	return uniform(-5,5)

def plt(n=25):
	x=[]
	y=[]
	qx=[]
	qy=[]
	for i in range(n):
		x.append(r())
		y.append(r())
		qx.append(sin(x[-1]))
		qy.append(cos(y[-1]))
	qxb=bisplrep(x,y,qx,s=0)
	qyb=bisplrep(x,y,qy,s=0)
	X=arange(-2,2,0.4)
	Y=arange(-2,2,0.4)
	cla()
	hold(True)
	quiver(x,y,qx,qy,pivot='tail',color='b')
	quiver2(X,Y,bisplev(X, Y,qxb),bisplev(X, Y,qyb),pivot='tail',color='r')
	hold(False)
