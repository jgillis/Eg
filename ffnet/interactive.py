execfile("/home/JG/eg/matplotlib/interactive.py")

from ffnet import ffnet

n=100

p=plot(range(2),visible=False)
p[0].visible=False
draw()

d=DragletInteractor()

def do(inp):
	conec = mlgraph((1,2,2,1))
	net = ffnet(conec)
	input=inp[0]
	target=inp[1]
	#net.train_genetic(input, target, individuals=20, generations=500)
	print "TRAINING NETWORK..."
	net.train_tnc(input, target, maxfun = 5000)
	x=mgrid[0:1:n*1j]
	y=map(lambda x: net([x])[0],x)
	return vstack((x,y))
	
d.setEvent('changed',d.lineDrawAgent(do))


# conclusies 2:10 inputs, (1,2:4,1) zeer onstabiel gedrag.
# genetic gedeelte is onschuldig

