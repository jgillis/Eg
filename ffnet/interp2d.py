########################################################################
##  Copyright (C) 2006 by Marek Wojciechowski
##  <mwojc@p.lodz.pl>
##
##  Distributed under the terms of the GNU General Public License (GPL)
##  http://www.gnu.org/copyleft/gpl.html
########################################################################

### Sine training example for ffnet ###

from ffnet import ffnet, mlgraph
from math import pi, sin, cos
from pylab import *
from numpy import *
import pylab as p
import matplotlib.axes3d as p3

# Let's define network connectivity by hand and then create network.
conec = mlgraph((2,4,1))
# Note 1: Biases in ffnet are handled as the connections
#         from special node numbered 0. Input nodes cannot be biased.
# Note 2: Node numbering and order of links in conec is meaningless,
#         but the connections have to be from source to target.
# Note 3: The same connectivity can be obtained using mlgraph function
#         provided with ffnet (layered architecture (1,4,1)).
net = ffnet(conec)

# Generate training data (sine values for x from 0 to 2*pi)

x,y=mgrid[0:1:5j,0:1:5j]
z=exp(x)+sin(6*y)

fig=p.figure()
ax = p3.Axes3D(fig)

ax.scatter3D(ravel(x),ravel(y),ravel(z))
p.show()

input=vstack((ravel(x),ravel(y))).transpose()
target=ravel(z)

# Train network
#first find good starting point with genetic algorithm (not necessary, but may be helpful)
print "FINDING STARTING WEIGHTS WITH GENETIC ALGORITHM..."
net.train_genetic(input, target, individuals=20, generations=500)
#then train with scipy tnc optimizer
print "TRAINING NETWORK..."
net.train_tnc(input, target, maxfun = 5000, messages=1)

# Test network
print
print "TESTNG NETWORK..."
output, regression = net.test(input, target, iprint = 2)

x,y=mgrid[0:1:25j,0:1:25j]
m=vstack((ravel(x),ravel(y))).transpose()
out=array([ net([s])[0] for s in m])

ax.scatter3D(ravel(x),ravel(y),out,color='r')


p.show()

