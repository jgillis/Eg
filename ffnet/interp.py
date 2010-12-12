########################################################################
##  Copyright (C) 2006 by Marek Wojciechowski
##  <mwojc@p.lodz.pl>
##
##  Distributed under the terms of the GNU General Public License (GPL)
##  http://www.gnu.org/copyleft/gpl.html
########################################################################

### Sine training example for ffnet ###

from ffnet import ffnet
from math import pi, sin, cos
from pylab import *

# Let's define network connectivity by hand and then create network.
conec = [(1, 2), (1, 3), (1, 4), (1, 5), (2, 6), (3, 6), (4, 6), (5, 6), \
         (0, 2), (0, 3), (0, 4), (0, 5), (0, 6)]
# Note 1: Biases in ffnet are handled as the connections
#         from special node numbered 0. Input nodes cannot be biased.
# Note 2: Node numbering and order of links in conec is meaningless,
#         but the connections have to be from source to target.
# Note 3: The same connectivity can be obtained using mlgraph function
#         provided with ffnet (layered architecture (1,4,1)).
net = ffnet(conec)

# Generate training data (sine values for x from 0 to 2*pi)

fixn=10
fixe=0.1

dinput=0.1/(fixn-1.0)
input=linspace(0,fixe,fixn)
target=input

a=linspace(0.4,0.8,20)
b=a*(1.2)-0.2

input = hstack((input,a))
target = hstack((target,b))

input = hstack((input,linspace(1-fixe,1,fixn)))
target = hstack((target,linspace(1-fixe,1,fixn)))


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

points = 128
xaxis = linspace(0,1,points)
sine = [ sin(x) for x in xaxis ]
cosine = [ cos(x) for x in xaxis ]
netsine = [ net([x])[0] for x in xaxis]
netcosine = [ net.derivative([x])[0][0] for x in xaxis ]
    
subplot(211)
plot(input, target, 'b-o', xaxis, netsine, 'k-')
grid(True)
title('Outputs of trained network.')
    
subplot(212)
plot(input[0:-1], diff(target)/dinput, 'b-o', xaxis, netcosine, 'k-')
grid(True)
show()
