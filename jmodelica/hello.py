from pylab import *
from numpy import *
from jmodelica.jmi import *

jmu_name=compile_jmu("world","helloworld.mop")
dyn = JMUModel(jmu_name)
opts=dyn.simulate_options()
res = dyn.simulate(final_time=1)

plot(res['time'],res['x'])
plot(res['time'],res['y'])

show()
