from pylab import *
from numpy import *
from jmodelica.jmi import *


jmu_name=compile_jmu("DynabeePower","Dynabee.mop")
dyn = JMUModel(jmu_name)


opts=dyn.simulate_options()
opts['solver'] = 'CVode'
opts['CVode_options']['rtol'] = 1e-7;
res = dyn.simulate(final_time=1)
