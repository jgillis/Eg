from casadi import *
from integrators import *
from numpy import *
from pylab import *
from timeit import timeit
import pickle

results = pickle.load(file('integratortest.dat','r'))

figure()
for r in results:
  loglog(r['N'],r['y'],label=r['name'])
xlabel('Number of steps')
ylabel('Absolute error')
legend()
grid(True)
title('Performance, number of steps')

figure()
for r in results:
  loglog(r['N'],r['te'],label=r['name'])
xlabel('Number of steps')
ylabel('Execution time [s]')
legend()
grid(True)
title('Execution time, evaluation')

figure()
for r in results:
  loglog(r['N'],r['tj'],label=r['name'])
xlabel('Number of steps')
ylabel('Execution time [s]')
legend()
grid(True)
title('Execution time, jacobian')

figure()
for r in results:
  loglog(r['N'],r['tc'],label=r['name'])
xlabel('Number of steps')
ylabel('Execution time [s]')
legend()
grid(True)
title('Construction time')

show()


