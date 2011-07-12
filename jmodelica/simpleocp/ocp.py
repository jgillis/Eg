#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
Why oh-why does this not work?
TypeError: 'list' object is not callable
"""


# Import library for path manipulations
import os.path

# Import numerical libraries
import numpy as N
from pylab import *

# Import the JModelica.org Python packages
from jmodelica.jmi import compile_jmu
from jmodelica.jmi import JMUModel


jmu_name = compile_jmu("opt","ocp.mop")

# Load the dynamic library and XML data
dyn = JMUModel(jmu_name)


dyn.initialize()

  
opts = dyn.optimize_options()
opts['n_e'] = 5	# number of elements
opts['n_cp'] = 3	# number of collocation points
opts['result_mode']='default'
opts['IPOPT_options']['max_iter']=100
opt = dyn.optimize(options=opts)

t=opt['time']

print t[-1]
plot(t,opt['sim.y'])

print opt['sim.y']
print opt['sim.a']
print opt['sim.b']
print opt['sim.x0']
print opt['sim.y0']
