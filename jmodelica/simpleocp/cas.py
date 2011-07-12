# JModelica
from jmodelica.jmi import compile_jmu
from jmodelica.jmi import JMUModel
import jmodelica

import zipfile
import os
from casadi import *
from numpy import *


# We're gonna parse modelica via jmodelica
try:
  # Try the old Jmodelica syntax
  pass
  #jmu_name = compile_jmu("opt", "ocp.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_xml':False})
except jmodelica.compiler.UnknownOptionError:
  # Try the new jmodelica syntax
  pass
  #jmu_name = compile_jmu("opt", "ocp.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_me_xml':False})
  
sfile = zipfile.ZipFile('opt.jmu','r')
mfile = sfile.extract('modelDescription.xml','.')
os.rename('modelDescription.xml','opt.xml')

# read the xml that jmodelica generated
parser = FMIParser('opt.xml')

# Obtain the symbolic representation of the OCP
ocp = parser.parse()

ocp.sortType()

print ocp

params = ocp.explicit_lhs_
values = ocp.explicit_rhs_

states = [i.var() for i in ocp.xd_]
print "states:", states

dstates = [i.der() for i in ocp.xd_]
print "derivatives:", dstates

print "These are trivial equations: " , zip(params,values)


print "These are the implicit dynamic equations: ", ocp.dynamic_eq_

t=SX("t")


f = SXFunction({'NUM': DAE_NUM_IN,
                 DAE_T: t,
                 DAE_Y: states,
                 DAE_YDOT: dstates,
                 DAE_P: params},
                 [ocp.dynamic_eq_])
                        
f.init()

integr = IdasIntegrator(f)
integr.init()

integr.input(INTEGRATOR_T0).set(ocp.t0)
integr.input(INTEGRATOR_TF).set(ocp.tf)
integr.input(INTEGRATOR_P).set([i.getValue() for i in values])
integr.input(INTEGRATOR_X0).set([0,1])

integr.evaluate()
print integr.output()

