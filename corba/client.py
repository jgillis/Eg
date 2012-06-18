#!/usr/bin/env python
import sys
from omniORB import CORBA
import CosNaming, Example

# Initialise the ORB
orb = CORBA.ORB_init(sys.argv, CORBA.ORB_ID)

# Obtain a reference to the root naming context
obj         = orb.resolve_initial_references("NameService")
rootContext = obj._narrow(CosNaming.NamingContext)

if rootContext is None:
    print "Failed to narrow the root naming context"
    sys.exit(1)
    

context = rootContext    
# Recursively list all directories

def explore(context, bindingNames = []):
  if isinstance(context,CosNaming._objref_NamingContext):
    for c in context.list(100)[0]:
      for v in explore(context.resolve(c.binding_name),bindingNames + c.binding_name):
        yield v
      # python 3.3: yield from explore(context.resolve(c.binding_name),bindingNames + c.binding_name)
  else:
    yield bindingNames
    
print "All objects registered on the NameServer:"
for namehierarchy in explore(context):
  print "/".join(map(lambda x : x.id + "." + x.kind,namehierarchy))

# Resolve the name "test.my_context/ExampleEcho.Object"
name = [CosNaming.NameComponent("test", "my_context"),
        CosNaming.NameComponent("ExampleEcho", "Object")]
try:
    obj = rootContext.resolve(name)

except CosNaming.NamingContext.NotFound, ex:
    print "Name not found"
    sys.exit(1)

# Narrow the object to an Example::Echo
eo = obj._narrow(Example.Echo)

if (eo is None):
    print "Object reference is not an Example::Echo"
    sys.exit(1)

# Invoke the echoString operation
message = "Hello from Python"
result  = eo.echoString(message)

print "I said '%s'. The object said '%s'." % (message,result)
