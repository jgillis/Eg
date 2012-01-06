import time
import os

print os.environ

f = file('helloworld.txt','w')
for i in range(5):
  time.sleep(1)
  f.write("Hello world @ %f!\n" % time.time())
  f.flush()


print "foo"

raise Exception('bar')
  
