import time

f = file('helloworld.txt','w')
while 1:
  time.sleep(1)
  f.write("Hello world @ %f!\n" % time.time())
  f.flush()
  
