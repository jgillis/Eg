from greetcpp import *
import sys

import tempfile
import os

def greet():
  greeting()
  sys.stdout.write('greeting on stdout python\n')
  sys.stderr.write('greeting on stderr python\n')


print "*" * 10, "sys.stdout"

sys.stdout = open("stdout.txt", "w") 

greet()

sys.stdout = sys.__stdout__

print "*" * 10, "std::cout"

t = tempfile.mktemp()

captureCout(t)

greet()

releaseCout()

print "*" * 10, "result"

print open(t,'r').read()

os.remove(t)

t = tempfile.mktemp()

captureStdout(t)

greet()

releaseStdout()

print "*" * 10, "result"

print open(t,'r').read()

os.remove(t)

