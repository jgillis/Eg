from StringIO import StringIO
import sys

class TeeString(StringIO):
  def __init__(self,stream):
    StringIO.__init__(self)
    self.stream = stream
        
  def write(self, data):
    StringIO.write(self,data)
    self.stream.write(data)

class Stdout():
    def __init__(self):
        self.stream = TeeString(sys.stdout)

    def __enter__(self):
        sys.stdout = self.stream
        return self.stream

    def __exit__(self, type, value, traceback):
        sys.stdout = self.stream.stream
        
class Stderr():
    def __init__(self):
        self.stream = TeeString(sys.stderr)

    def __enter__(self):
        sys.stderr = self.stream
        return self.stream

    def __exit__(self, type, value, traceback):
        sys.stderr = self.stream.stream
        


print "foo"
        
with Stdout() as stdout:
   print "waaw"

print "bar"

print stdout.getvalue()


