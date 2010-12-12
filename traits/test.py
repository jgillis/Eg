
class MyInstance(Instance):
	def __init__(self,myclass, foo):
		Instance.__init__(self,myclass)
		print "random stuff"
		self.garble=42
		
	def get(self,object,name):
		return (142,object.__dict__[name])
		
	def set(self,object,name,value):
		print "Set to ", value
		#object.__dict__[name]=value
		self.set_value(object,name,value)
		
class Baz():
	pass
	
class Foo(HasTraits):
	f=Instance(Baz)
	g=MyInstance(Baz, 123)
	
	def __init__(self):
		self.f=Baz()
		self.g=Baz()
	
f=Foo();

#f.f=Baz()
#f.g=Baz()

print Dict();