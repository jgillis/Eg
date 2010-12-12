class foo:
	bar=3.14
	def __init__(self):
		self.kreeft=3
		print "Huma"
class patat(foo):
	def __init__(self):
		foo.__init__(self)
		print self.bar
		print self.kreeft
		
		
p=patat()