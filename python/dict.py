print "With compare"

class Foo:
   def __init__(self,a):
     self.a = a
   def __repr__(self):
     return "(hash: %d, id: %d)" % (self.a , id(self))
   def __hash__(self):
     return self.a
   def __cmp__(self,a):
     return cmp(hash(self),hash(a))
   
     
a = Foo(1)
b = Foo(1)
c = Foo(2)

print a, hash(a), id(a)
print b, hash(b), id(b)
print c, hash(c), id(c)

d = {}

d[a] = 1
d[b] = 2
d[c] = 3


s = set([])
s.add(a)
s.add(b)
s.add(c)

print d
print s

print "Without compare"

class Foo:
   def __init__(self,a):
     self.a = a
   def __repr__(self):
     return "(hash: %d, id: %d)" % (self.a , id(self))
   def __hash__(self):
     return self.a
   
     
a = Foo(1)
b = Foo(1)
c = Foo(2)

print a, hash(a), id(a)
print b, hash(b), id(b)
print c, hash(c), id(c)

d = {}

d[a] = 1
d[b] = 2
d[c] = 3


s = set([])
s.add(a)
s.add(b)
s.add(c)

print d
print s
