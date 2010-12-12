from enthought.traits.api import TraitCompound, Any, BaseTraitHandler,Tuple,Either, Int, Expression,TraitHandler,TraitType , Str, Range, HasTraits, Instance, Dict, Float, Bool, on_trait_change, List, Trait,DelegatesTo





		
def TExpression(mytrait):
	if isinstance(mytrait,DelegatesTo):
		print "interesting"
		return TExpressionTraitDelegatesTo(mytrait.delegate,mytrait.prefix)
		print " or not"
	else:
		return TExpressionTrait(mytrait)

class Bar(HasTraits):
	foo=Float(15)
	
	
class ExpressionTraitListener(HasTraits):
	def Echanged(self, new):
		print self.name + "=" , new
		setattr(self,self.name,new)
		
	def changed(self, new):
		print self.name + "=" ,  new
		self.object.changed(self.name,new)
		
	def __init__(self,object,name):
		self.name = name
		self.object=object
		setattr(self,'_E_' + name + '_changed',self.Echanged)
		print "registering " +'_E_' + name + '_changed'
		print "registering " +'_' + name + '_changed'
		setattr(self,'_' +  name + '_changed',self.changed)
	def _width_changed(self):
		print "howdy\n\n"
		
class TExpressionTrait(TraitType):
	def validate(self,object,name,value):
		return value
		
	def __init__(self,mytrait):
		self.mytrait=mytrait
		TraitType.__init__(self)
	
class TExpressionTraitDelegatesTo(TraitType):
	def __init__(self,delegate,prefix):
		self.delegate=delegate
		self.prefix=prefix
		TraitType.__init__(self)

	def validate(self,object,name,value):
		return value
			
class HasExpressionTraits(HasTraits):
	

	
	def update(self):
		if not(hasattr(self,'_expressionDict')):
			self._expressionDict=dict()
		for k,v in self.traits().items():
			if isinstance(v.handler,TExpressionTrait) or isinstance(v.handler,TExpressionTraitDelegatesTo) :
				if not(self._expressionDict.has_key(k)):
					self._expressionDict[k]=dict()
				if isinstance(v.handler,TExpressionTrait) :
					self.updateExpressionTrait(k,v.handler)
				if isinstance(v.handler,TExpressionTraitDelegatesTo) :
					self.updateExpressionTraitDelegatesTo(k,v.handler)
		
	def updateExpressionTrait(self,name,handler):
		print "updating updateExpressionTrait"
		if not(self._expressionDict[name].has_key('hasmeta')):
			self._expressionDict[name]['hasmeta']=True
			try:
				self.add_class_trait('E_'+name,handler.mytrait)

				print "added meta trait"
			except Exception as e:
				print e
			try:
				#en=ExpressionTraitListener(self,name))
				#self.add_trait_listener(ExpressionTraitListener(self,name))
				#self._on_trait_change(lambda s,n: self.Echanged(name,n), 'E_'+name)
				#self._on_trait_change(lambda s,n: self.changed(name,n),name)
				print "Added listener for "+ 'E_'+name
			except Exception as e:
				print e
				
		self.updateExpressionTraitAll(name,handler)
		
	def updateExpressionTraitDelegatesTo(self,name,handler):
		print "updating updateExpressionTraitDelegatesTo"
		if not(self._expressionDict[name].has_key('hasmeta')):
			self._expressionDict[name]['hasmeta']=True
			try:
				myprefix=handler.prefix
				if myprefix is '':
					myprefix = name
				self.add_class_trait('E_'+name,DelegatesTo(handler.delegate,prefix=myprefix))
				print "added meta trait delegation " + handler.delegate + " -> " + myprefix
			except Exception as e:
				print e
			try:
				#self.add_trait_listener(ExpressionTraitListener(self,name))
				#self._on_trait_change(lambda s,n: self.Echanged(name,n), 'E_'+name)
				#self._on_trait_change(lambda s,n: self.changed(name,n),name)
				print "Added listener for " + 'E_'+name
			except Exception as e:
				print e
			
		self.updateExpressionTraitAll(name,handler)
	
	def updateExpressionTraitAll(self,name,handler):
		input = getattr(self,name)
		print name + " input: ", input
		#if self._expressionDict[name].has_key('Ecache'):
		#	if self._expressionDict[name]['Ecache']!=getattr(self,'E_'+name):
		#		self._expressionDict[name]['Ecache']=getattr(self,'E_'+name)
		#		del self._expressionDict[name]['expression']
		#		setattr(self,name,getattr(self,'E_'+name))
		#		return
		flag=False
		if isinstance(input,str) or isinstance(input,unicode):
			flag=True
			if self._expressionDict[name].has_key('expression'):
				#todo
				self._expressionDict[name]['expression'].set_expr(input)
			else :
				self._expressionDict[name]['expression']=self.variables.new_expression(input)
		elif isinstance(input,Expression):
			flag=True
			self._expressionDict[name]['expression']=input
		else :
			output=input
		
		if flag:
			output=self._expressionDict[name]['expression'].get_curr_value()
			if output is None:
				return
			
		#self._expressionDict[name]['cache']=output
		print "output: ", output
		if not(output==None):
			setattr(self,'E_'+name,output)
			
	def changed(self,name,value):
		print "We registered a change: " , name , " - " , value
		
	def Echanged(self,name,value):
		print "We registered a change: " , name , " - " , value
		setattr(self,name,value)
		
class Expression:
		def __init__(self,expr):
			self.set_expr(expr)
		def get_curr_value(self):
			return eval(self.expr)
		def set_expr(self,expr):
			self.expr=expr

class Variables:
		def new_expression(self,e):
			return Expression(e)

class Foo(HasExpressionTraits):
	variables=Variables()
	width=TExpression(Int)
	height=TExpression(Float)
	source=Instance(Bar)
	foo=TExpression(DelegatesTo('source'))
	mytrait = Either(List, Str)
	def __init__(self):
		self.source=Bar()
		

	

print Float
f=Foo()

f.width=5
f.update()
print  f.width
print  f.height
f.width=8

f.height=7
f.update()

print  f.width
f.update()

print f.foo

f.source.foo=17

print f.foo, f.E_foo	

f.E_foo= 15

print f.source.foo

print  f.width
f.update()

`

print "Other instance of Foo"
g=Foo()

g.update()

f.width=13
g.width=15
print f.width, g.width
	
