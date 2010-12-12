from pylab import *

# Er bestaat zoeits als de PolygonInteractor

# dragletInteractor
# definitions
#	* draglet: A clickable, dragabble, attribute-settable 'thing' that draws in an axes
# desired features:
#	* copy between different axes and figures, on a point-by-point basis
#	* bind to array.
# workflow:
# 	* i=dragletInteractor(fig=gcf()); put up the draglet in the toolbar; create scatter; 
#	* users clicks button, than clicks empty space. One draglet is created
#	* i.fixedpositions returns list of coordinates
# 


# demos: spline interp
# linear transformation point to point, same axes diferent axes, different figure
# cross




class DragletInteractor:
	def __init__(self,fig=None):
		print "Init"
		if fig == None:
			if gcf():
				fig=gcf()
			else:
				fig=figure()
		self.fig=fig
		self.tb=get_current_fig_manager().toolbar
		self.dragletCollection=[]
		self.button=self.tb._Button(text="draglet", file="/home/JG/eg/matplotlib/test.ppm",command=self.dragletButton)
		a=gca()
		self.axes=a
		self.dragging=None
		self.snapdist=35 # pixels
		self.tpChanged=None
		self.tpAdded=None
		self.tpRightclick=None
		self.tbMemory=None
		self.interpreter=lambda x: x
	def makeDraglet(self,pos,*args):
		figure(self.fig.number)
		draglet=Draglet(self,pos,*args)
		draglet.id=len(self.dragletCollection)
		self.dragletCollection.append(draglet)
		self.changed()
	def dragletButton(self,*args):
		mode='Draglet'
		if self.tb._active == mode:
		    self.tb._active = None
		else:
		    self.tb._active = mode
		self.tb.mode=mode
		if self.tb._idPress is not None:
		    self.tb._idPress=self.tb.canvas.mpl_disconnect(self.tb._idPress)
		    self.tb.mode = ''
		if self.tb._idRelease is not None:
		    self.tb._idRelease=self.tb.canvas.mpl_disconnect(self.tb._idRelease)
		    self.tb.mode = ''
		if  self.tb._active:
		    self.tb._idPress = self.tb.canvas.mpl_connect('button_press_event', self.press)
		    self.tb._idRelease = self.tb.canvas.mpl_connect('button_release_event', self.release)
		    self.tb._idMove = self.tb.canvas.mpl_connect('motion_notify_event', self.move)
		    self.tb._idResize = self.tb.canvas.mpl_connect('resize_event', self.resize)
		    self.tb.mode = mode
		    self.tb.canvas.widgetlock(self.tb)
		else:
		    self.tb.canvas.widgetlock.release(self.tb)
		self.tb.set_message(mode)
	def press(self,event):
		if event.inaxes:
			i,n=self.nearest((event.x,event.y))
			print "touched: #", i
			if event.button==3:
				if self.tpRightclick:
					changed=self.tpRightclick(n,event)
					if changed: self.changed()
			else:
				if n==None:
					self.makeDraglet((event.xdata,event.ydata))
				else:
					if self.dragging:
						self.dragging.dragEnd(event)
					self.dragging=n
					n.dragStart(event)
	def release(self,event):
		if self.dragging:
			self.dragging.dragEnd(event)
			self.dragging=None
	def move(self,event):
		if event.inaxes:
			if self.dragging:
				self.dragging.dragMove(event)
	def resize(self,event=None):
		#print "rescale"
		for d in self.dragletCollection:
			d.scale()
		draw()
	def nearest(self,screenPos):
		d=lambda a,b:  sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)
		dt=lambda a: d(a,screenPos)
		dist=map(dt,self.getScreenPositions())
		if len(dist)==0:
			return (None,None)
		m=min(dist)
		if m<=self.snapdist:
			i=dist.index(m)
			return (i,self.dragletCollection[i])
		else:
			return (None,None)
	def getPositions(self):
		return map(Draglet.getPosition,self.dragletCollection)
	def setPositions(self,pos):
		self.clear()
		for p in pos:
			self.makeDraglet(p)
	def clear(self):
		map(Draglet.destroy,self.dragletCollection)
		self.dragletCollection=[]
		self.tbMemory=None
		
	def getScreenPositions(self):
		return map(Draglet.getScreenPosition,self.dragletCollection)
	def scale(self):
		x=self.axes.get_xlim()
		y=self.axes.get_ylim()
		b=self.axes.get_window_extent()
		return (x[1]-x[0])/b.width(),(y[1]-y[0])/b.height()
	def changed(self):
		self.resize()
		if self.tpChanged:
			#try:
			self.tbMemory=self.tpChanged(self.getPositions(),self.tbMemory)
			#except:
			#	print "There was an error"
	def setEvent(self,type,method):
		knowntypes=('changed','added','rightclick')
		if type in knowntypes:
			if type=='changed':
				self.tpChanged=method
			if type=='added':
				self.tpChanged=method
			if type=='rightclick':
				self.tpRightclick=method
		else:
			raise "Unknow event type <" + str(type) +  "> \n known types: " + str(knowntypes)
	def setInterpreter(self,method):
		self.interpreter=method
	def lineDrawAgent(self,method,fig=None):
		if fig==None:
			fig=figure()
		def do(pos,p):
			if len(pos) < 2:
				return None
			figure(fig.number)
			b=method(array(pos).transpose())
			if p:
				p._x=b[0]
				p._y=b[1]
				draw()
			else:
				p=plot(b[0],b[1])[0]
			figure(self.fig.number)
			return p
		return do
		
	def wireframeDrawAgent(self,method,fig=None):
		import matplotlib.axes3d as p3
		if fig==None:
			fig=figure()
			ax = p3.Axes3D(fig)
		def do(pos,p):
			if len(pos) < 2:
				return None
			figure(fig.number)
			x,y,z=method(array(pos).transpose())
			if p:
				#p._x=x
				#p._y=y
				#p._y=z
				#draw()
				#cla()
				gca().collections.remove(p)
				p=ax.plot_wireframe(x,y,z)
				show()
			else:
				p=ax.plot_wireframe(x,y,z)
				show()
			figure(self.fig.number)
			return p
		return do
	def scatter3dDrawAgent(self,method,fig=None):
		import matplotlib.axes3d as p3
		if fig==None:
			fig=figure()
			ax = p3.Axes3D(fig)
		def do(pos,p):
			if len(pos) < 2:
				return None
			figure(fig.number)
			x,y,z=method(array(pos).transpose())
			if p:
				#p._x=ravel(x)
				#p._y=ravel(y)
				#p._y=ravel(z)
				#draw()
				#cla()
				gca().collections.remove(p._wrapped)
				del p
				p=ax.scatter3D(ravel(x),ravel(y),ravel(z))
				show()
			else:
				p=ax.scatter3D(ravel(x),ravel(y),ravel(z))
				show()
			figure(self.fig.number)
			return p
		return do

# todo: draw ellipse


class Draglet:
	def __init__(self,parent,pos,*args):
		print "Init Draglet"
		print "new", pos
		self.parent=parent
		self.axes=parent.axes
		self.pos=pos
		self.size=20 # pixels
		x=parent.scale()[0]*self.size
		y=parent.scale()[1]*self.size
		#self.patch=Circle(pos,self.size)
		self.patch=matplotlib.patches.Ellipse(pos,x,y)
		self.axes.add_patch(self.patch)
		self.dragMode=None
		self.draw()
		self.parent.changed()
	def destroy(self):
		if self in self.axes.patches:
			self.axes.patches.remove(self)
		del self
	def scale(self):
		self.patch.width=self.parent.scale()[0]*self.size
		self.patch.height=self.parent.scale()[1]*self.size
	def draw(self):
		draw()
	def setPos(self,pos):
		self.pos=pos
		self.patch.center=pos
		self.parent.changed()
	def getPosition(self):
		return self.pos
	def getScreenPosition(self):
		t=self.patch.get_transform()
		return t.xy_tup(self.getPosition())
	def dragStart(self,event):
		dragMode="simple"
		if event.key:
			dragMode=event.key
		self.dragMode=dragMode
		print "wh"
	def dragEnd(self,event):
		print "h!!"
		self.dragMode=None
	def dragMove(self,event):
		print "e"
		if self.dragMode:
			pos=(event.xdata,event.ydata)
			self.setPos(pos)
			self.draw()

# copy between different axes and figures, on a point-by-point basis
#~ clf()
#~ plot(range(10))

#~ d=DragletInteractor()


#~ def do(pos,p):
	#~ if len(pos) < 2:
		#~ return None
	#~ b=zip(*pos)
	#~ if p:
		#~ p._x=b[0]
		#~ p._y=b[1]
		#~ draw()
	#~ else:
		#~ p=plot(b[0],b[1])[0]
	#~ return p
	
	
#~ d.setEvent('changed',do)




#~ d.setEvent('changed',d.lineDrawAgent(lambda x: x*2))


# as a class (method definition), self.echo would have the same functionality as this lambda
# ppm is an image file format. use e.g. gimp to write it.
# the default image are located in: /usr/local/lib/python2.4/site-packages/matplotlib/mpl-data/images/


# c.draw(gcf().canvas.renderer)
