

from matplotlib.widgets import *

ax = subplot(111)

ax.plot(range[0,1])

def onselect(eclick, erelease):
	print 'startposition : (%f,%f)'%(eclick.xdata, eclick.ydata)
	print 'endposition : (%f,%f)'%(erelease.xdata, erelease.ydata)
	print 'used button : ', eclick.button

span = RectangleSelector(ax, onselect,drawtype='box')

show() 