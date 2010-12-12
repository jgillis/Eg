from enthought.traits.api import HasTraits, Code, Button
from enthought.traits.ui.api import View, Item, VGroup,\
HGroup, VSplit, HSplit, Tabbed

class Model(HasTraits):
	a = Code("print 'hello'")
	b = Button("click me")
 
	traits_view = View(HSplit(
	VGroup(
	Tabbed(
	Item('a'),
	Item('a'),
	Item('a')),
	Item('b')),
	VSplit(
	VGroup('b','b','b'),
	HGroup('a', show_border=True,
	label="traits is great")),
	dock="horizontal"
	),
	resizable=True,
	id="my.test.program.id")

m=Model()
m.configure_traits()
