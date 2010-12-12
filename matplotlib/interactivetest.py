plot(range(1000))

d=DragletInteractor()


def do(pos,p):
	if len(pos) < 2:
		return None
	b=zip(*pos)
	if p:
		p._x=b[0]
		p._y=b[1]
		draw()
	else:
		p=plot(b[0],b[1])[0]
	return p
	
	
d.setEvent('changed',do)




d.setEvent('changed',d.lineDrawAgent(lambda x: x*2))