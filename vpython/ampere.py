from visual import *

print """
Electromagnetism: Faraday Law (v2.52) 2007-02-27
Rob Salgado (salgado@physics.syr.edu)

Electric Field vectors are blue. Magnetic Field vectors are red.

  The thick yellow vector representing
d|B|/dt ("time-rate-of-change-of-the-magnitude-of-the-magnetic-field")
is associated with the spatial arrangement of the electric field according to
the FARADAY Law (as evaluated on the yellow loop).
[The sense of circulation on the yellow loop (by the RightHandRule) determines
the direction of change of the magnetic field... OPPOSITE to your thumb.]

      CLICK the mouse to start and stop the animation
      TOGGLE: (f)araday
              (d)im-fields (n) color-scheme  (v)erbose"""


colorScheme=0          #key n (negative)
colorBackground=[color.black,color.white]
colorEdimmed=[(0.0,0,0.4),(0.5,0.5,1)]
scene.ambient=0.4
Ecolor=[color.blue,(0.5,0.5,1),color.yellow]
Ecolor[1]=colorEdimmed[colorScheme]

scene.background=colorBackground[colorScheme]
scene.background=color.black; Ecolor=[color.blue,(0,0,.4),color.yellow]
#scene.background=color.white; Ecolor=[color.blue,(0.9,0.9,0.9),color.yellow]
scene.title="FARADAY: Changing-Bs are associated with Curly-Es"
scene.range=(2.5,2.5,2.5)

showFaraday=0
dimFields=0



B=[]
B.append( arrow(pos=vector(0.25,0,0),axis=vector(0,0,1e-3), shaftwidth=0.04,fixedwidth=1, color=color.red)  )
B.append( arrow(pos=vector(-0.25,0,0),axis=vector(0,0,1e-3),shaftwidth=0.04,fixedwidth=1, color=color.red)  )  
B.append( arrow(pos=vector(0,0.25,0),axis=vector(0,0,1e-3), shaftwidth=0.04,fixedwidth=1, color=color.red)  )
B.append( arrow(pos=vector(0,-0.25,0),axis=vector(0,0,1e-3), shaftwidth=0.04,fixedwidth=1, color=color.red)  )  
B.append( arrow(pos=vector(0.25,0,-2),axis=vector(0,0,1e-3), shaftwidth=0.04,fixedwidth=1, color=color.red)  )
B.append( arrow(pos=vector(-0.25,0,-2),axis=vector(0,0,1e-3),shaftwidth=0.04,fixedwidth=1, color=color.red)  )  
B.append( arrow(pos=vector(0,0.25,-2),axis=vector(0,0,1e-3), shaftwidth=0.04,fixedwidth=1, color=color.red)  )
B.append( arrow(pos=vector(0,-0.25,-2),axis=vector(0,0,1e-3), shaftwidth=0.04,fixedwidth=1, color=color.red)  )  


N=8
dBdt=0.2
E=[]
Ebox=[]


for z in [0]:
    for r in [0.5]:
        for i in arange(0,N):
            theta=2.*pi*i/N
            theta_hat=vector(-sin(theta), cos(theta), 0) 
            Efield= -dBdt*theta_hat/r
            A=arrow(pos=vector(r*cos(theta),r*sin(theta),z) , axis=Efield,shaftwidth=0.04,fixedwidth=1,color=color.blue)
            E.append(A)
            Ebox.append( box(pos=A.pos+A.axis/4.,axis=A.axis,length=mag(A.axis)/2.,height=0.04,width=0.04,color=color.blue) )
    for r in [1,1.5]:
        for i in arange(0,N):
            theta=2.*pi*i/N
            theta_hat=vector(-sin(theta), cos(theta), 0) 
            Efield= -dBdt*theta_hat/r
            A=arrow(pos=vector(r*cos(theta),r*sin(theta),z) , axis=Efield,shaftwidth=0.04,fixedwidth=1,color=color.blue)
            E.append(A)
            Ebox.append( box(pos=A.pos+A.axis/4.,axis=A.axis,length=mag(A.axis)/2.,height=0.04,width=0.04,color=color.blue) )

for z in [-0.5,0.5,-1,1]:
    for r in arange (.5,1.5,.5):
        for i in arange(0,N):
            theta=2.*pi*i/N
            theta_hat=vector(-sin(theta), cos(theta), 0) 
            Efield= -dBdt*theta_hat/r
            A=arrow(pos=vector(r*cos(theta),r*sin(theta),z) , axis=Efield,shaftwidth=0.04,fixedwidth=1,color=color.blue)
            E.append(A)
            Ebox.append( box(pos=A.pos+A.axis/4.,axis=A.axis,length=mag(A.axis)/2.,height=0.04,width=0.04,color=color.blue) )



hcolor=Ecolor[2]

Bp=[]
for b in B:
    Bp.append( arrow(pos=b.pos+b.axis,axis=b.axis,length=dBdt,fixedwidth=1, color=hcolor, shaftwidth=0.07,headwidth=0.14, visible=showFaraday)  )  


Eloop_rad=mag(E[0].pos)
FaradayLoop=curve(color=hcolor, x=Eloop_rad*cos(2.*pi*arange(40)/40.), y=Eloop_rad*sin(2.*pi*arange(40)/40.), visible=showFaraday )

    



#I=cylinder(radius=0.04,pos=vector(0,0,-2),axis=vector(0,0,4), color=color.yellow)
#chgpos=[]
#chg=[]
#for i in arange(0,N):
#    chgpos.append(vector(I.pos+I.axis*i/N))
#    chg.append(sphere(pos=chgpos[-1],radius=0.05,color=I.color))


    
t=0

#Now... WHEN AN OBJECT IS PICKED,
#TRANSLATE THE scene.center TO THE OBJECT'S POSITION        
while 1:
    rate(10)
    t +=1
##    for i in arange(0,N):
##        chg[i].pos = chgpos[i]+(t%4)*vector(0,0,.125)
    bcount=0
    for b in B:
        b.length = (t%20)/10.+1e-3
        Bp[bcount].pos=b.pos+b.axis; bcount +=1
        
        
    if scene.mouse.clicked:
        scene.mouse.getclick()
        newPick=scene.mouse.pick
        if newPick !=None:
            ### ANIMATE TO SELECTED POSITION
            tempcolor=newPick.color
            newPick.color=color.yellow
            target=newPick.pos
            step=(target-scene.center)/20.
            for i in arange(1,20,1):
                rate(10)
                scene.center +=step
                scene.scale *= 1.037  #(1.037**19=1.99)
            newPick.color=tempcolor

    if scene.kb.keys: # is there an event waiting to be processed?
        s = scene.kb.getkey() # obtain keyboard information
        if s=='f':
            showFaraday +=1; showFaraday %=2; FaradayLoop.visible=showFaraday
            for i in Bp:
                i.visible=showFaraday

            if showFaraday==1:
                for i in arange(0,N):
                    E[i].color=hcolor
                    Ebox[i].color=hcolor
            else:
                for i in arange(0,N):
                    E[i].color=color.blue
                    Ebox[i].color=color.blue

        if s=='d':
            dimFields +=1; dimFields %=2; 

            for i in arange(N,len(E)):
                E[i].color=Ecolor[dimFields]
                Ebox[i].color=Ecolor[dimFields]
            for i in arange(1,4*N+1):
                E[-i].visible=(1-dimFields)
                Ebox[-i].visible=(1-dimFields)



        if s=='n':
            colorScheme = (colorScheme+1)%2 #TOGGLE colorScheme
            scene.background=colorBackground[colorScheme]
            Ecolor[1]=colorEdimmed[colorScheme]
            scene.background=colorBackground[colorScheme]

            for i in arange(N,len(E)):
                E[i].color=Ecolor[dimFields]
                Ebox[i].color=Ecolor[dimFields]
            for i in arange(1,4*N+1):
                E[-i].visible=(1-dimFields)
                Ebox[-i].visible=(1-dimFields)
