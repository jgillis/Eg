r=300
a,b=mgrid[0:r,0:r]+0.0
p=sqrt(a**2+b**2)
matshow(p==floor(p))