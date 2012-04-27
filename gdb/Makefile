all:  array_bounds

show.o: show.c
	g++ -g2 -c show.c
	
	
# -lm for linking to math library

array_bounds: array_bounds.c show.o
	g++ -g2 array_bounds.c -o array_bounds show.o -lm
	
