all: prog

prog: prog.cpp fancy.hpp libfancy.so
	gcc prog.cpp -fPIC -o prog -L. -lstdc++ -lfancy
	
libfancy.so: fancy.f90
	gfortran -shared -fPIC -o libfancy.so fancy.f90
	
libfancy.a: fancy.f90
	gfortran -static -o libfancy.so fancy.f90
	
clean:
	rm -rf *.so *~ *.a
