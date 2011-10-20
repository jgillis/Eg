all: run run2
	ldd run
	LD_LIBRARY_PATH="foo" ldd run
	LD_LIBRARY_PATH="foo" ./run
	ldd run2
	./run2

foo/foo.o: foo/foo.cpp
	cd foo && gcc -c -fPIC foo.cpp -lstdc++

foo/libfoo.so: foo/foo.o
	cd foo && gcc -shared -Wl,-soname,libfoo.so -o libfoo.so  foo.o

run: foo/libfoo.so bar/bar.cpp
	cd bar && gcc bar.cpp -o bar -lstdc++ -I../foo -L../foo -lfoo
	cp bar/bar run 
	
run2: foo/libfoo.so bar/bar.cpp
	cd bar && gcc bar.cpp -o bar -lstdc++ -I../foo -L../foo -lfoo -Wl,-rpath,`pwd`/../foo
	cp bar/bar run2
	
clean:
	rm -rf foo/foo.o foo/libfoo.so run bar/bar
