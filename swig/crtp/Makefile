all: helloworld.py helloworld

fish.o: fish.cpp
	gcc -c -g -fPIC fish.cpp

animal.o: animal.cpp
	gcc -c -g -fPIC animal.cpp

vertebrate.o: vertebrate.cpp
	gcc -c -g -fPIC vertebrate.cpp

libhelloworld.so: fish.o animal.o vertebrate.o
	gcc -shared -Wl,-soname,libhelloworld.so -o libhelloworld.so  fish.o animal.o vertebrate.o
	
helloworld_wrap.cxx: libhelloworld.so helloworld.i animal.hpp fish.hpp
	swig -c++ -noexcept -python helloworld.i 
	
helloworld.py: setup.py helloworld_wrap.cxx libhelloworld.so
	python setup.py build_ext --inplace
	python testhelloworld.py
	
helloworld: helloworld.cpp
	gcc helloworld.cpp -o helloworld -L. -lhelloworld  -lstdc++

clean:
	rm -f helloworld.py helloworld.pyc helloworld_wrap.cxx helloworld *.o *.so
	rm -rf build
