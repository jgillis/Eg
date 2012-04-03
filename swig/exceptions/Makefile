all: helloworld.py

helloworld.o: helloworld.cpp
	gcc -c -g -fPIC helloworld.cpp -lstdc++

libhelloworld.so: helloworld.o
	gcc -shared -Wl,-soname,libhelloworld.so -o libhelloworld.so  helloworld.o
	
helloworld_wrap.cxx: helloworld.i helloworld.hpp
	swig -c++ -noexcept -python helloworld.i 
	
helloworld.py: setup.py helloworld_wrap.cxx libhelloworld.so
	python setup.py build_ext --inplace
	python testhelloworld.py
	
clean:
	rm -f helloworld.py helloworld.pyc helloworld_wrap.cxx helloworld.o _helloworld.so libhelloworld.so
	rm -rf build
