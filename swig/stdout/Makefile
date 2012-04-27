all: greetcpp.py

greetcpp.o: greetcpp.cpp
	gcc -c -g -fPIC greetcpp.cpp -lstdc++

libgreetcpp.so: greetcpp.o
	gcc -shared -Wl,-soname,libgreetcpp.so -o libgreetcpp.so  greetcpp.o
	
greetcpp_wrap.cxx: greetcpp.i greetcpp.hpp
	swig2.0 -c++ -noexcept -python greetcpp.i 
	
greetcpp.py: setup.py greetcpp_wrap.cxx libgreetcpp.so
	python setup.py build_ext --inplace
	python test.py
	
clean:
	rm -f greetcpp.py greetcpp.pyc greetcpp_wrap.cxx greetcpp.o _greetcpp.so libgreetcpp.so
	rm -rf build
