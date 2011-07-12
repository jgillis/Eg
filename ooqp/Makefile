SHELL = /bin/sh

OOQP=/usr/local
OOQPINCLUDEDIR=$(OOQP)/include
OOQPLIBDIR=$(OOQP)/lib


CXX      = c++

CXXFLAGS =-O 
CPPFLAGS =-I$(OOQPINCLUDEDIR)
LDFLAGS  =-L$(OOQPLIBDIR)
BLAS     = -lblas
MA27LIB  = -lMA27
MA57LIB  = -lMA57
FLIBS    = 
LIBS     = -looqpgensparse -looqpsparse  -looqpgondzio -looqpbase \
                $(BLAS) $(MA27LIB) $(MA57LIB)  $(FLIBS)

cxx_example.exe : example.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o

veryclean: clean
	rm -f *.exe

distclean: veryclean

.C.o:
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<
