all: client server

%: %.cpp
	g++ $< -o $@ -lboost_system -lpthread
	
clean:
	rm -f client server
