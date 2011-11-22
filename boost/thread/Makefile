all: unsafe doublesafe mutexsafe

unsafe: unsafe.cpp
	g++ $< -o $@ -lboost_system -lboost_thread -lpthread  -lstdc++
	
doublesafe: doublesafe.cpp
	g++ $< -o $@ -lboost_system -lboost_thread -lpthread  -lstdc++
	
mutexsafe: mutexsafe.cpp
	g++ $< -o $@ -lboost_system -lboost_thread -lpthread  -lstdc++
