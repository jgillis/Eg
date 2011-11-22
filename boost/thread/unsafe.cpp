#include <boost/thread.hpp>
#include <iostream>
#include <time.h>
#include <sys/time.h>

struct timeval time_temp;

/*
 Both threads will try to access time_temp concurrently
*/
double gettime() {
  gettimeofday(&time_temp, NULL);
  return time_temp.tv_sec + time_temp.tv_usec/1e6;
}

void task1() { 
    int max = 1e8;
    int dummy=0;
    for (int i=0;i<max;i++) {
      double time_here = gettime();
      while(gettime()-time_here < 5) {dummy ++;}
      std::cout << "foo " << i << ":" << (gettime()-time_here) << std::endl;
    };
}

void task2() { 
    int max = 1e8;
    int dummy=0;
    for (int i=0;i<max;i++) {
      double time_here = gettime();
      while(gettime()-time_here < 1) {dummy ++;}
      std::cout << " bar " << i << ":" << (gettime()-time_here) << std::endl;
    };
}

int main (int argc, char ** argv) {
    using namespace boost; 
    thread thread_1 = thread(task1);
    thread thread_2 = thread(task2);

    // do other stuff
    thread_2.join();
    thread_1.join();
    return 0;
}

