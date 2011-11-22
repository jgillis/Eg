#include <boost/thread.hpp>
#include <iostream>
#include <time.h>
#include <sys/time.h>

struct timeval time_temp1;
struct timeval time_temp2;

double gettime1() {
  gettimeofday(&time_temp1, NULL);
  return time_temp1.tv_sec + time_temp1.tv_usec/1e6;
}

double gettime2() {
  gettimeofday(&time_temp2, NULL);
  return time_temp2.tv_sec + time_temp2.tv_usec/1e6;
}

void task1() { 
    int max = 1e8;
    int dummy=0;
    for (int i=0;i<max;i++) {
      double time_here = gettime1();
      while(gettime1()-time_here < 5) {dummy ++;}
      std::cout << "foo " << i << ":" << (gettime1()-time_here) << std::endl;
    };
}

void task2() { 
    int max = 1e8;
    int dummy=0;
    for (int i=0;i<max;i++) {
      double time_here = gettime2();
      while(gettime2()-time_here < 1) {dummy ++;}
      std::cout << " bar " << i << ":" << (gettime2()-time_here) << std::endl;
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

