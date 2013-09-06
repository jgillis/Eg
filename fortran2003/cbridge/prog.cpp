#include <iostream>
#include <cstring>
#include <stdio.h>

#include "fancy.hpp"

// http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html

int main(int argc, char *argv[])
{
   int a = 3;
   int b;
   
   std::string ss = "foobar";
   char *s = new char [ss.length()+1];
   std::strcpy (s, ss.c_str());
   
   bool e = true;
   double d = 3.14;
   
   // Alternative: iso c binding
   fancy(&a,&b,s,&d,&e);
 
   std::cout << b << std::endl;
   return 0;
}
