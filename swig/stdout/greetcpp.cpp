#include <iostream>
#include <fstream>
#include <exception>
#include "greetcpp.hpp"
#include <stdio.h>

//FILE * stdout_;

void greeting() {
  std::cout << "greeting on stdout c++" << std::endl;
  std::cerr << "greeting on stderr c++" << std::endl;
  
  fprintf(stdout, "greeting on stdout c\n");
  fprintf(stderr, "greeting on stderr c\n");
}


std::ofstream filestr;
std::streambuf *backup;
void captureCout(const std::string & name) {
  std::streambuf *psbuf;
  filestr.open (name.c_str());
  backup = std::cout.rdbuf();     // back up cout's streambuf

  psbuf = filestr.rdbuf();   // get file's streambuf
  std::cout.rdbuf(psbuf);         // assign streambuf to cout
}

void releaseCout() {
  std::cout.rdbuf(backup);
  filestr.close();
}

fpos_t pos;
int fd;

void captureStdout(const std::string & name) {
  // Linux only
  fgetpos(stdout, &pos);
  fd = dup(fileno(stdout));
  freopen(name.c_str(), "w", stdout);
}

void releaseStdout() {
  fflush(stdout);
  dup2(fd, fileno(stdout));
  close(fd);
  clearerr(stdout);
  fsetpos(stdout, &pos);
}

