#include <iostream>
#include <exception>
#include "helloworld.hpp"

std::string greeting() {
  return greetme();
}


std::string greetme() {
  return hi();
}

std::string hi() {
  throw std::exception();
  return "false";
}

