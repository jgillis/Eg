%module helloworld
%include "std_string.i"
// Exceptions handling
%include "exception.i"
%exception {
  try {
  $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
    throw(e);
  } catch (const char* e) { // depreciated!!
    SWIG_exception(SWIG_RuntimeError, e);
    throw(e);
  }
}

%include "helloworld.hpp"
%{
#include "helloworld.hpp"
%}


