%module helloworld
%feature("autodoc","1");

%{
#include "fish.hpp"
#include "animal.hpp"
#include "vertebrate.hpp"
%}
%include "vertebrate.hpp"
%include "animal.hpp"
%template(species_fish) Animal<Fish>;
%include "fish.hpp"
