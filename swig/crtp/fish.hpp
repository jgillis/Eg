#ifndef FISH
#define FISH

#include <iostream>
#include "animal.hpp"
#include "vertebrate.hpp"

class Fish : public Animal<Fish>, public Vertebrate {
  public:
        Fish();
  	int swim();
  	int dnaBase(int base) { return 12; }
        using Vertebrate::numBones;
        int numBones(bool countAll) { return 12; }
};

#endif // FISH
