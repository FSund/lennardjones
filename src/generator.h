#ifndef GENERATOR_H
#define GENERATOR_H

#include <armadillo>

#include <src/state.h>
#include <src/atom.h>

using namespace std;
using namespace arma;

class Generator
{
public:
    Generator();
    State createCrystal(
            uint structure,
            uvec3 nUnitCells,
            vec3 unitCellSize,
            double interactionLength);
protected:

};

#endif // GENERATOR_H
