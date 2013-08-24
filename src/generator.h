#ifndef GENERATOR_H
#define GENERATOR_H

#include <armadillo>

#include <src/atom.h>
#include <src/state.h>
#include <src/lib.h>

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

    void setTemperature(State* state, const double &temperature, long *idum);

    State load(string filename);
    void saveState(State* state, const string &filename);
protected:

};

#endif // GENERATOR_H
