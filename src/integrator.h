#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <armadillo>

#include <src/atom.h>
#include <src/state.h>

using namespace std;
using namespace arma;

class Integrator
{
public:
    Integrator(State* state);
    void stepForward(double &dt);

protected:
    State* state;
};

#endif // INTEGRATOR_H

