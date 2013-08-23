#ifndef STATE_H
#define STATE_H

// forward declarations
class Integrator;

#include <armadillo>

#include <src/atom.h>
#include <src/box.h>
#include <src/inlines.h>

using namespace std;
using namespace arma;

class State
{
friend class Integrator;

public:
    State();
    State(vector<Atom*> atomVec, vec3 systemSize, double interactionLength);
//    ~State();
    void createBoxes();
    void sortAtoms();
    void putAtomInCorrectBox(Atom *atom);
    void addAtom(Atom* atom);

    void updateForces();
    void boundaryControl();

    Atom *getAtom(uint idx);

protected:
    vector<Atom*> atoms;
    vector<Box*> boxes;
    uint nAtoms;
    uint nMovingAtoms;
    uint nBoxes;
    vec3 size;
    vec3 boxSize;
    uvec3 nBoxesVec;
    double interactionLength;
};

#endif // STATE_H
