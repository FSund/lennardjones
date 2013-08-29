#ifndef STATE_H
#define STATE_H

// forward declarations
class Generator;
class Integrator;
class Processor;

#include <armadillo>

#include <src/atom.h>
#include <src/box.h>
#include <src/inlines.h>

using namespace std;
using namespace arma;

class State
{
friend class Generator;
friend class Integrator;
friend class Processor;

public:
    State();
    State(vector<Atom*> atomVec, vec3 systemSize, double interactionLength);
//    ~State();
    void printInfo();
    void createBoxes(double interactionLength);
    void sortAtoms();
    void putAtomInCorrectBox(Atom *atom);
    void addAtom(Atom* atom);

    void updateForces();
    void boundaryControl();
//    void testForces();

    inline const uint &getnAtoms() const;
    inline const Atom *readAtom(uint idx) const;

//    inline const vector<const Atom *> getAtoms() const;

protected:
    vector<Atom*> atoms;
    vector<Box*> boxes;
    uint nAtoms;
    uint nBoxes;
    vec3 size;
    vec3 boxSize;
    uvec3 nBoxesVec;
};

inline const uint &State::getnAtoms() const
{
    return nAtoms;
}

inline const Atom *State::readAtom(uint idx) const
{
    return atoms[idx];
}

//inline const vector<const Atom*> State::getAtoms() const
//{
//    vector<const Atom*> constAtoms;
//    constAtoms.reserve(nAtoms);
//    for (Atom* atom : atoms)
//    {
//        const Atom* constAtom = atom;
//        constAtoms.push_back(constAtom);
//    }
//    return constAtoms;
//}

#endif // STATE_H
