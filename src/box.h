#ifndef BOX_H
#define BOX_H

// forward declarations
class Generator;

#include <vector>

#include <armadillo>

#include <src/atom.h>
#include <src/linked_list.h>
#include <src/inlines.h>

using namespace std;
using namespace arma;

class Box
{
friend class Generator;

public:
    Box();
    Box(vec3 boxPos, vec3 boxSize, uvec3 boxIndex, uvec3 nBoxesVec);
    void findNeighbours(const vec3 systemSize, const vector<Box *> boxes);
    void addAtom(Atom *atom);
    void purgeAtoms(linkedList<Atom*> &purgedAtoms);
    void flush();

    void calculateForces();
    inline void resetForcesBool();

protected:
    const vec3 calculateForceFromBox(const vec3 &rvec, const bool &isMatrixAtom);
    const vec3 calculateForceFromSelf(const linkedList<Atom *> *runner);
    inline const bool &forcesAreCalculated() const;

    vec3 pos;
    vec3 size; // in state?
    uvec3 index;
    uvec3 nBoxesVec;
    mat displacementVectors; // the displacement vectors are used to make sure we follow the minimum image convention
//    bool empty; // not needed?
    linkedList<Atom*> atomList;
    uint nAtoms;
    vector<Box*> neighbours;
    uint nNeighbours;

    bool m_forcesAreCalculated;
};

inline const bool &Box::forcesAreCalculated() const
{
    return m_forcesAreCalculated;
}

inline void Box::resetForcesBool()
{
    m_forcesAreCalculated = false;
}

#endif // CBOX_H
