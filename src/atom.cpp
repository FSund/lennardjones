#include "atom.h"

Atom::Atom()
{
    cout << "! Atom default constructor !" << endl;
}

Atom::Atom(
        const string &atomType_,
        const vec3 &position,
        const vec3 &velocity,
        const vec3 &force):
    position(position),
    velocity(velocity),
    force(force),
//    potEn(0.0),
//    pressure(0.0),
//    boundaryCrossings(zeros<ivec>(3,1)),
    matrixAtom(0),
    atomType(atomType_)
{
//    cout << "Atom custom constructor" << endl;
    setAtomType(atomType_);
}

Atom::~Atom()
{
    // default destructor
}

//double Atom::getPotEn() const
//{
//    return potEn;
//}

//double Atom::getPressure() const
//{
//    return pressure;
//}

//ivec3 Atom::getBoundaryCrossings() const
//{
//    return boundaryCrossings;
//}

void Atom::setPosition(const vec3 &position_)
{
    position = position_;
}

void Atom::setVelocity(const vec3 &velocity_)
{
    velocity = velocity_;
}

void Atom::setForce(const vec3 &force_)
{
    force = force_;
}

void Atom::setAtomType(const string &atomType_)
{
    if (atomType_ == "Ar")
    {
        atomType = atomType_;
        matrixAtom = false;
    }
    else if (atomType_ == "Ar_m")
    {
        atomType = atomType_;
        matrixAtom = true;
        velocity.zeros();
    }
    else // default
    {
        cout << endl << "! Unknown atom type, exiting !" << endl << endl;
        exit(1);
    }
}

void Atom::setMatrixAtom()
{
    atomType = "Ar_m";
    matrixAtom = true;
}

void Atom::setMovingAtom()
{
    atomType = "Ar";
    matrixAtom = false;
}

string Atom::getAtomType() const
{
    return atomType;
}
