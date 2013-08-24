#ifndef ATOM_H
#define ATOM_H

#include <armadillo>

//const double SIGMA = 3.405; // Angstrom
//const double MASS = 39.948; // amu
//const double EPSILON = 0.010318; // eV

const double L0 = 3.405;    // Angstrom
const double t0 = 2156.9;   // fs
const double F0 = .30303;   // eV
const double E0 = 0.01038;  // eV
const double T0 = 119.74;   // K
const double pi = atan(1)*4;
const double forcemax = 1e4; // MD units

using namespace std;
using namespace arma;

class Atom
{
public:
    Atom();
    Atom(
            const string &atomType,
            const vec3 &position,
            const vec3 &velocity = zeros<vec>(3),
            const vec3 &force = zeros<vec>(3));
    ~Atom();

    inline const vec3 &getPosition() const;
    inline const vec3 &getVelocity() const;
    inline const vec3 &getForce() const;
//    inline void getData(vec3 &position_, vec3 &velocity_, vec3 &force_) const;

//    double getPotEn() const;
//    double getPressure() const;
//    ivec3 getBoundaryCrossings() const;

    void setPosition(const vec3 &position_);
    void setVelocity(const vec3 &velocity_);
    void setForce(const vec3 &newForce);

    inline void addToForce(const vec3 &addForce);
//    inline void addToBoundaryCrossings(ivec3 &addBoundaryCrossings);
//    inline void addToStatistics(double &addPot, double &addPressure);

    void setAtomType(const string &atomType_);
    void setMatrixAtom();
    void setMovingAtom();
    string getAtomType() const;
    inline const bool &isMatrixAtom() const;

protected:
    vec3 position;
    vec3 velocity;
    vec3 force;

//    double potEn;
//    double pressure;
//    ivec3 boundaryCrossings;

    bool matrixAtom;
    string atomType;
};

inline const vec3 &Atom::getPosition() const
{
    return position;
}

inline const vec3 &Atom::getVelocity() const
{
    return velocity;
}

inline const vec3 &Atom::getForce() const
{
    return force;
}

//inline void Atom::getData(vec3 &position_, vec3 &velocity_, vec3 &force_) const
//{
//    position_ = position;
//    velocity_ = velocity;
//    force_ = force;
//}

inline void Atom::addToForce(const vec3 &addForce)
{
    force(0) += addForce(0);
    force(1) += addForce(1);
    force(2) += addForce(2);
}

//inline void Atom::addToStatistics(double &addPot, double &addPressure)
//{
//    potEn    += addPot;
//    pressure += addPressure;
//}

//inline void Atom::addToBoundaryCrossings(ivec3 &addBoundaryCrossings)
//{
//    boundaryCrossings += addBoundaryCrossings;
//}

inline const bool &Atom::isMatrixAtom() const
{
    return matrixAtom;
}

#endif // ATOM_H
