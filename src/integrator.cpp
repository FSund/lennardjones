#include <src/integrator.h>

Integrator::Integrator()
{
    cout << "! Integrator default constructor !" << endl;
}

Integrator::Integrator(State *state):
    state(state)
{
    cout << "Integrator custom constructor" << endl;
}

void Integrator::stepForward(double &dt)
{
//    cout << "Integrator::stepForward()" << endl;

    vec3 position, velocity, force;
    for (Atom* atom : state->atoms)
    {
        if (atom->isMatrixAtom())
            continue;

        position = atom->getPosition();
        velocity = atom->getVelocity();
        force = atom->getForce();

        velocity += 0.5*force*dt;
        position += velocity*dt;
        atom->setPosition(position);
        atom->setVelocity(velocity);

        ////
//        cout << "position " << position.t();
//        cout << "force    " << force.t();
//        cout << "velocity " << velocity.t();
        ////
    }

    state->updateForces();

    for (Atom* atom : state->atoms)
    {
        if (atom->isMatrixAtom())
            continue;

        velocity = atom->getVelocity();
        force = atom->getForce();

        velocity += 0.5*force*dt;
        atom->setVelocity(velocity);

        ////
//        cout << "position " << atom->getPosition().t();
//        cout << "force    " << force.t();
//        cout << "velocity " << velocity.t();
        ////
    }

//    cout << "Exiting Integrator::stepForward()" << endl;
}
