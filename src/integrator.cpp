#include <src/integrator.h>

Integrator::Integrator(State *state):
    state(state)
{
}

void Integrator::stepForward(double &dt)
{
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
    }
}
