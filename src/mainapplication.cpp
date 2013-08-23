#include <src/mainapplication.h>

MainApplication::MainApplication()
{
}

void MainApplication::run(int argc, char *argv[])
{
    Generator generator;
    uint structure = 0;
    uvec3 nUnitCells;
    vec3 unitCellSize;
    double interactionLength;

    nUnitCells << 2 << endr
               << 2 << endr
               << 2;
    unitCellSize << 10.0 << endr
                 << 10.0 << endr
                 << 10.0;
    interactionLength = 3.0;
    State state = generator.createCrystal(structure, nUnitCells, unitCellSize, interactionLength);

    cout << state.getAtom(31)->getPosition().t() << endl;

//    double dt = 0.01;
//    Integrator integrator(&state);
//    integrator.stepForward(dt);
}
