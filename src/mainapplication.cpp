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

    nUnitCells << 5 << endr
               << 5 << endr
               << 5;
    unitCellSize << 3.0 << endr
                 << 3.0 << endr
                 << 3.0;
    interactionLength = 3.0;
    State state = generator.createCrystal(structure, nUnitCells, unitCellSize, interactionLength);

    generator.saveState(&state, "state.0000.xyz");

    long idum = -1;
    generator.setTemperature(&state, 1, &idum);

    double dt = 0.005;
    Integrator integrator(&state);

    ostringstream filename;
    uint nCycles = 200;
    for (uint i = 1; i < 1+nCycles; i++)
    {
        cout << "Cycle " << i << " of " << nCycles << endl;
        integrator.stepForward(dt);
        filename.str(string());
        filename << "./state." << setfill('0') << setw(4) << i << ".xyz";
        generator.saveState(&state, filename.str());
    }
}
