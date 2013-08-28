#include <src/mainapplication.h>

MainApplication::MainApplication()
{
}

int MainApplication::run(int argc, char *argv[])
{
    Generator generator;
    uint structure = 0;
    uvec3 nUnitCells;
    vec3 unitCellSize;
    double interactionLength;

    nUnitCells << 12 << endr
               << 12 << endr
               << 12;
    double a = 5.72; // Angstrom
    unitCellSize << a << endr
                 << a << endr
                 << a;
    unitCellSize /= L0; // converting to MD units
    interactionLength = 3.0; // MD units
    State state = generator.createCrystal(structure, nUnitCells, unitCellSize, interactionLength);

    state.printInfo();
    generator.saveStateBox(&state, "state-0000.xyz");

    long idum = -1;
    double temperature = 5.0; // MD units
    generator.setTemperature(&state, temperature, &idum, 1);

    double dt = 0.005; // MD units
    Integrator integrator(&state);

    ostringstream filename;
    uint nCycles = 200;
    for (uint i = 1; i < 1+nCycles; i++)
    {
        cout << "Cycle " << i << " of " << nCycles << endl;
        integrator.stepForward(dt);
        filename.str(string());
        filename << "./state-" << setfill('0') << setw(4) << i << ".xyz";
        generator.saveStateBox(&state, filename.str());
    }

    return 0;
}

int MainApplication::benchmark(int argc, char *argv[])
{
    Generator generator;
    uint structure = 0;
    uvec3 nUnitCells;
    vec3 unitCellSize;
    double interactionLength;

    nUnitCells << 12 << endr
               << 12 << endr
               << 12;
    unitCellSize << 5.720 << endr
                 << 5.720 << endr
                 << 5.720; // Angstrom
    unitCellSize /= L0; // converting to MD units
    interactionLength = 3.0; // MD units
    State state = generator.createCrystal(structure, nUnitCells, unitCellSize, interactionLength);

//    state.printInfo();

//    generator.saveStateBox(&state, "state-0000.xyz");

    long idum = -1;
    double temperature = 1.0; // MD units
    generator.setTemperature(&state, temperature, &idum, 1);

    double dt = 0.005; // MD units
    Integrator integrator(&state);

//    ostringstream filename;
    uint nCycles = 70;
    for (uint i = 1; i < 1+nCycles; i++)
    {
//        cout << "Cycle " << i << " of " << nCycles << endl;
        integrator.stepForward(dt);
//        filename.str(string());
//        filename << "./state-" << setfill('0') << setw(4) << i << ".xyz";
//        generator.saveStateBox(&state, filename.str());
    }

    return 0;
}
