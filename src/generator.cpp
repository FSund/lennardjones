#include <src/generator.h>

Generator::Generator()
{
}

State Generator::createCrystal(
        uint structure,
        uvec3 nUnitCellsVec,
        vec3 unitCellSize,
        double interactionLength)
{
    // setting up the unit cell and lattice points
    mat latticePoints;
    switch (structure)
    {
    case 0:
        latticePoints << 0.0 << 0.5 << 0.0 << 0.5 << endr
                      << 0.0 << 0.5 << 0.5 << 0.0 << endr
                      << 0.0 << 0.0 << 0.5 << 0.5;

        latticePoints.row(0) *= unitCellSize(0);
        latticePoints.row(1) *= unitCellSize(1);
        latticePoints.row(2) *= unitCellSize(2);

        break;
    default:
        cout << "Unknown structure, exiting now!" << endl;
        exit(1);
    }

    // allocating and constructing atoms, and adding them to the state
    vec3 cellPos, atomPos;
    vector<Atom*> atoms;
    uint nAtoms = prod(nUnitCellsVec)*latticePoints.n_cols;
    atoms.reserve(nAtoms);
    for (uint ix = 0; ix < nUnitCellsVec(0); ix++)
    {
        for (uint iy = 0; iy < nUnitCellsVec(1); iy++)
        {
            for (uint iz = 0; iz < nUnitCellsVec(2); iz++)
            {
                for (uint i = 0; i < latticePoints.n_cols; i++) // number of atoms in each unit cell
                {
                    cellPos << unitCellSize(0)*ix << endr
                            << unitCellSize(1)*iy << endr
                            << unitCellSize(2)*iz;

                    atomPos = cellPos + latticePoints.col(i);
                    Atom* atom = new Atom("Ar", atomPos);
                    atoms.push_back(atom);

                    cout << "atom position = " << atom->getPosition().t() << endl;
                }
            }
        }
    }

    vec3 systemSize = nUnitCellsVec%unitCellSize;
    State state(atoms, systemSize, interactionLength);
    return state;
}
