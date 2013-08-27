#include <src/generator.h>

Generator::Generator()
{
    cout << "Generator default constructor." << endl;
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
                      << 0.0 << 0.0 << 0.5 << 0.5; // relative units

        // converting to MD-units
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
                }
            }
        }
    }

    vec3 systemSize = nUnitCellsVec%unitCellSize;
    State state(atoms, systemSize, interactionLength);
    return state;
}

void Generator::setTemperature(
        State* state,
        const double &temperature,
        long* idum,
        uint distribution)
{
    uint nAtoms = state->getnAtoms();
    mat velocities(3, nAtoms);

    // generating random velocities
    switch (distribution)
    {
    case 0:
        velocities = randu<mat>(3, nAtoms);
        velocities = (velocities - 0.5)*sqrt(temperature);
        break;
    case 1:
        velocities = randn<mat>(3, nAtoms);
        velocities *= sqrt(temperature);
        break;
    default:
        cout << "! Generator::setTemperature - Unknown random distribution !" << endl;
        exit(1);
    }

    // removing any linear momentum from the system
    vec3 momentum = sum(velocities, 1)/double(nAtoms);
    for(uint i = 0; i < nAtoms; i++)
    {
        velocities.col(i) -= momentum;
    }

    // sending the velocities to the atoms
    for (uint i = 0; i < nAtoms; i++)
    {
        state->atoms[i]->setVelocity(velocities.col(i));
    }
}

void Generator::saveState(State *state, const string &filename)
{
//    cout << "Generator::saveState" << endl;

    // // check and fix the extension of the filename
    // if (filename.substr(filename.find_last_of(".")) != ".xyz")
    // {
    //     filename.append(".xyz");
    // }

    bool indexing = true;
    bool saveForces = true;

    // open file stream
    ofstream ofile;
    ofile.open(filename.c_str());

    // check if we managed to open the file on disk
    if (!ofile)
    {
        cout << endl;
        cout << "! It seems like the '.xyz'-file couldn't be created. ";
        cout << "No .xyz-files will be saved. Please fix the folder or filename ";
        cout << "and rerun if you want to save the states." << endl;
        return;
    }

    uint nAtoms = state->getnAtoms();
    const Atom* atom;
    vec3 pos, vel, force;

    ofile << nAtoms << endl;
    ofile << "Comment" << endl;
    for (uint i = 0; i < nAtoms; i++)
    {
        atom = state->readAtom(i);
        ofile << atom->getAtomType();
        //ofile << scientific; // uncomment if you want "0.000000e+00" formatting
        ofile << setw(16) << setprecision(8);

        pos = atom->getPosition()*L0; // converting to SI units
        for (uint j = 0; j < 3; j++)
        {
            ofile << setw(18) << setprecision(8) << pos(j);
        }
        vel = atom->getVelocity(); // MD units (for now)
        for (uint j = 0; j < 3; j++)
        {
            ofile << setw(16) << setprecision(8) << vel(j);
        }
        if (saveForces)
        {
            force = atom->getForce();  // MD units (for now)
            for (uint j = 0; j < 3; j++)
            {
                ofile << setw(16) << setprecision(8) << force(j);
            }
            ofile << setw(16) << setprecision(8) << norm(force, 2);
        }
        if (indexing)
        {
            ofile << setw(16) << setprecision(8) << i+1;
        }
        ofile << endl;
    }
    ofile.close();

    //    cout << "Exiting Generator::saveState" << endl;
}

void Generator::saveStateBox(State *state, const string &filename)
{
    bool indexing = true;
    bool saveForces = true;

    // open file stream
    ofstream ofile;
    ofile.open(filename.c_str());

    // check if we managed to open the file on disk
    if (!ofile)
    {
        cout << endl;
        cout << "! It seems like the '.xyz'-file couldn't be created. ";
        cout << "No .xyz-files will be saved. Please fix the folder or filename ";
        cout << "and rerun if you want to save the states." << endl;
        return;
    }

    uint nAtoms = state->getnAtoms();
    const Atom* atom;
    vec3 pos, vel, force;

    ofile << nAtoms << endl;
    ofile << "Comment" << endl;

    const linkedList<Atom*>* runner;
    uint i = 0;
    for (Box* box : state->boxes)
    {
        runner = &box->atomList;
        while (runner->readNext() != 0)
        {
            atom = runner->readItem();
            ofile << atom->getAtomType();
            //ofile << scientific; // uncomment if you want "0.000000e+00" formatting
            ofile << setw(16) << setprecision(8);

            pos = atom->getPosition()*L0; // converting to SI units
            for (uint j = 0; j < 3; j++)
            {
                ofile << setw(18) << setprecision(8) << pos(j);
            }
            vel = atom->getVelocity(); // MD units (for now)
            for (uint j = 0; j < 3; j++)
            {
                ofile << setw(16) << setprecision(8) << vel(j);
            }
            if (saveForces)
            {
                force = atom->getForce();  // MD units (for now)
                for (uint j = 0; j < 3; j++)
                {
                    ofile << setw(16) << setprecision(8) << force(j);
                }
                ofile << setw(16) << setprecision(8) << norm(force, 2);
            }
            if (indexing)
            {
                ofile << setw(16) << setprecision(8) << i+1;
                ofile << setw(16) << setprecision(8) << sub2ind3d(box->index, state->nBoxesVec);
            }
            ofile << endl;
            i++;
            runner = runner->readNext();
        }
    }

    ofile.close();
}

