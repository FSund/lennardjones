#include <src/state.h>

State::State()
{
    cout << "! State default constructor !" << endl;
}

State::State(
        vector<Atom*> atomVec,
        vec3 systemSize,
        double interactionLength):
    atoms(atomVec),
    nAtoms(atomVec.size()),
    size(systemSize)
{
//    cout << "State custom constructor" << endl;
    createBoxes(interactionLength);
    sortAtoms();
}

void State::printInfo()
{
    cout << "| — — — — — — — — — — — — — — — — — — — — — — — — — — — " << endl;
    cout << "| System info" << endl;
    cout << "| — — — — — — — — — — — — — — — — — — — — — — — — — — — " << endl;
    cout << "| nAtoms = " << nAtoms << endl;
    cout << "| Size of boxes  (MD): " << boxSize.t()
         << "|                (SI): " << boxSize.t()*L0
         << "| Size of system (MD): " << size.t()
         << "|                (SI): " << size.t()*L0
         << "| Number of boxes    : " << nBoxesVec.t();
    cout << "| — — — — — — — — — — — — — — — — — — — — — — — — — — — " << endl;
}

//State::~State()
//{
//    for (Box* box : boxes)
//    {
//        delete box;
//    }
//    for (Atom* atom : atoms)
//    {
//        delete atom;
//    }
//}

void State::createBoxes(double interactionLength)
{
    for (uint i = 0; i < 3; i++)
        nBoxesVec(i) = uint(floor(size(i)/interactionLength));

    boxSize = size/nBoxesVec;
    nBoxes = prod(nBoxesVec);

    boxes = vector<Box*>();
    boxes.reserve(nBoxes);
    boxes.resize(nBoxes); // using resize since we aren't using push_back to add items

    uvec3 boxIndex;
    vec3 boxPos;
    for (uint ix = 0; ix < nBoxesVec(0); ix++)
    {
        for (uint iy = 0; iy < nBoxesVec(1); iy++)
        {
            for (uint iz = 0; iz < nBoxesVec(2); iz++)
            {
                boxIndex << ix << endr
                         << iy << endr
                         << iz;
                boxPos = boxIndex%boxSize;
                Box* box = new Box(boxPos, boxSize, boxIndex, nBoxesVec);
                boxes[sub2ind3d(boxIndex, nBoxesVec)] = box;
            }
        }
    }

    // finding all the neighbouring boxes of each box
    for (uint i = 0; i < nBoxes; i++)
        boxes[i]->findNeighbours(size, boxes);
}

void State::sortAtoms()
{
    /* Sorts all atoms in "atoms" into the correct boxes */

    // flush all boxes before (probably not necessary)
    for (Box* box : boxes)
        box->flush();

    // sort atoms into boxes
    for (Atom* atom : atoms)
        putAtomInCorrectBox(atom);
}

void State::putAtomInCorrectBox(Atom *atom)
{
//    cout << "State::putAtomInCorrectBox" << endl;

    vec3 atomPos = atom->getPosition();
    uvec3 boxIndex;
    for (uint j = 0; j < 3; j++)
        boxIndex(j) = uint(floor(atomPos(j)/boxSize(j)));

    boxes[sub2ind3d(boxIndex, nBoxesVec)]->addAtom(atom);
}

void State::updateForces()
{
//    cout << "State::updateForces()" << endl;

    boundaryControl(); // PBC

    // checking if any atoms have moved outside their boxes
    vector<Atom*> purgedAtoms;
    // "purgedAtoms" will have a random Atom* as "item", but since next == 0,
    // this won't cause any trouble
    for (Box* box : boxes)
    {
        box->purgeAtoms(purgedAtoms);
    }

    // putting boxless atoms into their correct boxes
    for (auto it = purgedAtoms.begin(); it != purgedAtoms.end(); ++it)
    {
        putAtomInCorrectBox(*it);
    }

    // resetting stuff
    for (Box* box : boxes)
    {
        box->resetForcesBool();
    }
    for (Atom* atom : atoms)
    {
        atom->resetForce();
    }

    // calculating the forces between the atoms
    for (Box* box : boxes)
    {
        box->calculateForces();
    }

//    cout << "Exiting State::updateForces()" << endl;
}

void State::boundaryControl()
{
    vec3 pos;
    bool moved;
    for (Atom* atom : atoms)
    {
        moved = false;
        pos = atom->getPosition();
        for (uint i = 0; i < 3; i++)
        {
            if (pos(i) < 0.0 || pos(i) >= size(i))
            {
                pos(i) -= floor(pos(i)/size(i))*size(i);
                moved = true;
            }
        }
        if (moved)
        {
            atom->setPosition(pos);
        }
    }
}
