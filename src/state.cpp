#include <src/state.h>

State::State()
{
    cout << "State::State() -- Using default constructor!" << endl;
}

State::State(
        vector<Atom*> atomVec,
        vec3 systemSize,
        double interactionLength):
    atoms(atomVec),
    nAtoms(atomVec.size()),
    size(systemSize),
    interactionLength(interactionLength)
{
    createBoxes();
    sortAtoms();
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

void State::createBoxes()
{
    for (uint i = 0; i < 3; i++)
        nBoxesVec(i) = floor(size(i)/interactionLength);

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

//    cout << "box dimensions (MD): " << boxDimensions.t()
//         << "               (SI): " << boxDimensions.t()*L0
//         << "system size    (MD): " << systemSize.t()
//         << "               (SI): " << systemSize.t()*L0
//         << "number of boxes    :    " << nBoxes << endl << endl;
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
    vec3 atomPos = atom->getPosition();
    uvec3 boxIndex;
    for (uint j = 0; j < 3; j++)
        boxIndex(j) = uint(floor(atomPos(j)/boxSize(j)));

    boxes[sub2ind3d(boxIndex, nBoxesVec)]->addAtom(atom);
}

void State::updateForces()
{
    boundaryControl();

    // checking if any atoms have moved outside their boxes
    linkedList<Atom*> purgedAtoms;
    for (Box* box : boxes)
    {
        box->purgeAtoms(purgedAtoms);
    }

    // putting boxless atoms into their correct boxes
    Atom* atom;
    while (purgedAtoms.readNext() != 0)
    {
        atom = purgedAtoms.readItem();
        putAtomInCorrectBox(atom);

        purgedAtoms = *purgedAtoms.next;
    }

    // calculating the forces between the atoms
    for (Box* box : boxes)
    {
        box->resetForcesBool();
    }
    for (Box* box : boxes)
    {
        box->calculateForces();
    }
}

void State::boundaryControl()
{
    vec3 pos;
    for (Atom* atom : atoms)
    {
        pos = atom->getPosition();
        for (uint i = 0; i < 3; i++)
        {
            pos(i) -= floor(pos(i)/size(i))*size(i);
        }
    }
}

Atom *State::getAtom(uint idx)
{
    cout << "nAtoms = " << nAtoms << endl;
    cout << "nAtoms2 = " << atoms.size() << endl;
    return atoms[idx];
}
