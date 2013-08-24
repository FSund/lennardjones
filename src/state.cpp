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
    cout << "State custom constructor" << endl;
    createBoxes(interactionLength);
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

void State::createBoxes(double interactionLength)
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

//    ////
//    for (Atom* atom : atoms)
//        for (uint i = 0; i < 3; i++)
//            if (atom->getPosition()(i) < 0 || atom->getPosition()(i) > size(i))
//                cout << "! OUT OF BOUNDS !" << endl;
//    ////

    // checking if any atoms have moved outside their boxes
    linkedList<Atom*> purgedAtoms;
    // "purgedAtoms" will have a random Atom* as "item", but since next == 0,
    // this won't cause any trouble
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

//    ////
//    cout << "After" << endl;
//    for (Atom* atom : atoms)
//    {
//        cout << atom->getPosition().t();
//        cout << atom->getVelocity().t();
//        cout << atom->getForce().t();
//    }
//    cout << endl;
//    for (Box* box : boxes)
//        cout << "Atoms in box #" << sub2ind3d(box->index, nBoxesVec) << " = " << box->nAtoms << endl;
//    ////

    // calculating the forces between the atoms
    for (Box* box : boxes)
    {
        box->resetForcesBool();
    }
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

//void State::testForces()
//{
//    const uint atomNr = 12;
//    const vec3 rvec = atoms[atomNr]->getPosition();

//    vec3 force = zeros<vec>(3);
//    vec3 drvec;
//    double dr2, dr6, scalarForce;

//    for (uint i = 0; i < nAtoms; i++)
//    {
//        if (i == atomNr) continue;

//        drvec = rvec - atoms[i]->getPosition();

//        // minimum image convention
//        for (uint j = 0; j < 3; j++)
//            if (drvec(j) > size(j)/2.0)
//                drvec(j) = drvec(j)/2.0;

//        force += norm(drvec, 2)*drvec;
//    }

//    cout << "total 'force'' on atom #" << atomNr << " = " << force.t() << endl;
//}
