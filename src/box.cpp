#include <src/box.h>

Box::Box()
{
    cout << "! Box default constructor !" << endl;
}

Box::Box(vec3 boxPos, vec3 boxSize, uvec3 boxIndex, uvec3 nBoxesVec):
    pos(boxPos),
    size(boxSize),
    index(boxIndex),
    nBoxesVec(nBoxesVec)
{
//    cout << "Box custom constructor" << endl;
    m_forcesAreCalculated = false;
}

void Box::findNeighbours(const vec3 systemSize, const vector<Box*> boxes)
{
    /* This method finds all the neighbouring boxes of the box, and adds them to
     * a vector with pointers to the neighbours. We do NOT add ourselves to the list
     * of neighbours */

//    cout << "Finding neighbours of box #" << sub2ind3d(this->index, nBoxesVec) << ", " << this << endl;

    ivec3 tempIndex, dIndex; // boxIndex needs to be integer, since di,.. can be negative
    uvec3 boxIndex;
    vec3 displacement;
    Box* box;

    displacementVectors.clear();
    nNeighbours = 0;
    neighbours.reserve(26); // reserve space for max number of neighbouring boxes
    for (int di = -1; di <= 1; di++)
    {
        for (int dj = -1; dj <= 1; dj++)
        {
            for (int dk = -1; dk <= 1; dk++)
            {
                dIndex(0) = di;
                dIndex(1) = dj;
                dIndex(2) = dk;
                tempIndex = index + dIndex;

                // making sure the index is within the box-space we have,
                // i.e. applying periodic boundary conditions for the boxes
                for (uint i = 0; i < 3; i++)
                {
                    boxIndex(i) = uint( tempIndex(i) - int(floor(double(tempIndex(i))/double(nBoxesVec(i))))*nBoxesVec(i) );
                }
                box = boxes[sub2ind3d(boxIndex, nBoxesVec)];

                // don't want to add ourselves to the list of neighbours
                // testing (di == 0) etc. isn't good enough if we have a small system!!
                if (box == this) continue;

                // check to see if we already have this box/index in the list
                // (only possible for systems with less than 27 boxes, meaning
                // smaller than 6*6*6*(3 sigma))
                bool inList = false;
                for (uint i = 0; i < nNeighbours; i++)
                {
                    if (neighbours[i] == box)
                    {
                        inList = true;
                        break;
                    }
                }
                if (inList) continue;

                // adding the neighbour to the list
                neighbours.push_back(box);

                // the displacement vectors are used to make sure we follow the minimum image convention
                displacement(0) = (-(int(boxIndex(0)) < int(index(0))+di) + (int(boxIndex(0)) > int(index(0))+di))*systemSize(0);
                displacement(1) = (-(int(boxIndex(1)) < int(index(1))+dj) + (int(boxIndex(1)) > int(index(1))+dj))*systemSize(1);
                displacement(2) = (-(int(boxIndex(2)) < int(index(2))+dk) + (int(boxIndex(2)) > int(index(2))+dk))*systemSize(2);
                displacementVectors = join_rows(displacementVectors, displacement);

                nNeighbours++; // keeping track of the number of neighbours we have added
                // this is mainly important when we have small systems (with less
                // than 6x6x6 unit cells for the Argon system), when we have less
                // than 26 neighbours
            }
        }
    }

//    ////
//    cout << displacementVectors << endl;
//    cout << "We found " << nNeighbours << "/" << neighbours.size() << " neighbours." << endl;
//    cout << endl;
//    cout << "the neighbours of box #" << sub2ind3d(index, nBoxesVec) << " are " << endl;
//    for (Box* box : neighbours)
//        cout << sub2ind3d(box->index, nBoxesVec) << endl;
//    cout << endl;
//    ////
}

void Box::addAtom(Atom* atom)
{
    atomList.insertFirstItem(atom);
    nAtoms++;
}

void Box::purgeAtoms(linkedList<Atom*> &purgedAtoms)
{
    /* Runs through the linked list of atoms and checks which (if any) atoms
     * have moved outside the box. Those atoms are removed from the list, and
     * added to the list "purgedAtoms", which is returned */

    linkedList<Atom*>* runner = &atomList; // this will edit the atomlist directly
    vec3 atomPos;

    // running through the list
    while (runner->readNext() != 0)
    {
        atomPos = runner->readItem()->getPosition();
        // checking if the atom is outside the box
        if (
                atomPos(0) <  pos(0) ||
                atomPos(1) <  pos(1) ||
                atomPos(2) <  pos(2) ||
                atomPos(0) >= pos(0) + size(0) ||
                atomPos(1) >= pos(1) + size(1) ||
                atomPos(2) >= pos(2) + size(2)
            )
        {
            // dropping the atom from the list
            purgedAtoms.insertFirstItem(runner->item);
            runner->dropFirstItem();
            nAtoms--;

            // don't want to advance to the next item in the list, since we just
            // removed an item from the list
            continue;
        }
        runner = runner->next;
    }
}

void Box::flush()
{
    /* Removes all atoms from the box.
     * The box only has a linked list of pointers to atoms, so we only need to
     * reset this list, and reset the # of atoms to do this */

    atomList.next = 0;
    atomList.item = 0;
    nAtoms = 0;
}

void Box::calculateForces()
{
    vec3 forceOnAtom, atomPos, rvec, forceFromBox;
    Box* box;
    const linkedList<Atom*>* runner = &atomList; // pointer to constant linkedList<Atom*>, NOT constant pointer
    bool isMatrixAtom;

    while (runner->readNext() != 0)
    {
        atomPos = runner->readItem()->getPosition();
        isMatrixAtom = runner->readItem()->isMatrixAtom();
        forceOnAtom = zeros<vec>(3);

        // forces from all neighbouring boxes
        for (uint i = 0; i < nNeighbours; i++)
        {
            box = neighbours[i];

            // TODO: replace this if-test with sorted vector of neighbours, since
            // the same boxes have already been calculated each time
            if (box->forcesAreCalculated()) continue;

            // we use the displacement-vectors to ensure that we obey the minimum
            // image convention
            rvec(0) = atomPos(0) + displacementVectors(0,i);
            rvec(1) = atomPos(1) + displacementVectors(1,i);
            rvec(2) = atomPos(2) + displacementVectors(2,i);

            forceFromBox = box->calculateForceFromBox(rvec, isMatrixAtom);
            forceOnAtom(0) += forceFromBox(0);
            forceOnAtom(1) += forceFromBox(1);
            forceOnAtom(2) += forceFromBox(2);
        }
        // force from all atoms in this box
        forceOnAtom += calculateForceFromSelf(runner);

        // adding the force from the neighbouring boxes on this atom, to the atom
        runner->readItem()->addToForce(forceOnAtom);

        // advancing to the next atom
        runner = runner->readNext();
    }

    // since we use N2L when calculating the forces from the neighbouring boxes,
    // we don't need to calculate any forces from atoms in this box again --
    // so we mark this box as calculated, so that the loop above skips it
    m_forcesAreCalculated = true;
}

const vec3 Box::calculateForceFromBox(const vec3 &rvec, const bool &isMatrixAtom)
{
    vec3 drvec, forceFromAtom;
    vec3 forceFromBox = zeros<vec>(3);
    double dr2, dr6, scalarForce;

    const linkedList<Atom*>* runner = &atomList;

    while (runner->readNext() != 0)
    {
        if (isMatrixAtom && runner->readItem()->isMatrixAtom()) continue;

        drvec = rvec - runner->readItem()->getPosition();

//        ////
//        if (norm(drvec,2) < 0.01)
//            cout << "! drvec == 0" << endl;
//        ////

        dr2 = drvec(0)*drvec(0) + drvec(1)*drvec(1) + drvec(2)*drvec(2);
        dr6 = dr2*dr2*dr2;
        scalarForce = 24.0*(2.0 - dr6)/(dr6*dr6*dr2);

        forceFromAtom(0) = scalarForce*drvec(0);
        forceFromAtom(1) = scalarForce*drvec(1);
        forceFromAtom(2) = scalarForce*drvec(2);

        forceFromBox(0) += forceFromAtom(0);
        forceFromBox(1) += forceFromAtom(1);
        forceFromBox(2) += forceFromAtom(2);

        runner->readItem()->addToForce(-forceFromAtom);

        runner = runner->readNext();
    }

    return forceFromBox;
}

const vec3 Box::calculateForceFromSelf(const linkedList<Atom*>* runner)
{
    /* Calculates the force on at atom inside *this, from all other atoms inside
     * the box, except the ones that are before it in the list of atoms. We skip
     * the first part of the list with N2L.
     * "runner" will be a copy of the linked list we input to the function, so
     * we can operate freely on it (it's const anyway though...) */

    vec3 drvec, forceFromAtom;
    vec3 forceFromSelf = zeros<vec>(3);
    double scalarForce, dr2, dr6;

    const vec3 rvec = runner->readItem()->getPosition();
    runner = runner->readNext(); // don't want to calculate force between the same atom

    while (runner->readNext() != 0)
    {
        drvec = rvec - runner->readItem()->getPosition();

//        ////
//        if (norm(drvec,2) < 0.01)
//            cout << "! drvec == 0" << endl;
//        ////

        dr2 = drvec(0)*drvec(0) + drvec(1)*drvec(1) + drvec(2)*drvec(2);
        dr6 = dr2*dr2*dr2;
        scalarForce = 24.0*(2.0 - dr6)/(dr6*dr6*dr2);

        forceFromAtom(0) = scalarForce*drvec(0);
        forceFromAtom(1) = scalarForce*drvec(1);
        forceFromAtom(2) = scalarForce*drvec(2);

        forceFromSelf(0) += forceFromAtom(0);
        forceFromSelf(1) += forceFromAtom(1);
        forceFromSelf(2) += forceFromAtom(2);

        runner->readItem()->addToForce(-forceFromAtom);

        runner = runner->readNext();
    }

    return forceFromSelf;
}
