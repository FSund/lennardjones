#include <src/box.h>

Box::Box(vec3 boxPos, vec3 boxSize, uvec3 boxIndex, uvec3 nBoxesVec):
    pos(boxPos),
    size(boxSize),
    index(boxIndex),
    nBoxesVec(nBoxesVec)
{
    neighbours.reserve(27); // reserve space for max number of neighbouring boxes
    m_forcesAreCalculated = false;
}

void Box::findNeighbours(const vec3 systemSize, const vector<Box*> boxes)
{
    nNeighbours = 0;
    ivec3 tempIndex, dIndex; // boxIndex needs to be integer, since di,.. can be negative
    uvec3 boxIndex;
    vec3 displacement;
    Box* box;
    displacementVectors.clear();

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
                    boxIndex(i) = tempIndex(i) - int(floor(double(tempIndex(i))/double(nBoxesVec(i))))*nBoxesVec(i);
                }
                box = boxes[sub2ind3d(boxIndex, nBoxesVec)];

                // check to see if we already have this box/index in the list
                // (only possible for systems with less than 27 boxes, meaning
                // smaller than 6*6*6*(3 sigma))
                bool inList = false;
                for (uint i = 0; i < nNeighbours; i++)
                {
                    if (neighbours[i] == box) // this will make sure "this" is included in the vector
                    {
                        inList = true;
                        break;
                    }
                }
                if (inList) continue;

                neighbours[nNeighbours] = box;

                // the displacement vectors are used to make sure we follow the minimum image convention
                displacement(0) = (-(boxIndex(0) < index(0)+di) + (boxIndex(0) > index(0)+di))*systemSize(0);
                displacement(1) = (-(boxIndex(1) < index(1)+dj) + (boxIndex(1) > index(1)+dj))*systemSize(1);
                displacement(2) = (-(boxIndex(2) < index(2)+dk) + (boxIndex(2) > index(2)+dk))*systemSize(2);
                displacementVectors = join_rows(displacementVectors, displacement);

                nNeighbours++; // keeping track of the number of neighbours we have added
                // this is most important when we have small systems (with less
                // than 6x6x6 unit cells for the Argon system), when we have less
                // than 26 neighbours
            }
        }
    }
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
    const linkedList<Atom*>* runner = &atomList;
    bool isMatrixAtom;

    // finding the forces from all neighbouring boxes (and itself)
    while (runner->readNext() != 0)
    {
        atomPos = runner->readItem()->getPosition();
        isMatrixAtom = runner->readItem()->isMatrixAtom();
        forceOnAtom.zeros();

        for (uint i = 0; i < nNeighbours; i++)
        {
            box = neighbours[i];

            // if we have already calculated the forces for this box (and its
            // neighbours) we just skip the force-calculation
            if (box->forcesAreCalculated()) continue;

            // else we calculate the force from the neighbouring box on the
            // current atom, and add the force from this atom to the forces on
            // the atoms in the neighbouring box

            rvec(0) = atomPos(0) + displacementVectors(0,i);
            rvec(1) = atomPos(1) + displacementVectors(1,i);
            rvec(2) = atomPos(2) + displacementVectors(2,i);

            forceFromBox = box->forceFromBox(rvec, isMatrixAtom);
            forceOnAtom(0) += forceFromBox(0);
            forceOnAtom(1) += forceFromBox(1);
            forceOnAtom(2) += forceFromBox(2);
        }
        // adding the force from the neighbouring boxes on this atom, to the atom
        runner->readItem()->addToForce(forceOnAtom);

        // advancing to the next atom
        runner = runner->readNext();
    }
}

const vec3 Box::forceFromBox(const vec3 &rvec, const bool &isMatrixAtom)
{
    vec3 forceFromBox, drvec, forceFromAtom;
    double dr2, dr6, scalarForce;
    const linkedList<Atom*>* runner = &atomList;
    while (runner->readNext() != 0)
    {
        if (runner->readItem()->isMatrixAtom() && isMatrixAtom) continue;

        drvec = rvec - runner->readItem()->getPosition();
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
    m_forcesAreCalculated = true;
    return forceFromBox;
}
