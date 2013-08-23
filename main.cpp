#include <src/mainapplication.h>

int main(int argc, char* argv[])
{
    MainApplication program;
    program.run(argc, argv);

//    string atomType;
//    vec pos;
//    long seed;
//    uint nAtoms;
//    linkedList<Atom*> atomList;

//    nAtoms = 1;
//    seed = -1;
//    atomType = "Ar";
//    Atom* atom;

//    // filling the atom list with atoms with random positions
//    for (uint i = 0; i < nAtoms; i++)
//    {
//        pos << ran0(&seed) << endr
//            << ran0(&seed) << endr
//            << ran0(&seed);
//        atom = new Atom(atomType, pos);
//        atomList.insertFirstItem(atom);

////        cout << "atom = " << atom << endl;
//        cout << "atom->getPosition().t() = " << atom->getPosition().t();
//    }
//    cout << endl;

//    linkedList<Atom*> newAtomList(atomList);
//    linkedList<Atom*>* runner = new linkedList<Atom*>;
//    runner = &newAtomList;

//    runner->dropFirstItem();

//    uint i = 0;
//    while (runner->readNext() != 0)
//    {
//        if (i == 4)
//        {
//            runner->dropFirstItem();
//            i++;
//            continue;
//        }
//        cout << "runner->readItem()->getPosition().t()" << runner->readItem()->getPosition().t();

//        runner = runner->next;
//        i++;
//    }

//    return 0;
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
