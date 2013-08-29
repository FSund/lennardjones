#include <mpi.h>

#include <src/mainapplication.h>

int main(int argc, char* argv[])
{
    // MPI initialization
    int myRank, nProcs;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);

    MainApplication program(myRank, nProcs);
    program.run(argc, argv);
//    program.test();

    // MPI
    MPI_Finalize();

    return 0;
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
