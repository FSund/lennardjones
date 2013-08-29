#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <mpi.h>

#include <src/state.h>

class Processor
{
public:
    Processor(uint myRank, uint nProcs);
private:
    uint myRank;
    uint nProcs;
};


#endif // PROCESSOR_H
