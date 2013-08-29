#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <mpi.h>

#include <src/atom.h>
#include <src/box.h>
#include <src/state.h>

class Processor
{
public:
    Processor(State* state, uvec3 nProcsVec, uint myRank, uint nProcs);

private:
    void setUpProcessors();

    State* state;
    uvec3 nProcsVec;
    uint myRank;
    uint nProcs;

    vector<Box*> myBoxes; // boxes I find forces on
    vector<Box*> myNeighbours; // boxes I need positions from to find forces on myBoxes
};


#endif // PROCESSOR_H
