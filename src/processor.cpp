#include <src/processor.h>

Processor::Processor(State *state, uvec3 nProcsVec, uint myRank, uint nProcs):
    state(state),
    nProcsVec(nProcsVec),
    myRank(myRank),
    nProcs(nProcs)
{
}

void Processor::setUpProcessors()
{

}

