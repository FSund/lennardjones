#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H

#include <armadillo>

#include <src/state.h>
#include <src/generator.h>
#include <src/integrator.h>

using namespace std;
using namespace arma;

class MainApplication
{
public:
    MainApplication(uint myRank, uint nProcs);
    int run(int argc, char *argv[]);
    int benchmark(int argc, char *argv[]);

protected:
    uint myRank;
    uint nProcs;
};

#endif // MAINAPPLICATION_H
