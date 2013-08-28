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
    MainApplication();
    int run(int argc, char *argv[]);
    int benchmark(int argc, char *argv[]);
};

#endif // MAINAPPLICATION_H
