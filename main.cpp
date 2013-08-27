#include <src/mainapplication.h>

int main(int argc, char* argv[])
{
    MainApplication program;
    program.run(argc, argv);
//    program.test();

    return 0;
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
