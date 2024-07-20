#include <string>
#include <fstream>

#include "../../Libraries/Simulator/SOURCE/system.hpp"
#include "../../Libraries/library.hpp"
#include "../../Libraries/fmtlib.hpp"


int main(int argc, char* argv[]) {

    if(argc != 2) {
        fmt::print("\nERROR! USAGE: {} Solid/Liquid/Gas\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    const std::string path = "./" + std::string(argv[1]) + "/";
    const std::string rnd_path = "../../Libraries/RandomGen/";
    fmt::print("\n");
    
    System SYS;
    SYS.initialize(path, rnd_path);
    ofstream out(path + "OUTPUT/set_acceptance.txt");
    if(strcmp(argv[1], "Gas") == 0) SYS.set_delta(2, 1e4, out);
    else SYS.set_acceptance(0.5, 0.01, 1e4, out, path);
    out.close();
    SYS.initialize_properties(path);

    fmt::print("\nSIMULATION\n");
    int steps = SYS.get_nsteps();
    int blocks = SYS.get_nbl();

    for(int i=0; i<blocks; i++) {
        for(int j=0; j<steps; j++) {
            Progress_Bar(i*steps + j, blocks*steps -1);
            SYS.step();
            SYS.measure();
        }
    SYS.averages(i+1, path);
    SYS.block_reset(i+1, path);
    }
    
    fmt::print("\n\n");
    SYS.finalize(path);

    return 0;
}