#include <string>

#include "../../Libraries/Simulator/SOURCE/system.hpp"
#include "../../Libraries/library.hpp"
#include "../../Libraries/fmtlib.hpp"

int main(int argc, char* argv[]) {

    if(argc != 2) {
        fmt::print("\nERROR! USAGE: {} Solid/Liquid/Gas\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    const std::string path = "./" + std::string(argv[1]) + "/Both/";
    
    System SYS;
    SYS.initialize(path);
    SYS.initialize_properties(path);

    SYS.equilibration(path);
    SYS.initialize_properties(path);

    fmt::print("\nSIMULATION\n");
    int blocks = SYS.get_nbl();
    int steps = SYS.get_nsteps();

    for(int i=0 ; i<blocks ; i++) {
        for(int j=0 ; j<steps ; j++){
            Progress_Bar(i*steps + j, blocks*steps - 1);
            SYS.step();
            SYS.measure();
        }
        SYS.averages(i+1, path);
        SYS.block_reset(i+1, path);
    }
    
    fmt::print("\nSIMULATION COMPLETED\n\n");
    SYS.finalize(path);

    return 0;
}