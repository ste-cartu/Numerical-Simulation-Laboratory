#include <string>

#include "../../Libraries/Simulator/SOURCE/system.hpp"
#include "../../Libraries/library.hpp"
#include "../../Libraries/fmtlib.hpp"

int main() {

    fmt::print("\nSIMULATION\n");

    const std::string path = "./";
    System SYS;
    SYS.initialize(path);
    SYS.initialize_properties(path);
    SYS.block_reset(0);

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
    
    SYS.finalize(path);
    fmt::print("\nSIMULATION COMPLETED\n\n");

    return 0;
}