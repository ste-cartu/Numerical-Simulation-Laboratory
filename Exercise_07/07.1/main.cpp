#include <string>
#include <fstream>

#include "../../Libraries/Simulator/SOURCE/system.hpp"
#include "../../Libraries/library.hpp"
#include "../../Libraries/fmtlib.hpp"

int main() {

    const std::string path = "./";
    const std::string rnd_path = "../../Libraries/RandomGen/";
    
    System SYS;
    SYS.initialize(path, rnd_path);
    SYS.initialize_properties(path);

    SYS.equilibration(path);

    fmt::print("\nPotential tail correction: {:.5f}\n", SYS.get_vtail());
    fmt::print("Pressure tail correction: {:.5f}\n\n", SYS.get_ptail());
    ofstream out("OUTPUT/tail.txt");
    out << "Potential tail correction: " << SYS.get_vtail() << endl;
    out << "Pressure tail correction: " << SYS.get_ptail() << endl;
    out.close();

    fmt::print("SIMULATION\n");
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

    SYS.rename_files(path, path + "OUTPUT/SIMULATION/");
    fmt::print("\n\n");
    SYS.finalize(path);

    return 0;
}