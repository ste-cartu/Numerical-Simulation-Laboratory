#include <string>

#include "../../../../Libraries/Simulator/SOURCE/system.hpp"
#include "../../../../Libraries/library.hpp"

int main() {

    const std::string path = "./Equilibration/";
    const std::string rnd_path = "../../../Libraries/RandomGen/";

    System SYS;
    SYS.initialize(path, rnd_path);
    SYS.initialize_properties(path);

    SYS.equilibration(path);
    SYS.finalize(path);

    return 0;
}