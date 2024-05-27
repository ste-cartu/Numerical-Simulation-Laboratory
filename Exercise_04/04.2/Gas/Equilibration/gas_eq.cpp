#include <string>

#include "../../../../Libraries/Simulator/SOURCE/system.hpp"
#include "../../../../Libraries/library.hpp"

int main() {

    const std::string path = "./Equilibration/";
    System SYS;
    SYS.initialize(path);
    SYS.initialize_properties(path);

    SYS.equilibration(path);
    SYS.finalize(path);

    return 0;
}