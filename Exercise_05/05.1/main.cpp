#include <fstream>

#include "../../Libraries/metropolis.hpp"
#include "../../Libraries/blockaverage.hpp"


int main() {

    int nstep = static_cast<int>(1e6);
    unsigned int extr = static_cast<unsigned int>(1e7);
    unsigned int blocks = static_cast<unsigned int>(5e2);
    int unif = 0, gaus = 1;

    // equilibration
    double step_100 = 1., step_210 = 5.;
    double target = 0.5, prec = 0.001;

    hydrogen_100 prob_100;
    hydrogen_210 prob_210;
    arma::vec x_100 = {1.5,0,0};
    arma::vec x_210 = {5,0,0};
    Metropolis metro_100(step_100, unif, prob_100, x_100);
    Metropolis metro_210(step_210, unif, prob_210, x_210);
    //std::shared_ptr<Metropolis> metro_100 = std::make_shared<Metropolis>(step_100, unif, prob_100);
    //std::shared_ptr<Metropolis> metro_210 = std::make_shared<Metropolis>(step_210, unif, prob_210);


    /*––––––––––––––––––––––––––––––– GROUND STATE - UNIFORM SAMPLING –––––––––––––––––––––––––––––––*/

    fmt::print("\nGROUND STATE - UNIFORM SAMPLING\n");
    std::ofstream out("acceptance_100_unif.txt");
    out << "Target: " << target << "\tPrecision: " << prec << std::endl << std::endl;
    double step_100_unif = metro_100.SetAcceptance(target, prec, step_100, nstep, &out);
    out.close();

    // equilibration
    // out.open("eq_100_unif.csv");
    // out << "blocks,extractions,radius,error" << std::endl;
    // BA_Metro eq_100_unif(extr, blocks, metro_100, false);
    // eq_100_unif.Progressive(out);
    // out.close();

    // blocking average
    out.open("radius_100_unif.csv");
    out << "blocks,extractions,radius,error" << std::endl;
    std::ofstream sampl("sampling.csv");
    sampl << "x,y,z" << std::endl;
    BA_Metro radius_100_unif(extr, blocks, metro_100, true);
    // radius_100_unif.Equilib(extr);
    // fmt::print("\n");
    radius_100_unif.Progressive(out);
    out.close();

    // saving sampling
    rename("sampling.csv", "sampling_100_unif.csv");


    /*–––––––––––––––––––––––––––––– GROUND STATE - GAUSSIAN SAMPLING ––––––––––––––––––––––––––––––*/

    fmt::print("\n\nGROUND STATE - GAUSSIAN SAMPLING\n");
    out.open("acceptance_100_gaus.txt");
    out << "Target: " << target << "\tPrecision: " << prec << std::endl << std::endl;
    metro_100.Reset(step_100, gaus);
    metro_100.SetAcceptance(target, prec, step_100, nstep, &out);
    out.close();

    // equilibration
    // out.open("eq_100_gaus.csv");
    // out << "blocks,extractions,radius,error" << std::endl;
    // BA_Metro eq_100_gaus(extr, blocks, metro_100, false);
    // eq_100_gaus.Progressive(out);
    // out.close();

    // blocking average
    out.open("radius_100_gaus.csv");
    out << "blocks,extractions,radius,error" << std::endl;
    sampl.open("sampling.csv");
    sampl << "x,y,z" << std::endl;
    // radius_100_gaus.Equilib(extr);
    // fmt::print("\n");
    BA_Metro radius_100_gaus(extr, blocks, metro_100, true);
    radius_100_gaus.Progressive(out);
    out.close();

    // saving sampling
    rename("sampling.csv", "sampling_100_gaus.csv");
 

    /*––––––––––––––––––––––––––––––– EXCITED STATE - UNIFORM SAMPLING –––––––––––––––––––––––––––––––*/

    fmt::print("\n\nEXCITED STATE - UNIFORM SAMPLING\n");
    out.open("acceptance_210_unif.txt");
    out << "Target: " << target << "\tPrecision: " << prec << std::endl << std::endl;
    double step_210_unif = metro_210.SetAcceptance(target, prec, step_210, nstep, &out);
    out.close();

    // equilibration
    // out.open("eq_210_unif.csv");
    // out << "blocks,extractions,radius,error" << std::endl;
    // BA_Metro eq_210_unif(extr, blocks, metro_210, false);
    // eq_210_unif.Progressive(out);
    // out.close();

    // blocking average
    out.open("radius_210_unif.csv");
    out << "blocks,extractions,radius,error" << std::endl;
    sampl.open("sampling.csv");
    sampl << "x,y,z" << std::endl;
    // radius_210_unif.Equilib(extr);
    // fmt::print("\n");
    BA_Metro radius_210_unif(extr, blocks, metro_210, true);
    radius_210_unif.Progressive(out);
    out.close();

    // saving sampling
    rename("sampling.csv", "sampling_210_unif.csv");


    /*–––––––––––––––––––––––––––––– EXCITED STATE - GAUSSIAN SAMPLING ––––––––––––––––––––––––––––––*/

    fmt::print("\n\nEXCITED STATE - GAUSSIAN SAMPLING\n");
    out.open("acceptance_210_gaus.txt");
    out << "Target: " << target << "\tPrecision: " << prec << std::endl << std::endl;
    metro_210.Reset(step_210, gaus);
    metro_210.SetAcceptance(target, prec, step_210, nstep, &out);
    out.close();

    // equilibration
    // out.open("eq_210_gaus.csv");
    // out << "blocks,extractions,radius,error" << std::endl;
    // BA_Metro eq_210_gaus(extr, blocks, metro_210, false);
    // eq_210_gaus.Progressive(out);
    // out.close();

    // blocking average
    out.open("radius_210_gaus.csv");
    out << "blocks,extractions,radius,error" << std::endl;
    sampl.open("sampling.csv");
    sampl << "x,y,z" << std::endl;
    // radius_210_gaus.Equilib(extr);
    // fmt::print("\n");
    BA_Metro radius_210_gaus(extr, blocks, metro_210, true);
    radius_210_gaus.Progressive(out);
    out.close();

    // saving sampling
    rename("sampling.csv", "sampling_210_gaus.csv");


    /*–––––––––––––––––––––––––––––––––––– FAR FROM THE ORIGIN ––––––––––––––––––––––––––––––––––––*/

    fmt::print("\n\nFAR FROM THE ORIGIN\n");
    arma::vec x_far = {100., 100., 100.};
    Metropolis far_100(step_100_unif, unif, prob_100, x_far);
    Metropolis far_210(step_210_unif, unif, prob_210, x_far);

    out.open("far_100_unif.csv");
    out << "blocks,extractions,radius,error" << std::endl;
    sampl.open("sampling.csv");
    sampl << "x,y,z" << std::endl;
    BA_Metro far_100_unif(extr, blocks, far_100, true);
    far_100_unif.Progressive(out);
    out.close();
    rename("sampling.csv", "sampl_far_100_unif.csv");

    out.open("far_210_unif.csv");
    out << "blocks,extractions,radius,error" << std::endl;
    sampl.open("sampling.csv");
    sampl << "x,y,z" << std::endl;
    BA_Metro far_210_unif(extr, blocks, far_210, true);
    far_210_unif.Progressive(out);
    out.close();
    rename("sampling.csv", "sampl_far_210_unif.csv");


    fmt::print("\n\n");
 
    return 0;

}