#include <fstream>
#include <vector>
#include <cmath>

#include "../../Libraries/fmtlib.hpp"
#include "../../Libraries/RandomGen/random.hpp"
#include "../../Libraries/randomwalk.hpp"

using namespace std;


int main(int argc, char** argv){

    fmt::print("\n");

    // initial check
    if (argc != 5) {
        fmt::print("ERROR! Program usage: {} <walk_step> <n_steps> <n_extractions> <n_blocks>\n\n", argv[0]);
        return -1;
    }

    double l = stod(argv[1]);                   // step length of the random walk
    unsigned int steps = stod(argv[2]);         // number of random walk steps 
    unsigned int extr = stod(argv[3]);          // number of extractions
    unsigned int blocks = stod(argv[4]);        // number of blocks

    if (extr%blocks != 0) {
        fmt::print("ERROR! <n_extractions> must be a mutiple of <n_blocks>\n\n");
        return -1;
    }

    unsigned int dim = extr/blocks;

    // output files with discrete random walks
    ofstream disc_out("discrete_distance.csv");
    disc_out << "step,distance,error" << endl;
    ofstream disc_ex("discrete_example.csv");
    disc_ex << "step,x,y,z" << endl;

    // output files with continue random walks
    ofstream cont_out("continue_distance.csv");
    cont_out << "step,distance,error" << endl;
    ofstream cont_ex("continue_example.csv");
    cont_ex << "step,x,y,z" << endl;

    // initialization of variables
    Walker w(l);
    vector <Walker> disc_walks(extr, w);
    vector <Walker> cont_walks(extr, w);
    Random rnd("../../Libraries/RandomGen/");

    double r, the, phi;
    double disc_block_dist = 0, disc_block_dist2 = 0;
    double disc_prog_dist = 0, disc_prog_dist2 = 0;
    double cont_block_dist = 0, cont_block_dist2 = 0;
    double cont_prog_dist = 0, cont_prog_dist2 = 0;
    double disc_err, cont_err;

    fmt::print("Random walk\n");
    for (int s=0 ; s<steps ; s++) {
        for (int i=0 ; i<blocks ; i++) {
            for (int j=0 ; j<dim ; j++) {
                r = rnd.Rannyu();
                the = rnd.Rannyu();
                phi = rnd.Rannyu();

                disc_walks[dim*i + j].Discrete_Step(r);
                disc_block_dist2 += disc_walks[dim*i + j].Distance2();
                cont_walks[dim*i + j].Continue_Step(the, phi);
                cont_block_dist2 += cont_walks[dim*i + j].Distance2();
            }

            disc_block_dist2 /= dim;
            disc_block_dist = sqrt(disc_block_dist2);
            cont_block_dist2 /= dim;
            cont_block_dist = sqrt(cont_block_dist2);

            disc_prog_dist += disc_block_dist;
            disc_prog_dist2 += disc_block_dist2;
            cont_prog_dist += cont_block_dist;
            cont_prog_dist2 += cont_block_dist2;
        }

        disc_prog_dist /= blocks;
        disc_prog_dist2 /= blocks;
        cont_prog_dist /= blocks;
        cont_prog_dist2 /= blocks;

        disc_err = Error(disc_prog_dist, disc_prog_dist2, blocks);
        cont_err = Error(cont_prog_dist, cont_prog_dist2, blocks);

        disc_out << s+1 << "," << disc_prog_dist << "," << disc_err << endl;
        disc_ex << s+1 << "," << disc_walks[0].X() << "," << disc_walks[0].Y() << "," << disc_walks[0].Z() << endl;
        cont_out << s+1 << "," << cont_prog_dist << "," << cont_err << endl;
        cont_ex << s+1 << "," << cont_walks[0].X() << "," << cont_walks[0].Y() << "," << cont_walks[0].Z() << endl;

        disc_block_dist = 0, disc_block_dist2 = 0;
        disc_prog_dist = 0, disc_prog_dist2 = 0;
        cont_block_dist = 0, cont_block_dist2 = 0;
        cont_prog_dist = 0, cont_prog_dist2 = 0;

        Progress_Bar(s, steps-1);

    }
    
    fmt::print("\n\n");
    
    return 0;
}