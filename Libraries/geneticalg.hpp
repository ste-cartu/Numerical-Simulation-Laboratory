#pragma once

#include <armadillo>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <iomanip>

#include "RandomGen/random.hpp"
#include "fmtlib.hpp"
#include "library.hpp"

using namespace std;
using namespace arma;


/****************************************************************************/
/*                        TRAVELING SALESMAN PROBLEM                        */
/****************************************************************************/



/*––––––––––––––––––––––––––––––– INDIVIDUAL –––––––––––––––––––––––––––––––*/


class Individual {

    public :
        Individual() {
            rnd_ = nullptr;
            len_ = 0, loss_ = 0;
        }
        Individual(Random* rnd) : rnd_(rnd) {len_ = 0, loss_ = 0;}
        Individual(Random* rnd, uvec index) : rnd_(rnd) {
            len_ = index.n_elem;
            index_ = index;
            loss_ = 0;
        }
        Individual(Random* rnd, int len) : rnd_(rnd) {
            len_ = len;
            index_ = linspace<uvec>(0, len_ - 1, len_);
            loss_ = 0;
        }
        ~Individual() {;}

        Individual& operator =(const Individual& ind);
        unsigned long long& operator[](size_t i) {return index_(i);}

        bool Check();                               // check individuals constraints
        void Permutation(int blk_length = 0);       // performs the permutation of two elements or two blocks within an individual
        void Shift();                               // shifts a subset of individual elements
        void Inversion();                           // inverts the order of an individual subset
        void SetLoss(double loss) {loss_ = loss;}
        void SetIndex(uvec index) {index_ = index;}

        Random* GetRnd() const {return rnd_;}
        int GetLen() const {return len_;}
        uvec GetIndex() const {return index_;}
        double GetLoss() const {return loss_;}

    private :
        Random* rnd_;           // random numbers generator
        int len_;               // size of the individual - number of cities to visit
        uvec index_;            // core of the individual: indices to its elements (cities): it is easier to carry and to do operations only with indices
        double loss_;           // evaluation of the loss function - length of the path through all the cities
}; 



/*––––––––––––––––––––––––––––––– POPULATION –––––––––––––––––––––––––––––––*/


class Population {

    public :
        Population() {
            rnd_ = nullptr;
            dim_ = 0, len_ = 0, norm_ = 0;        
        }
        Population(Random* rnd) : rnd_(rnd) {dim_ = 0, len_ = 0, norm_ = 0;}
        Population(Random* rnd, int dim, mat genes, int order = 2) : rnd_(rnd) {        // constructor that generates a random population
            dim_ = dim;
            len_ = genes.n_cols;
            norm_ = order;
            genes_ = genes;
            dist_.set_size(len_, len_);
            Individual ind(rnd_, len_);
            
            for(int i=0 ; i<dim_ ; i++) {
                dna_[i] = ind;
                for(int j=0 ; j<2*len_ ; j++) {
                    ind.Permutation(1);
                    ind.Permutation();
                    ind.Shift();
                    ind.Inversion();
                }
            }
        }
        ~Population() {;}

        Population& operator =(const Population& pop);
        Individual& operator[](size_t i) {return dna_.at(i);}

        void Distances();                       // stores all distances of the elements to speed up the calculation of the loss values
        void Losses();                          // evaluates the losses of all the individuals of the population
        void Crossover(int i1);                 // performs the crossover between two individuals of the population
        void Order(bool crescent = true);       // reorders the population individuals with respect to the loss values
        void SetDna(unordered_map<int, Individual> dna) {dna_ = dna;}

        Random* GetRnd() const {return rnd_;}
        int GetDim() const {return dim_;}
        int GetLen() const {return len_;}
        int GetNorm() const {return norm_;}
        mat GetGenes() const {return genes_;}
        mat GetDist() const {return dist_;}
        unordered_map<int, Individual> GetDna() const {return dna_;}

    private :
        Random* rnd_;                           // random numbers generator
        int dim_, len_, norm_;                  // size of the population, size of an individual, order of the norm in the loss function (distance)
        mat genes_, dist_;                      // matrices with all the cities and all the distances between them
        unordered_map<int, Individual> dna_;    // true population: a map of individuals with indices
};



/*––––––––––––––––––––––––––––––– EVOLUTION –––––––––––––––––––––––––––––––*/


class TSP {

    public:
        ~TSP() {
            rnd_ = nullptr;
            sel_ = 0, p_ = 0, n_gens_ = 0, n_mut_ = 0;
            dim_ = 0, len_ = 0, norm_ = 0;
            rank_ = 0, n_migr_ = 0;
            temp_i_ = 1., temp_f_ = temp_i_, temp_ = temp_i_;
        }
        TSP(Random* rnd, int rank = 0) : rnd_(rnd) {
            rank_ = rank;
            sel_ = 0, p_ = 0, n_gens_ = 0, n_mut_ = 0;
            dim_ = 0, len_ = 0, norm_ = 0;
            n_migr_ = 0;
            temp_i_ = 1., temp_f_ = temp_i_, temp_ = temp_i_;
        }

        Individual& operator[](size_t i) {return pop_[i];}

        void Init(const string input_file = "input.txt");
        void InitPopulation();                                          // initializes the starting population with respect to the problem type
        void SetPopulation(Population pop);                             // changes the actual population
        void Circle(double x_c = 0, double y_c = 0, double r = 1);      // generation of cities on a circle (given the coordinates of the centre and the radius)
        void Square(double x_v = 0, double y_v = 0, double l = 1);      // generation of cities on a circle (given the coordinates of the lower left vertex and of the side)
        void FromFile(const string input_file);                         // read the cities from an input file
        void Order() {pop_.Order();}                                    // calls the population reordering method: reorders the population individuals with respect to the loss values
        void Selection(int type = 0);                                   // generates a new population selecting with higher probabilities individuals with small losses
        void Mutations();                                               // performs crossover and other mutations on all the population, each mutation has its probability to happen
        // int Selectionn(int type = 0);
        // void Mutationss(int inde);

        string GetType() {return type_;}
        int GetLen() {return len_;}
        int GetDim() {return dim_;}
        int GetNorm() {return norm_;}
        int GetNGens() {return n_gens_;}
        int GetNMigr() {return n_migr_;}
        Population GetPop() {return pop_;}
        string GetCitiesFile() {return cities_file_;}
        double GetTempI() {return temp_i_;}
        double GetTempF() {return temp_f_;}
        double GetTemp() {return temp_;}
        void SetType(string type) {type_ = type; cities_file_ = "cities_" + type_ + ".tsv";}
        void SetTemp(double temp) {temp_ = temp;}


    private:
        Random* rnd_;                           // random numbers generator
        Population pop_;                        // current population
        vec probs_;                             // probabilities for crossover and mutations
        string type_, cities_file_;             // type of the problem: circle, square or italy and file with cities
        int sel_, p_, n_gens_, n_mut_;          // type of selection algorithm, exponent of the order-based selection (type 0), number of generations, number of possible mutations
        int dim_, len_, norm_;                  // size of the population, size of an individual (number of cities), order of the norm in the loss function (distance)
        int rank_, n_migr_;                     // rank of parallel process and migration number
        double temp_i_, temp_f_, temp_;         // initial, final and current temperatures for parallel tempering
};
