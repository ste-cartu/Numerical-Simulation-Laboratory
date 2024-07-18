#include "geneticalg.hpp"


/****************************************************************************/
/*                        TRAVELING SALESMAN PROBLEM                        */
/****************************************************************************/




/*––––––––––––––––––––––––––––––– INDIVIDUAL –––––––––––––––––––––––––––––––*/


Individual& Individual :: operator =(const Individual& ind) {
    this->rnd_ = ind.GetRnd();
    this->len_ = ind.GetLen();
    this->index_ = ind.GetIndex();
    this->loss_ = ind.GetLoss();
    
    return *this;
}


bool Individual :: Check() {
    bool test = true;
    for(int i=0 ; i<len_-1 ; i++) {
        for(int j=i+1 ; j<len_ ; j++) {
            if(index_(i) == index_(j)) {
                test = false;
                fmt::print("\nERROR! Elements {} and {} are equal!\n", i, j);
            }
            if(!test) break;            // adding breaks to speed the cycles, I only need to know if there is an error or not
        }
        if(!test) break;
    }
    return test;
}


void Individual :: Permutation(int blk_length) {
    if(blk_length == 0) blk_length = static_cast<int>(rnd_->Rannyu(1, len_/2));      // static_cast<int> cuts off the decimals, so in this way I obtain an int value in the range [1,2,...,len_/2-1]
    int i1, i2, period_1 = 0, period_2 = 0;
    i1 = static_cast<int>(rnd_->Rannyu(1, len_));      
    do {
        i2 = static_cast<int>(rnd_->Rannyu(1, len_));
    } while(i2 == i1 or max(i1, i2) - min(i1, i2) < blk_length or (max(i1, i2) + blk_length >= len_ and max(i1, i2) + blk_length - len_ + 1 >= min(i1, i2)));
    
    uvec temp = zeros<uvec>(blk_length);
    for(int i=0 ; i<blk_length ; i++) {
        if(i1+i >= len_) period_1 = len_ - 1;
        if(i2+i >= len_) period_2 = len_ - 1;
        temp(i) = index_(i1 + i - period_1);
        index_(i1 + i - period_1) = index_(i2 + i - period_2);
        index_(i2 + i - period_2) = temp(i);
    }

    if(!Check()) {fmt::print("Error in Permutation!\n\n"); exit(1);}
}


void Individual :: Shift() {
    int blk_length, step, init;
    do {
        blk_length = static_cast<int>(rnd_->Rannyu(2, len_));        // int value in range [2,3,...,len_ - 1]
        step = static_cast<int>(rnd_->Rannyu(1, blk_length));        // int value in range [1,2,...,blk_length - 1]
        if(rnd_->Rannyu() < 0.5) step *= (-1);
        step %= blk_length;
    } while(step == 0);
    init = static_cast<int>(rnd_->Rannyu(1, len_));                  // int value in range [1,2,...,len_ - 1]

    uvec v = zeros<uvec>(blk_length);
    int period = 0;
    for(int i=0 ; i<blk_length ; i++) {
        if(init + i >= len_) period = len_ - 1;
        v(i) = index_(init + i - period);
    }
    v = shift(v, step);
    period = 0;
    for(int i=0 ; i<blk_length ; i++) {
        if(init + i >= len_) period = len_ - 1;
        index_(init + i - period) = v(i);
    }

    if(!Check()) {fmt::print("Error in Shift!\n\n"); exit(1);}
}


void Individual :: Inversion() {
    int init, blk_length;
    blk_length = static_cast<int>(rnd_->Rannyu(2, len_));    // int value in range [2,3,...,len_ - 1]
    init = static_cast<int>(rnd_->Rannyu(1, len_));          // int value in range [1,2,...,len_ - 1]
    
    uvec v = zeros<uvec>(blk_length);
    int period = 0;
    for(int i=0 ; i<blk_length ; i++) {
        if(init + i >= len_) period = len_ - 1;
        v(i) = index_(init + i - period);
    }
    period = 0;
    for(int i=0 ; i<blk_length ; i++) {
        if(init + i >= len_) period = len_ - 1;
        index_(init + i - period) = v(blk_length - 1 - i);
    }
    
    if(!Check()) {fmt::print("Error in Inversion!\n\n"); exit(1);}
}



/*––––––––––––––––––––––––––––––– POPULATION –––––––––––––––––––––––––––––––*/


Population& Population :: operator =(const Population& pop) {
    this->rnd_ = pop.GetRnd();
    this->dim_ = pop.GetDim();
    this->len_ = pop.GetLen();
    this->norm_ = pop.GetNorm();
    this->genes_ = pop.GetGenes();
    this->dist_ = pop.GetDist();
    this->dna_ = pop.GetDna();

    return *this;
}


void Population :: Distances() {
    for(int i=0 ; i<len_ ; i++)
        for(int j=0 ; j<len_ ; j++) {dist_(i,j) = norm(genes_.col(i) - genes_.col(j), norm_);}
}


void Population :: Losses() {
    double loss;
    for(int i=0 ; i<dim_ ; i++) {
        loss = dist_(dna_.at(i)[len_-1], dna_.at(i)[0]);         // distance between last and first elements, to close the circle
        for(int j=0 ; j<len_-1 ; j++) {
            loss += dist_(dna_.at(i)[j], dna_.at(i)[j+1]);
        }
        dna_.at(i).SetLoss(loss);
    }
}


void Population :: Crossover(int i1) {
    int i2;
    do {i2 = static_cast<int>(rnd_->Rannyu(0, dim_));} while(i2 == i1);      // int value in range [0,1,...,dim_ - 1]

    if(i1 >= dim_ or i1 < 0) {fmt::print("ERROR! Index {} out of population of size {}!\n\n", i1, dim_); return;}
    if(i2 >= dim_ or i2 < 0) {fmt::print("ERROR! Index {} out of population of size {}!\n\n", i2, dim_); return;}

    int cut = static_cast<int>(rnd_->Rannyu(2, len_));           // int value in range [2,3,...,len_ - 1]
    Individual copy1 = dna_.at(i1), copy2 = dna_.at(i2);
    int i_fill = cut, i_check = 1, i = 1;

    // filling the first individual
    while(i_fill < len_) {
        while(copy2[i_check] != copy1[i] and i < cut) {i++;}
        if(i == cut) {
            dna_.at(i1)[i_fill] = copy2[i_check]; 
            i_fill++;
        }
        i_check++;
        i = 1;
    }

    i_fill = cut, i_check = 1, i = 1;
    // filling the second individual
    while(i_fill < len_) {
        while(copy1[i_check] != copy2[i] and i < cut) {i++;}
        if(i == cut) {
            dna_.at(i2)[i_fill] = copy1[i_check]; 
            i_fill++;
        }
        i_check++;
        i = 1;
    }
    
    if(!dna_.at(i1).Check() or !dna_.at(i2).Check()) {fmt::print("Error in Crossover!\n\n"); exit(1);}
}


void Population :: Order(bool crescent) {
    Losses();
    vector<pair<int,double>> losses;
    for(int i=0 ; i<dim_ ; i++) {losses.push_back(make_pair(i, dna_.at(i).GetLoss()));}

    if(crescent) {sort(losses.begin(), losses.end(), [](const pair<int, double>& a, const pair<int, double>& b) {return a.second < b.second;});}
    else {sort(losses.begin(), losses.end(), [](const pair<int, double>& a, const pair<int, double>& b) {return a.second > b.second;});}

    unordered_map<int, Individual> sorted_dna;
    for(int i=0 ; i<dim_ ; i++) {sorted_dna[i] = dna_.at(losses[i].first);}

    dna_ = sorted_dna;
}




/*––––––––––––––––––––––––––––––– EVOLUTION –––––––––––––––––––––––––––––––*/


void TSP :: Init(const string input_file) {
    string param;
    int i_mut = 0;
    ifstream in(input_file);
    if(in.is_open()) {
        while(!in.eof()) {
            in >> param;
            if(param == "problem_type") {in >> type_; cities_file_ = "cities_" + type_ + ".tsv";}
            else if(param == "cities_number") {in >> len_;}
            else if(param == "distance_order") {in >> norm_;}
            else if(param == "selection_power") {in >> p_;}
            else if(param == "population_size") {in >> dim_;}
            else if(param == "generations_number") {in >> n_gens_;}
            else if(param == "migration_number") {in >> n_migr_;}
            else if(param == "mutations_number") {in >> n_mut_; probs_.resize(n_mut_);}
            else if(param == "crossover_prob") {in >> probs_(i_mut); i_mut++;}
            else if(param == "pair_permut_prob") {in >> probs_(i_mut); i_mut++;}
            else if(param == "block_permut_prob") {in >> probs_(i_mut); i_mut++;}
            else if(param == "shift_prob") {in >> probs_(i_mut); i_mut++;}
            else if(param == "inversion_prob") {in >> probs_(i_mut); i_mut++;}
            else if(param == "ENDINPUT") {break;}
            else {fmt::print("ERROR! Unknown input parameter: {}!\n\n", param); exit(1);}
        }
        in.close();
    }
    else {fmt::print("ERROR! Unable to open {}!\n\n", input_file); exit(1);}
    
    if(n_gens_ == 0) {fmt::print("ERROR! Program must initialize the number of generations!\n\n"); exit(1);}
    else if(n_mut_ == 0) {fmt::print("ERROR! Program must initialize the number of mutations!\n\n"); exit(1);}
    else if(dim_ == 0) {fmt::print("ERROR! Program must initialize the size of populations!\n\n"); exit(1);}
    else if(len_ == 0) {fmt::print("ERROR! Program must initialize the number of cities!\n\n"); exit(1);}
    else if(norm_ == 0) {fmt::print("ERROR! Program must initialize the order of the distances!\n\n"); exit(1);}
}


void TSP :: InitPopulation() {
    if(type_ == "circle") {this->Circle();}
    else if(type_ == "square") {this->Square();}
    else if(type_ == "italy") {this->FromFile(cities_file_);}
    else {fmt::print("ERROR! Unknown type of problem!\n\n"); exit(1);}
}


void TSP :: SetPopulation(Population pop) {
    dim_ = pop.GetDim();
    len_ = pop.GetLen();
    norm_ = pop.GetNorm();
    pop_ = pop;
}


void TSP :: Circle(double x_c, double y_c, double r) {
    ofstream out(cities_file_);
    if(!out.is_open()) {fmt::print("ERROR! Unable to open file {}!\n\n", cities_file_); exit(1);}

    out << fixed << setprecision(6) << "n\tx       \ty" << endl;
    mat cities(2, len_);
    double theta, x, y;
    for(int i=0 ; i<len_ ; i++) {
        theta = rnd_->Rannyu(0, 2*M_PI);
        x = r*cos(theta) + x_c;
        y = r*sin(theta) + y_c;
        cities(0,i) = x;
        cities(1,i) = y;
        out << i+1 << "\t" << x << "\t" << y << endl;
    }
    out.close();

    Population pop(rnd_, dim_, cities, norm_);
    pop.Distances();            // computing the distances between all cities only once
    pop.Losses();
    this->SetPopulation(pop);
}


void TSP :: Square(double x_v, double y_v, double l) {
ofstream out(cities_file_);
    if(!out.is_open()) {fmt::print("ERROR! Unable to open file {}!\n\n", cities_file_); exit(1);}

    out << fixed << setprecision(6) << "n\tx       \ty" << endl;
    mat cities(2, len_);
    double x, y;
    for(int i=0 ; i<len_ ; i++) {
        x = rnd_->Rannyu(x_v, l);
        y = rnd_->Rannyu(y_v, l);
        cities(0,i) = x;
        cities(1,i) = y;
        out << i+1 << "\t" << x << "\t" << y << endl;
    }
    out.close();

    Population pop(rnd_, dim_, cities, norm_);
    pop.Distances();            // computing the distances between all cities only once
    pop.Losses();
    this->SetPopulation(pop);
}


void TSP :: FromFile(const string input_file) {
    if(rnd_ == nullptr or dim_ == 0 or norm_ == 0) {
        fmt::print("ERROR! Before reading elements from a file, this program must initialize:\n");
        if(rnd_ == nullptr) {fmt::print("- The random numbers generator\n");}
        if(dim_ == 0) {fmt::print("- The population size\n");}
        if(norm_ == 0) {fmt::print("- The norm norm_er\n");}
        fmt::print("\n");
        exit(1);
    }

    ifstream in(input_file);
    if(!in.is_open()) {fmt::print("ERROR! Unable to open file {}!\n\n", input_file); exit(1);}
    int n;
    double x, y;
    mat cities(2, 0);
    if(in.is_open()) {
        vec col;
        string firstline;
        getline(in, firstline);
        while(!in.eof()) {
            in >> n >> x >> y;
            col = {x,y};
            cities.insert_cols(cities.n_cols, col);
        }
        in.close();
        cities.shed_col(cities.n_cols - 1);
    }
    else {fmt::print("ERROR! Unable to open {}!\n\n", input_file); exit(1);}
    if(len_ != n) {
        fmt::print("ERROR! Individual length mismatch between initialization file 'input.txt' and cities file '{}'!\n", input_file);
        fmt::print("In 'input.txt': {}   In '{}': {}\n\n", len_, input_file, n);
        exit(1);
    }
    Population pop(rnd_, dim_, cities, norm_);
    pop.Distances();            // computing the distances between all cities only once
    pop.Losses();
    this->SetPopulation(pop);    
}


void TSP :: Selection(int type) {
    int index;
    this->Order();
    unordered_map<int, Individual> dna;
    for(int i=0 ; i<dim_ ; i++) {       // norm_er-based selection
        if(type == 0) {index = static_cast<int>((double)dim_ * (pow(rnd_->Rannyu(), p_)));}
        else {fmt::print("ERROR! Unknown selection type!\n\n"); exit(1);}
        dna[i] = pop_[index];
    }
    pop_.SetDna(dna);
    this->Order();
}


int TSP :: Selectionn(int type) {
    int index;
    if(type == 0) {         // norm_er-based selection
        pop_.Order();   
        index = static_cast<int>((double)dim_ * (pow(rnd_->Rannyu(), p_)));
    }
    else {fmt::print("ERROR! Unknown selection type!\n\n"); exit(1);}
    return index;
}


void TSP :: Mutations() {
    for(int i=0 ; i<dim_ ; i++) {
        for(int j=0 ; j<n_mut_ ; j++) {
            if(rnd_->Rannyu() < probs_(j)) {
                if(j == 0) {pop_.Crossover(i);}
                else if (j == 1) {pop_[i].Permutation(1);}
                else if (j == 2) {pop_[i].Permutation();}
                else if (j == 3) {pop_[i].Shift();}
                else if (j == 4) {pop_[i].Inversion();}
                else {fmt::print("ERROR! Unknown mutation!\n\n"); exit(1);}
            }
        }
    }
}


void TSP :: Mutationss(int i) {
    for(int j=0 ; j<n_mut_ ; j++) {
        if(rnd_->Rannyu() < probs_(j)) {
            if(j == 0) {pop_.Crossover(i);}
            else if (j == 1) {pop_[i].Permutation(1);}
            else if (j == 2) {pop_[i].Permutation();}
            else if (j == 3) {pop_[i].Shift();}
            else if (j == 4) {pop_[i].Inversion();}
            else {fmt::print("ERROR! Unknown mutation!\n\n"); exit(1);}
        }
    }

}

