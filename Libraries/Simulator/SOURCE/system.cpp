/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/



#include <cmath>
#include <cstdlib>
#include <string>
#include <format>
#include "system.hpp"

using namespace std;
using namespace arma;

void System :: step(){ // Perform a simulation step
  if(_sim_type == 0) this->Verlet();  // Perform a MD step
  else for(int i=0; i<_npart; i++) {
    double r = _rnd.Rannyu();
    double index = r*_npart;
    int ind = static_cast<int>(index); 
    this->move(ind);
    } // Perform a MC step on a randomly choosen particle
  _nattempts += _npart; //update number of attempts performed on the system
  return;
}

void System :: Verlet(){
  double xnew, ynew, znew;
  for(int i=0; i<_npart; i++){ //Force acting on particle i
    _fx(i) = this->Force(i,0);
    _fy(i) = this->Force(i,1);
    _fz(i) = this->Force(i,2);
  }
  for(int i=0; i<_npart; i++){ //Verlet integration scheme
    xnew = this->pbc( 2.0 * _particle(i).getposition(0,true) - _particle(i).getposition(0,false) + _fx(i) * pow(_delta,2), 0);
    ynew = this->pbc( 2.0 * _particle(i).getposition(1,true) - _particle(i).getposition(1,false) + _fy(i) * pow(_delta,2), 1);
    znew = this->pbc( 2.0 * _particle(i).getposition(2,true) - _particle(i).getposition(2,false) + _fz(i) * pow(_delta,2), 2);
    _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0 * _delta));
    _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0 * _delta));
    _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0 * _delta));
    _particle(i).acceptmove(); // xold = xnew
    _particle(i).setposition(0, xnew);
    _particle(i).setposition(1, ynew);
    _particle(i).setposition(2, znew);
  }
  _naccepted += _npart;
  return;
}

double System :: Force(int i, int dim){
  double f=0.0, dr;
  vec distance;
  distance.resize(_ndim);
  for (int j=0; j<_npart; j++){
    if(i != j){
      distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
      distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
      dr = sqrt( dot(distance,distance) );
      if(dr < _r_cut){
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }
  return f;
}

void System :: move(int i){ // Propose a MC move for particle i
  if(_sim_type == 3){ //Gibbs sampler for Ising
    double delta_E = 2.0 * ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin()) + _H);
    double acceptance = 1.0 / (1.0 + exp(-_beta*delta_E));
    if(_rnd.Rannyu() < acceptance) _particle(i).setspin(1);
    else _particle(i).setspin(-1);
    _naccepted++;
  } else {           // M(RT)^2
    if(_sim_type == 1){       // LJ system
      vec shift(_ndim);       // Will store the proposed translation
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
      }
      _particle(i).translate(shift, _side);  //Call the function Particle::translate
      if(this->metro(i)){ //Metropolis acceptance evaluation
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback(); //If translation is rejected, restore the old configuration
    } else {                  // Ising 1D
      if(this->metro(i)){     //Metropolis acceptance evaluation for a spin flip involving spin i
        _particle(i).flip();  //If accepted, the spin i is flipped
        _naccepted++;
      }
    }
  }
  return;
}

bool System :: metro(int i){ // Metropolis algorithm
  bool decision = false;
  double delta_E, acceptance;
  if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  else delta_E = 2.0 * _particle(i).getspin() * ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  acceptance = exp(-_beta*delta_E);
  if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
  return decision;
}

void System :: set_delta(double delta, int nstep, ofstream& out){

  _naccepted = 0, _nattempts = 0;
  double fraction;
  _delta = delta;
  for(int j=0 ; j<nstep ; j++) this->step();
  fraction = (double)_naccepted/(double)_nattempts;
  fmt::print("Delta: {:.1f}\tAcceptance: {:.4f}\n", _delta, fraction);
  out << "Delta: " << setprecision(1) << _delta << "\tAcceptance: " << setprecision(4) << fraction << endl;

}

double System :: set_acceptance(double target, double prec, int nstep, ofstream& out, const string path){

  fmt::print("SETTING ACCEPTANCE:\n");
  int counter = 1;
  double fraction, step = _delta, corr = _delta;

  do {
    this->initialize(path);
    _naccepted = 0, _nattempts = 0;
    _delta = step;

    for(int j=0 ; j<nstep ; j++) this->step();
    fraction = (double)_naccepted/(double)_nattempts;
    fmt::print("{}) Delta: {:.5f}\tAcceptance: {:.4f}\n", counter, _delta, fraction);
    out << counter << ") Delta: " << setprecision(5) << _delta << "\tAcceptance: " << setprecision(4) << fraction << endl;
    counter++;

    if(corr > prec*step) corr /= 2.0;
    else corr = prec*step;

    if(fraction < target - prec) step -= corr;
    else if(fraction > target + prec) step += corr;

} while((fraction < target - prec) or (fraction > target + prec));

return _delta;

}

double System :: Boltzmann(int i, bool xnew){
  double energy_i=0.0;
  double dx, dy, dz, dr;
  for (int j=0; j<_npart; j++){
    if(j != i){
      dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      if(dr < _r_cut){
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }
  return 4.0 * energy_i;
}

double System :: pbc(double position, int i){ // Enforce periodic boundary conditions
  return position - _side(i) * rint(position / _side(i));
}

int System :: pbc(int i){ // Enforce periodic boundary conditions for spins
  if(i >= _npart) i = i - _npart;
  else if(i < 0)  i = i + _npart;
  return i;
} 

void System :: initialize(const string path, const string rnd_path){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from rnd_path/primes64001.in a pair of numbers to be used to initialize the RNG
  ifstream Primes(rnd_path + "primes64001.in");
  if(!Primes.is_open()) {fmt::print("PROBLEM: unable to open primes64001.in!\n\n");}
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed(rnd_path + "seed.in");
  if(!Seed.is_open()) {fmt::print("PROBLEM: unable to open seed.in!\n");}
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);

  ofstream couta(path + "OUTPUT/acceptance.csv"); // Set the heading line in file ../OUTPUT/acceptance.csv
  couta << "N_BLOCK,ACCEPTANCE" << endl;
  couta.close();

  ifstream input(path + "INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open(path + "OUTPUT/output.csv");
  string property;
  double delta;
  while ( !input.eof() ){
    input >> property;
    if( property == "SIMULATION_TYPE" ){
      input >> _sim_type;
      if(_sim_type > 1){
        input >> _J;
        input >> _H;
      }
      if(_sim_type > 3 or _sim_type < 0){
        cerr << "PROBLEM: unknown simulation type" << endl;
        exit(EXIT_FAILURE);
      }
      if(_sim_type == 0)      coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"   << endl;
      else if(_sim_type == 1) coutf << "LJ MONTE CARLO (NVT) SIMULATION"          << endl;
      else if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION"  << endl;
      else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION"  << endl;
    } else if( property == "RESTART" ){
      input >> _restart;
    } else if( property == "TEMP" ){
      input >> _temp;
      _beta = 1.0/_temp;
      coutf << "TEMPERATURE = " << _temp << endl;
    } else if( property == "NPART" ){
      input >> _npart;
      _fx.resize(_npart);
      _fy.resize(_npart);
      _fz.resize(_npart);
      _particle.set_size(_npart);
      for(int i=0; i<_npart; i++){
        _particle(i).initialize();
        if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // to randomize the spin configuration
      }
      coutf << "NPART = " << _npart << endl;
    } else if( property == "RHO" ){
      input >> _rho;
      _volume = _npart/_rho;
      _side.resize(_ndim);
      _halfside.resize(_ndim);
      double side = pow(_volume, 1.0/3.0);
      for(int i=0; i<_ndim; i++) _side(i) = side;
      _halfside=0.5*_side;
      coutf << "SIDE = ";
      for(int i=0; i<_ndim; i++){
        coutf << setw(12) << _side[i];
      }
      coutf << endl;
    } else if( property == "R_CUT" ){
      input >> _r_cut;
      coutf << "R_CUT = " << _r_cut << endl;
    } else if( property == "DELTA" ){
      input >> delta;
      coutf << "DELTA = " << delta << endl;
      _delta = delta;
    } else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS = " << _nblocks << endl;
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS = " << _nsteps << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  this->read_configuration(path);
  this->initialize_velocities(path);
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}

void System :: initialize_velocities(const string path){
  if(_restart and _sim_type==0){
    ifstream cinf;
    cinf.open(path + "INPUT/CONFIG/velocities.dat");
    if(cinf.is_open()){
      double vx, vy, vz;
      for(int i=0; i<_npart; i++){
        cinf >> vx >> vy >> vz;
        _particle(i).setvelocity(0,vx);
        _particle(i).setvelocity(1,vy);
        _particle(i).setvelocity(2,vz);
      }
    } else cerr << "PROBLEM: Unable to open INPUT file velocities.dat"<< endl;
    cinf.close();
  } else {
    vec vx(_npart), vy(_npart), vz(_npart);
    vec sumv(_ndim);
    sumv.zeros();
    for (int i=0; i<_npart; i++){
      vx(i) = _rnd.Gauss(0.,sqrt(_temp));
      vy(i) = _rnd.Gauss(0.,sqrt(_temp));
      vz(i) = _rnd.Gauss(0.,sqrt(_temp));
      sumv(0) += vx(i);
      sumv(1) += vy(i);
      sumv(2) += vz(i);
    }
    for (int idim=0; idim<_ndim; idim++) sumv(idim) = sumv(idim)/double(_npart);
    double sumv2 = 0.0, scalef;
    for (int i=0; i<_npart; i++){
      vx(i) = vx(i) - sumv(0);
      vy(i) = vy(i) - sumv(1);
      vz(i) = vz(i) - sumv(2);
      sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
    }
    sumv2 /= double(_npart);
    scalef = sqrt(3.0 * _temp / sumv2);   // velocity scale factor 
    for (int i=0; i<_npart; i++){
      _particle(i).setvelocity(0, vx(i)*scalef);
      _particle(i).setvelocity(1, vy(i)*scalef);
      _particle(i).setvelocity(2, vz(i)*scalef);
    }
  }
  if(_sim_type == 0){
  double xold, yold, zold;
    for (int i=0; i<_npart; i++){
      xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0);
      yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
      zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);
      _particle(i).setpositold(0, xold);
      _particle(i).setpositold(1, yold);
      _particle(i).setpositold(2, zold);
    }
  }
  return;
}

void System :: initialize_properties(const string path){ // Initialize data members used for measurement of properties

  string property;
  int index_property = 0;
  _nprop = 0;

  _measure_penergy  = false; //Defining which properties will be measured
  _measure_kenergy  = false;
  _measure_tenergy  = false;
  _measure_pressure = false;
  _measure_gofr     = false;
  _measure_magnet   = false;
  _measure_cv       = false;
  _measure_chi      = false;

  ifstream input(path + "INPUT/properties.dat");
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "POTENTIAL_ENERGY" ){
        ofstream coutp(path + "OUTPUT/potential_energy.csv");
        coutp << "BLOCK,ACTUAL_PE,PE_AVE,ERROR" << endl;
        coutp.close();
        _nprop++;
        _measure_penergy = true;
        _index_penergy = index_property;
        index_property++;
        _vtail = ((8.0 * M_PI * _rho ) / 3.0) * ((1.0/(3.0 * pow(_r_cut, 9))) - (1.0/pow(_r_cut, 3)));
      } else if( property == "KINETIC_ENERGY" ){
        ofstream coutk(path + "OUTPUT/kinetic_energy.csv");
        coutk << "BLOCK,ACTUAL_KE,KE_AVE,ERROR" << endl;
        coutk.close();
        _nprop++;
        _measure_kenergy = true;
        _index_kenergy = index_property;
        index_property++;
      } else if( property == "TOTAL_ENERGY" ){
        ofstream coutt(path + "OUTPUT/total_energy.csv");
        coutt << "BLOCK,ACTUAL_TE,TE_AVE,ERROR" << endl;
        coutt.close();
        _nprop++;
        _measure_tenergy = true;
        _index_tenergy = index_property;
        index_property++;
      } else if( property == "TEMPERATURE" ){
        ofstream coutte(path + "OUTPUT/temperature.csv");
        coutte << "BLOCK,ACTUAL_T,T_AVE,ERROR" << endl;
        coutte.close();
        _nprop++;
        _measure_temp = true;
        _index_temp = index_property;
        index_property++;
      } else if( property == "PRESSURE" ){
        ofstream coutpr(path + "OUTPUT/pressure.csv");
        coutpr << "BLOCK,ACTUAL_P,P_AVE,ERROR" << endl;
        coutpr.close();
        _nprop++;
        _measure_pressure = true;
        _index_pressure = index_property;
        index_property++;
        _ptail = ((32.0 * M_PI * _rho ) / 3.0) * ((1.0/(3.0 * pow(_r_cut, 9))) - (1.0/(2.0 * pow(_r_cut, 3))));
      } else if( property == "GOFR" ){
        ofstream coutgr(path + "OUTPUT/gofr_blocks.csv");
        coutgr << "BLOCK,GOFR_AVE,ERROR" << endl;
        coutgr.close();
        coutgr.open(path + "OUTPUT/gofr_final.csv");
        coutgr << "STEP,DISTANCE,AVE_GOFR,ERROR" << endl;
        coutgr.close();
        input>>_n_bins;
        _nprop+=_n_bins;
        _bin_size = (_halfside.min() )/(double)_n_bins;
        _measure_gofr = true;
        _index_gofr = index_property;
        index_property+= _n_bins;
      } else if( property == "MAGNETIZATION" ){
        ofstream coutpr(path + "OUTPUT/magnetization.csv");
        coutpr << "BLOCK,ACTUAL_M,M_AVE,ERROR" << endl;
        coutpr.close();
        _nprop++;
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property == "SPECIFIC_HEAT" ){
        ofstream coutpr(path + "OUTPUT/specific_heat.csv");
        coutpr << "BLOCK,ACTUAL_CV,CV_AVE,ERROR" << endl;
        coutpr.close();
        _nprop++;
        _measure_cv = true;
        _index_cv = index_property;
        index_property++;
      } else if( property == "SUSCEPTIBILITY" ){
        ofstream coutpr(path + "OUTPUT/susceptibility.csv");
        coutpr << "BLOCK,ACTUAL_X,X_AVE,ERROR" << endl;
        coutpr.close();
        _nprop++;
        _measure_chi = true;
        _index_chi = index_property;
        index_property++;
      } else if( property == "ENDPROPERTIES" ){
        ofstream coutf;
        coutf.open(path + "OUTPUT/output.csv",ios::app);
        coutf << "Reading properties completed!" << endl;
        coutf.close();
        break;
      } else cerr << "PROBLEM: unknown property" << endl;
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl;

  // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2,_global_av3
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}

void System :: finalize(const string path){
  this->write_configuration(path);
  this->write_velocities(path);
  const string outpath = path + "OUTPUT/";
  _rnd.SaveSeed(outpath);
  ofstream coutf;
  coutf.open(path + "OUTPUT/output.csv",ios::app);
  coutf << "Simulation completed!" << endl;
  coutf.close();
  return;
}

void System :: write_configuration(const string path){ // Write current configuration as a .xyz file in directory ../OUTPUT/CONFIG/
  ofstream coutf;
  if(_sim_type < 2){
    coutf.open(path + "OUTPUT/CONFIG/config.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  " 
              << setw(16) << _particle(i).getposition(0,true)/_side(0)          // x
              << setw(16) << _particle(i).getposition(1,true)/_side(1)          // y
              << setw(16) << _particle(i).getposition(2,true)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    this->write_velocities(path);
  } else {
    coutf.open(path + "OUTPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++) coutf << _particle(i).getspin() << " ";
    coutf.close();
  }
  return;
}

void System :: write_XYZ(int nconf, const string path){ // Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
  ofstream coutf;
  coutf.open(path + "OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");
  if(coutf.is_open()){
    coutf << _npart << endl;
    coutf << "#Comment!" << endl;
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  " 
            << setw(16) << _particle(i).getposition(0,true)          // x
            << setw(16) << _particle(i).getposition(1,true)          // y
            << setw(16) << _particle(i).getposition(2,true) << endl; // z
    }
  } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
  coutf.close();
  return;
}

void System :: write_velocities(const string path){
  ofstream coutf;
  coutf.open(path + "OUTPUT/CONFIG/velocities.dat");
  if(coutf.is_open()){
    for(int i=0; i<_npart; i++){
      coutf << setw(16) << _particle(i).getvelocity(0)          // vx
            << setw(16) << _particle(i).getvelocity(1)          // vy
            << setw(16) << _particle(i).getvelocity(2) << endl; // vz
    }
  } else cerr << "PROBLEM: Unable to open velocities.dat" << endl;
  coutf.close();
  return;
}

void System :: read_configuration(const string path){ // Read configuration from a .xyz file in directory ../OUTPUT/CONFIG/
  ifstream cinf;
  if(_sim_type < 2) cinf.open(path + "INPUT/CONFIG/config.xyz");
  else cinf.open(path + "INPUT/CONFIG/config.ising");
  if(cinf.is_open()){
    string comment;
    string particle;
    double x, y, z;
    int ncoord;
    cinf >> ncoord;
    if (ncoord != _npart){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE);
    }
    cinf >> comment;
    for(int i=0; i<_npart; i++){
      cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
      _particle(i).setposition(0, this->pbc(_side(0)*x, 0));
      _particle(i).setposition(1, this->pbc(_side(1)*y, 1));
      _particle(i).setposition(2, this->pbc(_side(2)*z, 2));
      _particle(i).acceptmove(); // _x_old = _x_new
    }
  } else cerr << "PROBLEM: Unable to open INPUT file config.xyz"<< endl;
  cinf.close();
  if(_restart and _sim_type > 1){
    int spin;
    cinf.open(path + "INPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++){
      cinf >> spin;
      _particle(i).setspin(spin);
    }
    cinf.close();
  }
  return;
}

void System :: block_reset(int blk, const string path){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open(path + "OUTPUT/output.csv",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av.zeros();
  return;
}

void System :: measure(){ // Measure properties
  _measurement.zeros();
  // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
  int bin;
  vec distance;
  distance.resize(_ndim);
  double penergy_temp=0.0, dr; // temporary accumulator for potential energy
  double kenergy_temp=0.0; // temporary accumulator for kinetic energy
  double tenergy_temp=0.0;
  double pressure_temp=0.0;
  double magnetization=0.0;
  //double virial=0.0;
  if (_measure_penergy or _measure_tenergy or _measure_pressure or _measure_gofr) {
    for (int i=0; i<_npart-1; i++){
      for (int j=i+1; j<_npart; j++){
        distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
        distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
        distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
        dr = sqrt( dot(distance,distance) );
        if (_measure_gofr){
          bin = int(dr/_bin_size);
          if(bin < _n_bins){
            _measurement(_index_gofr+bin) += 2.0;
          }
        }
        if(dr < _r_cut){
          if (_measure_penergy)   penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6); // POTENTIAL ENERGY
          if (_measure_pressure) pressure_temp += 1.0/pow(dr,12) - 0.5/pow(dr,6); // PRESSURE
        }
      }
    }
  }
  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy or _measure_tenergy or _measure_cv){
    penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
    if (_measure_penergy) _measurement(_index_penergy) = penergy_temp;
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy or _measure_tenergy or _measure_temp or _measure_pressure or _measure_cv){
    for (int i=0; i<_npart; i++) kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() ); 
    kenergy_temp /= double(_npart);
    if (_measure_kenergy) _measurement(_index_kenergy) = kenergy_temp;
  }
  // TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
  if (_measure_tenergy or _measure_cv){
    if (_sim_type < 2) tenergy_temp = kenergy_temp + penergy_temp;
    else {
      double s_i, s_j;
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
      }
      tenergy_temp /= double(_npart);
    }
    if (_measure_tenergy) _measurement(_index_tenergy) = tenergy_temp; 
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp) _measurement(_index_temp) = (2.0/3.0) * kenergy_temp;
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure){
    if (!_measure_temp){
      cerr << "PROBLEM: pressure needs to measure also the temperature!" << endl;
      exit(EXIT_FAILURE);
    }
    _measurement(_index_pressure) = _ptail + _rho * _measurement(_index_temp) + 48.0 * pressure_temp / (3.0 * _volume * double(_npart));
  }  
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  if (_measure_magnet or _measure_chi){
    for (int i=0; i<_npart; i++) magnetization += double(_particle(i).getspin());
    magnetization /= double(_npart);
    if (_measure_magnet) _measurement(_index_magnet) = magnetization;
  }
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
  if (_measure_cv) _measurement(_index_cv) = pow(tenergy_temp, 2);
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  if (_measure_chi) _measurement(_index_chi) = _beta * pow(magnetization, 2) * double(_npart);

  _block_av += _measurement; //Update block accumulators

  return;
}

void System :: averages(int blk, const string path){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  _average = _block_av / double(_nsteps);
  if (_measure_cv) _average(_index_cv) = _beta * _beta * (_average(_index_cv) - pow(_average(_index_tenergy), 2)) * double(_npart); // Measuring specific heat
  if (_measure_gofr){ // Normalising the radial distribution function
      for(int i=0; i<_n_bins; i++){
        double dV = (4.0/3.0) * M_PI * (pow(static_cast<double>(i+1) * _bin_size, 3) - pow(static_cast<double>(i) * _bin_size, 3)); // Volume shell
        _average(_index_gofr + i) /= ( static_cast<double>(_npart) * _rho * dV);
    }
  }
  _global_av += _average;
  _global_av2 += _average % _average; // % -> element-wise multiplication

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    coutf.open(path + "OUTPUT/potential_energy.csv",ios::app);
    average  = _average(_index_penergy);
    sum_average = _global_av(_index_penergy);
    sum_ave2 = _global_av2(_index_penergy);
    coutf << blk 
          << "," << average
          << "," << sum_average/double(blk)
          << "," << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    coutf.open(path + "OUTPUT/kinetic_energy.csv",ios::app);
    average  = _average(_index_kenergy);
    sum_average = _global_av(_index_kenergy);
    sum_ave2 = _global_av2(_index_kenergy);
    coutf << blk
          << "," << average
          << "," << sum_average/double(blk)
          << "," << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // TOTAL ENERGY //////////////////////////////////////////////////////////////
  if (_measure_tenergy){
    coutf.open(path + "OUTPUT/total_energy.csv",ios::app);
    average  = _average(_index_tenergy);
    sum_average = _global_av(_index_tenergy);
    sum_ave2 = _global_av2(_index_tenergy);
    coutf << blk
          << "," << average
          << "," << sum_average/double(blk)
          << "," << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp){
    coutf.open(path + "OUTPUT/temperature.csv",ios::app);
    average  = _average(_index_temp);
    sum_average = _global_av(_index_temp);
    sum_ave2 = _global_av2(_index_temp);
    coutf << blk
          << "," << average
          << "," << sum_average/double(blk)
          << "," << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure){
    coutf.open(path + "OUTPUT/pressure.csv",ios::app);
    average  = _average(_index_pressure);
    sum_average = _global_av(_index_pressure);
    sum_ave2 = _global_av2(_index_pressure);
    coutf << blk
          << "," << average
          << "," << sum_average/double(blk)
          << "," << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // GOFR //////////////////////////////////////////////////////////////////////
  if (_measure_gofr) {
    coutf.open(path + "OUTPUT/gofr_final.csv",ios::app);
    ofstream coutf2(path + "OUTPUT/gofr_blocks.csv",ios::app);
    double sum = 0.0, sum2 = 0.0;

    for (int i=0 ; i<_n_bins ; i++){
      sum += _average(_index_gofr + i);
      sum2 += pow(_average(_index_gofr + i), 2);

      if (blk == _nblocks){
        average = _average(_index_gofr + i);
        sum_average = _global_av(_index_gofr + i);
        sum_ave2 = _global_av2(_index_gofr + i);
        coutf << i+1
              << "," << double(i)*_bin_size
              << "," << sum_average/double(blk)
              << "," << this->error(sum_average, sum_ave2, blk) << endl;
      }
    }

    coutf2 << blk
          << "," << sum / double(_n_bins)
          << "," << this->error(sum, sum2, _n_bins) << endl;
    coutf.close();
    coutf2.close();
  }
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  if (_measure_magnet){
    coutf.open(path + "OUTPUT/magnetization.csv",ios::app);
    average  = _average(_index_magnet);
    sum_average = _global_av(_index_magnet);
    sum_ave2 = _global_av2(_index_magnet);
    coutf << blk
          << "," << average
          << "," << sum_average/double(blk)
          << "," << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
  if (_measure_cv){
    coutf.open(path + "OUTPUT/specific_heat.csv",ios::app);
    average  = _average(_index_cv);
    sum_average = _global_av(_index_cv);
    sum_ave2 = _global_av2(_index_cv);
    coutf << blk
          << "," << average
          << "," << sum_average/double(blk)
          << "," << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  if (_measure_chi){
    coutf.open(path + "OUTPUT/susceptibility.csv",ios::app);
    average  = _average(_index_chi);
    sum_average = _global_av(_index_chi);
    sum_ave2 = _global_av2(_index_chi);
    coutf << blk
          << "," << average
          << "," << sum_average/double(blk)
          << "," << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // ACCEPTANCE ////////////////////////////////////////////////////////////////
  double fraction;
  coutf.open(path + "OUTPUT/acceptance.csv",ios::app);
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; 
  coutf << blk << "," << fraction << endl;
  coutf.close();
  
  return;
}

void System :: instantaneous(int extr, const string path){

  ofstream coutf;
  if(_measure_cv) _measurement(_index_cv) = _beta * _beta * (_measurement(_index_cv) - pow(_measurement(_index_tenergy), 2)) * double(_npart); // Measuring specific heat

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    coutf.open(path + "OUTPUT/potential_energy.csv",ios::app);
    coutf << extr << "," << _measurement(_index_penergy) << endl;
    coutf.close();
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    coutf.open(path + "OUTPUT/kinetic_energy.csv",ios::app);
    coutf << extr << "," << _measurement(_index_kenergy) << endl;
    coutf.close();
  }
  // TOTAL ENERGY //////////////////////////////////////////////////////////////
  if (_measure_tenergy){
    coutf.open(path + "OUTPUT/total_energy.csv",ios::app);
    coutf << extr << "," << _measurement(_index_tenergy) << endl;
    coutf.close();
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp){
    coutf.open(path + "OUTPUT/temperature.csv",ios::app);
    coutf << extr << "," << _measurement(_index_temp) << endl;
    coutf.close();
  }
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure){
    coutf.open(path + "OUTPUT/pressure.csv",ios::app);
    coutf << extr << "," << _measurement(_index_pressure) << endl;
    coutf.close();
  }
  // GOFR //////////////////////////////////////////////////////////////////////
  if (_measure_gofr) {
    coutf.open(path + "OUTPUT/gofr_final.csv",ios::app);
    for(int i=0; i<_n_bins; i++){ // Normalising the radial distribution function
      double dV = (4.0/3.0) * M_PI * (pow(static_cast<double>(i+1) * _bin_size, 3) - pow(static_cast<double>(i) * _bin_size, 3)); // Volume shell
      _measurement(_index_gofr + i) /= ( static_cast<double>(_npart) * _rho * dV);
      coutf << i+1 << "," << double(i)*_bin_size << "," << _measurement(_index_gofr + i) << endl;
    }
    coutf.close();
  }
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  if (_measure_magnet){
    coutf.open(path + "OUTPUT/magnetization.csv",ios::app);
    coutf << extr << "," << _measurement(_index_magnet) << endl;
    coutf.close();
  }
  // SPECIFIC HEAT ///////////////////////////////////////////////////////////// TO BE FIXED
  if (_measure_cv){
    coutf.open(path + "OUTPUT/specific_heat.csv",ios::app);
    coutf << extr << "," << _measurement(_index_cv) << endl;
    coutf.close();
  }
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  if (_measure_chi){
    coutf.open(path + "OUTPUT/susceptibility.csv",ios::app);
    coutf << extr << "," << _measurement(_index_chi) << endl;
    coutf.close();
  }
  double fraction;
  // ACCEPTANCE ////////////////////////////////////////////////////////////////
  coutf.open(path + "OUTPUT/acceptance.csv",ios::app);
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; 
  coutf << extr << "," << fraction << endl;
  coutf.close();
  
  return;
}

void System :: equilibration(const string path){ // Perform the equilibration of the system

  fmt::print("\nEQUILIBRATION\n");

  if (_sim_type == 0){
    if (!_measure_temp){
      cerr << "PROBLEM: equilibration needs to measure the temperature!" << endl;
      exit(EXIT_FAILURE);
    }

    double temp_target = _temp;
    int div = 1;
    double temp, corr;
    double val = 0.0, sigma = 0.0;
    ofstream coutf(path + "OUTPUT/equilibration_" + format("{:.3f}", temp_target) + ".csv");
    coutf << "Equilibration started!" << endl << endl;
      

    do {
      if(temp_target/double(div) < sigma) corr = sigma;
      else corr = temp_target/double(div);

      if(val != 0.0 and val < temp_target) _temp = temp + corr;
      else if(val != 0.0 and val > temp_target) _temp = temp - corr;
      temp = _temp;
      div *= 2;

      this->read_configuration(path);
      this->initialize_velocities(path);
      this->initialize_properties(path);
      this->block_reset(0);

      fmt::print("\nInitial temperature: {:.3f}\n", _temp);
      coutf << "Initial temperature: " << _temp << endl;
      
      for(int i=0; i<_nblocks; i++) {
        for(int j=0; j<_nsteps; j++) {
          Progress_Bar(i*_nsteps + j, _nblocks*_nsteps -1);
          this->step();
          this->measure();
        }
        this->averages(i+1, path);
        this->block_reset(i+1, path);
      }

      val = _global_av(_index_temp)/double(_nblocks);
      sigma = this->error(_global_av(_index_temp), _global_av2(_index_temp), _nblocks);


      fmt::print("\nMeasured T*: {:.5f} ± {:.5f}\tTarget T*: {:.3f}\n", val, sigma, temp_target); 
      coutf << "\tMeasured T*: " << val << " ± " << sigma << "\tTarget T*: " << temp_target << endl;

    } while(temp_target < val - 3.0*sigma or temp_target > val + 3.0*sigma);

    fmt::print("\nEQUILIBRATION COMPLETED\n");
    coutf << "\nEquilibration completed!\n\nT_INPUT: " << _temp << endl << endl << endl;
    coutf.close();

  } else {
    if(!_measure_tenergy) {
      cerr << "PROBLEM: equilibration needs to measure the internal energy!" << endl;
      exit(EXIT_FAILURE);
    }

    for(int i=0; i<_nblocks; i++) {
      for(int j=0; j<_nsteps; j++) {
        Progress_Bar(i*_nsteps + j, _nblocks*_nsteps -1);
        this->step();
        this->measure();
      }
      this->averages(i+1, path);  
      this->block_reset(i+1, path);
    }

  }

  this->rename_files(path, path + "OUTPUT/EQUILIBRATION/");
  fmt::print("\n");
  this->initialize_properties(path);
  this->block_reset(0);
  
}

void System :: rename_files(const string oldpath, const string newpath){

  ostringstream temp_stream;
  temp_stream << fixed << setprecision(3) << _temp;
  string temp_str = temp_stream.str();
  ostringstream h_stream;
  h_stream << fixed << setprecision(2) << _H;
  string h_str = h_stream.str();

  string oldname, newname;
  const char* old_name = oldname.c_str();
  const char* new_name = newname.c_str();

  if(_measure_temp) {
    oldname = oldpath + "OUTPUT/temperature.csv";
    newname = newpath + "temperature_H=" + h_str + "_t=" + temp_str + ".csv";
    old_name = oldname.c_str();
    new_name = newname.c_str();
    rename(old_name, new_name);
  }
  if(_measure_penergy) {
    oldname = oldpath + "OUTPUT/potential_energy.csv";
    newname = newpath + "potential_energy_H=" + h_str + "_t=" + temp_str + ".csv";
    old_name = oldname.c_str();
    new_name = newname.c_str();
    rename(old_name, new_name);
  }
  if(_measure_kenergy) {
    oldname = oldpath + "OUTPUT/kinetic_energy.csv";
    newname = newpath + "kinetic_energy_H=" + h_str + "_t=" + temp_str + ".csv";
    old_name = oldname.c_str();
    new_name = newname.c_str();
    rename(old_name, new_name);
  }
  if(_measure_tenergy) {
    oldname = oldpath + "OUTPUT/total_energy.csv";
    newname = newpath + "total_energy_H=" + h_str + "_t=" + temp_str + ".csv";
    old_name = oldname.c_str();
    new_name = newname.c_str();
    rename(old_name, new_name);
  }
  if(_measure_pressure) {
    oldname = oldpath + "OUTPUT/pressure.csv";
    newname = newpath + "pressure_H=" + h_str + "_t=" + temp_str + ".csv";
    old_name = oldname.c_str();
    new_name = newname.c_str();
    rename(old_name, new_name);
  }
  if(_measure_gofr) {
    oldname = oldpath + "OUTPUT/gofr_blocks.csv";
    newname = newpath + "gofr_blocks_H=" + h_str + "_t=" + temp_str + ".csv";
    old_name = oldname.c_str();
    new_name = newname.c_str();
    rename(old_name, new_name);

    oldname = oldpath + "OUTPUT/gofr_final.csv";
    newname = newpath + "gofr_final_H=" + h_str + "_t=" + temp_str + ".csv";
    old_name = oldname.c_str();
    new_name = newname.c_str();
    rename(old_name, new_name);
  }
  if(_measure_magnet) {
    oldname = oldpath + "OUTPUT/magnetization.csv";
    newname = newpath + "magnetization_H=" + h_str + "_t=" + temp_str + ".csv";
    old_name = oldname.c_str();
    new_name = newname.c_str();
    rename(old_name, new_name);
  }
  if(_measure_cv) {
    oldname = oldpath + "OUTPUT/specific_heat.csv";
    newname = newpath + "specific_heat_H=" + h_str + "_t=" + temp_str + ".csv";
    old_name = oldname.c_str();
    new_name = newname.c_str();
    rename(old_name, new_name);
  }
  if(_measure_chi) {
    oldname = oldpath + "OUTPUT/susceptibility.csv";
    newname = newpath + "susceptibility_H=" + h_str + "_t=" + temp_str + ".csv";
    old_name = oldname.c_str();
    new_name = newname.c_str();
    rename(old_name, new_name);
  }

}

double System :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

int System :: get_nbl(){
  return _nblocks;
}

int System :: get_nsteps(){
  return _nsteps;
}

double System :: get_vtail(){
  return _vtail;
}

double System :: get_ptail(){
  return _ptail;
}

int System :: get_sim_type(){
  return _sim_type;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
