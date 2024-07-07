#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "../../Libraries/fmtlib.hpp"
#include "../../Libraries/RandomGen/random.hpp"
#include "../../Libraries/library.hpp"
#include "../../Libraries/blockaverage.hpp"

using namespace std;



//double Errore(vector <double>, vector <double>, int);



int main(int argc, char** argv){

    // initial check
    fmt::print("\n");
    if (argc != 3) {
        cout << "ERROR! Program usage: " << argv[0]  << " <n_extractions> <n_blocks>" << endl << endl;
        return -1;
    }

    // setting parameters
    unsigned int m = stod(argv[1]);      // number of extractions
    unsigned int n = stod(argv[2]);      // number of blocks
    unsigned int l = (int)m/n;           // number of extraction per block
    if (m%n != 0) {
        cout << "ERROR! <n_extractions> must be a mutiple of <n_blocks>" << endl << endl;
        return -1;
    }

    // initializing random numbers generator
    Random rnd("../../Libraries/RandomGen/");



    // 1ST METHOD: one big for cycle
    /*–––––––––––––––––––––––––––––––––– PART 1, 2, 3 ––––––––––––––––––––––––––––––––––*/

    double mean, mean2;
    double sigma, sigma2;

    double sum = 0, sum2 = 0;
    double sum_s = 0, sum_s2 = 0;
    double err_mean = 0;

    double prog_mean = 0, prog_mean2 = 0;
    double prog_sigma = 0, prog_sigma2 = 0;
    double err_sigma = 0;

    double chi2 = 0;

    ofstream out("output.csv");
    out << "# extr,#blocks,mean-0.5,error,variance,error,chi squared" << endl;

    for (int i=0 ; i<n ; i++){
        // computing blocking average, variance (and squares) and observed values for the chi squared
        mean = 0, sigma = 0, chi2 = 0;
        vector <int> hit(n, 0);
        for (int j=0 ; j<l ; j++){
            double r = rnd.Rannyu();
            mean += r;
            sigma += pow((r-0.5), 2);

            for (int c=0 ; c<n ; c++) { 
                if (r >= (double)c/n && r < (double)(c+1)/n) {hit[c]++;} 
            }
        }
        mean /= l;
        sigma /= l;
        mean2 = mean*mean;
        sigma2 = sigma*sigma;

        // computing progressive mean
        sum += mean;
        sum2 += mean2;
        prog_mean = (double)sum/(i+1);
        prog_mean2 = (double)sum2/(i+1);
        err_mean = Error(prog_mean, prog_mean2, i);

        // computing progressive variance
        sum_s += sigma;
        sum_s2 += sigma2;
        prog_sigma = (double)sum_s/(i+1);
        prog_sigma2 = (double)sum_s2/(i+1);
        err_sigma = Error(prog_sigma, prog_sigma2, i);

        //computing chi squared value
        for (int j=0 ; j<n ; j++) {chi2 += (pow(((double)hit[j] - l/n), 2)/(double)(l/n));}

        // output file: number of extractions, of blocks mean values, errors and chi squared
        out << (i+1)*l << "," << i+1 << "," << prog_mean-0.5 << "," << err_mean << "," << prog_sigma-(1./12) << "," << err_sigma << "," << chi2 << endl;
    }
    out.close();



    // 2ND METHOD: using classes
    /*–––––––––––––––––––––––––––––––––– PART 1 ––––––––––––––––––––––––––––––––––*/

    fmt::print("Computing mean\n");
    BA_Mean mea(m,n);
    out.open("mean.csv");
    out << "blocks,extractions,mean,error" << endl;
    mea.Progressive(out);
    out.close();
    fmt::print("\n\n");



    /*–––––––––––––––––––––––––––––––––– PART 2 ––––––––––––––––––––––––––––––––––*/

    fmt::print("Computing variance\n");
    BA_Variance var(m,n);
    out.open("variance.csv");
    out << "blocks,extractions,mean,error" << endl;
    var.Progressive(out);
    out.close();
    fmt::print("\n\n");



    return 0;
}




// computing the progressive mean with vectors
/*
int main(int argc, char** argv){

    // initial check
    cout << endl;
    if (argc != 3) {
        cout << "ERROR! Program usage: " << argv[0]  << " <n_extractions> <n_batches>" << endl << endl;
        return -1;
    }

    int m = atoi(argv[1]);      // number of extractions
    int n = atoi(argv[2]);      // number of batches
    int l = (int)m/n;           // nummber of extraction per batch
    if (m%n != 0) {
        cout << "ERROR! <n_extractions> must be a mutiple of <n_batches>" << endl << endl;
        return -1;
    }

    // random numbers generator
    Random rnd("../../Libraries/RandomGen/");
    // calcolo delle medie dei blocchi e dei loro quadrati
    vector<double> media;
    vector<double> media2;
    for (int i=0 ; i<n ; i++){
        double somma = 0;

        for (int j=0 ; j<l ; j++){
            somma += rnd.Rannyu();
        }
        media.push_back(somma/l);
        media2.push_back(media[i]*media[i]);
        
        //cout << media[i] << " " << media2[i] << endl;
    }

    // calcolo delle medie progressive
    vector<double> media_vera;
    vector<double> media2_vera;
    vector<double> errore;
    for (int i=0 ; i<n ; i++){
        double somma = 0;
        double somma2 = 0;

        for (int j=0 ; j<i+1 ; j++){
            somma += media[j];
            somma2 += media2[j];
        }

        media_vera.push_back((double)somma/(i+1));
        media2_vera.push_back((double)somma2/(i+1));
        errore.push_back(Errore(media_vera, media2_vera, i));

    }

    //stampo su file i numeri di lanci per ogni blocco, le medie e gli errori
    ofstream fout("media.csv");

    fout << "n_lanci\tmedia-0.5\terrore" << endl;
    for (int i=0 ; i<n ; i++){
        fout << (i+1)*l << "\t" << media_vera[i]-0.5 << "\t" << errore[i] << endl;
    }


    return 0;
}

// errore: deviazione standard della media: sqrt[(<x^2> - <x>^2) / (n-1)]
double Errore(vector <double> media, vector <double> media2, int n){
    if (n == 0) {return 0;}
    else {return sqrt((media2[n] - pow(media[n],2)) / n);}
}
*/

