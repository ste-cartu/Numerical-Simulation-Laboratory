#include "blockaverage.hpp"




/*––––––––––––––––––––––––––––––– BLOCKING AVERAGE –––––––––––––––––––––––––––––––*/

void BlockingAverage :: Progressive(std::ofstream& out) {
    
    double val, val2;
    double sum = 0, sum2 = 0;
    double prog_val, prog_val2;
    double err;

    for (int i=0 ; i<blocks_ ; i++) {
        val = Increase();
        val2 = val*val;

        sum += val;
        sum2 += val2;

        prog_val = (double)sum/(i+1);
        prog_val2 = (double)sum2/(i+1);
        err = Error(prog_val, prog_val2, i);
        out << (i+1) << "," << (i+1)*dim_ << "," << prog_val << "," << err << std::endl;
    }

}


double BA_Mean :: Increase() {

    double val = 0;
    for (int i=0 ; i<dim_ ; i++) {val += rnd_.Rannyu();}
    return val/dim_;

}


double BA_Variance :: Increase() {

    double val = 0;
    for (int i=0 ; i<dim_ ; i++) {val += pow(rnd_.Rannyu() - 0.5, 2);}
    return val/dim_;

}


double BA_Buffon :: Increase() {

    int n_hit = 0, n_tot = 0;
    while (n_tot < dim_) {

        // sampling uniformly y1 in range [0, d)
        double Y1 = rnd_.Rannyu(0, distance_);

        // sampling uniformly x2 and y2 in range [-1,1)
        double x2 = rnd_.Rannyu(-1, 1);
        double y2 = rnd_.Rannyu(-1, 1);

        // slecting values inside the circle
        double d = sqrt(x2*x2 + y2*y2);
        if (d < 1) {
            double Y2 = Y1 + length_*(y2/d);

            if(Y1 == 0 || Y2 <= 0 || Y2 > distance_) {n_hit++;}
            n_tot++;
        }
    }

    double val = (2*length_*n_tot)/(n_hit*distance_);

    return val;

}


double BA_Integral :: Increase() {

    double val = 0;
    for(int i=0 ; i<dim_ ; i++) {
        double u = rnd_.Rannyu();
        double x = Inv_(u);
        val += f_(x)/Sampl_(x);
    }

    return val/dim_;

}


void BA_Option :: SetSampling(double tinit, double tfin, unsigned int nstep) {
    
    start_ = tinit;
    stop_ = tfin;
    nstep_ = nstep;

}


double BA_Option :: Increase() {
    
    double val = 0, price = initial_;
    double dt = static_cast<double>(stop_ - start_)/nstep_;

    for(int i=0 ; i<dim_ ; i++) {
        for(int j=0 ; j<nstep_ ; j++) {
            price *= std::exp( (drift_ - 0.5*std::pow(volatility_,2))*dt + volatility_*rnd_.Gauss(0.,1.)*std::sqrt(dt));
        }
        val += std::exp(-drift_*stop_) * Profit_(price, strike_);
        price = initial_;
    }

    return val/dim_;

}



