#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <RcppClock.h>
#include <math.h>
#include <random>
#include "utils.h"
#include "variance.h"
#include "binarynodeClass.h"

// [[Rcpp::depends(RcppEigen, RcppClock)]]

//' @export
// [[Rcpp::export]]
Rcpp::List ncl(
        Eigen::MatrixXd &DATA,
        Eigen::VectorXd &THETA,
        std::vector<bool> &CONSTRAINTS,
        const bool VERBOSEFLAG = false
){
    // Identify dimensions
    const unsigned int n = DATA.rows();
    const unsigned int d = THETA.size();
    const unsigned int p = (- 1 + sqrt(1 + 8*d))/2;

    binarynodeClass Node;
    double cl = 0;
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
    for(unsigned int i = 0; i < n; i++){

        Eigen::VectorXd data_i = DATA.row(i);
        for(unsigned int node = 0; node < p; node++){

            Node.setup_(data_i, THETA, CONSTRAINTS, p, node);
            cl -= Node.ll_();
            gradient -= Node.gradient_();
        }
    }

    Rcpp::List output =
        Rcpp::List::create(
            Rcpp::Named("nll") = cl,
            Rcpp::Named("ngradient") = gradient
        );

    return output;

}


//' @export
// [[Rcpp::export]]
Rcpp::List isingGraph(
        const Eigen::MatrixXd &DATA,
        const Eigen::VectorXd &THETA_INIT,
        const std::vector<bool> &CONSTRAINTS,
        const unsigned int MAXT,
        const unsigned int BURN,
        const double STEPSIZE,
        const double NU,
        const int METHODFLAG,
        const unsigned int SEED = 123,
        const bool VERBOSEFLAG = false

){
    // Set up clock monitor to export to R session trough RcppClock
    Rcpp::Clock clock;
    clock.tick("main");

    // Identify dimensions
    const unsigned int n = DATA.rows();
    const unsigned int d = THETA_INIT.size();
    const unsigned int p = (- 1 + sqrt(1 + 8*d))/2;
    const unsigned int kk= p;

    // scaling constant according to sampling_scheme
    // Compute scaling constant
    double scale;
    switch(METHODFLAG){
    case 0:
        scale = 1/static_cast<double>(n) ;
        break;
    case 1:
        scale = 1/static_cast<double>(NU);
        break;
    case 2:
        scale = 1/static_cast<double>(NU);
        break;
    case 3:
        scale = 1/static_cast<double>(NU);
        break;
    }
    // Rcpp::Rcout << " Final scale = " << scale<< "\n";

    // Initialize storage for iterations quantities
    Eigen::MatrixXd path_theta    = Eigen::MatrixXd::Zero(MAXT + 1, d); path_theta.row(0)    = THETA_INIT;
    Eigen::MatrixXd path_av_theta = Eigen::MatrixXd::Zero(MAXT + 1, d); path_av_theta.row(0) = THETA_INIT;
    Eigen::MatrixXd path_grad     = Eigen::MatrixXd::Zero(MAXT,     d);
    // Eigen::VectorXd path_ncl      = Eigen::VectorXd::Zero(MAXT);

    // Initialize vector of indexes for entries in pairs_table
    // Set-up the randomizer
    std::vector<int> sublik_pool(n*p) ;
    std::iota (std::begin(sublik_pool), std::end(sublik_pool), 0);
    std::vector<int> vector_weights(n*p);


    ///////////////////////
    /* OPTIMISATION LOOP */
    ///////////////////////

    binarynodeClass Node;
    Eigen::VectorXd theta_t = THETA_INIT;
    for(unsigned int t = 1; t <= MAXT; t++){
        // check user interruption
        Rcpp::checkUserInterrupt();
        Rcpp::Rcout << "\r Iter:" << t << " ";

        /////////////////////
        // SAMPLING SCHEME //
        /////////////////////
        clock.tick("sampling_step");
        Rcpp::NumericMatrix sampling_weights(n,kk);
        std::fill(sampling_weights.begin(), sampling_weights.end(), 0) ;

        double prob;
        switch(METHODFLAG){
        case 0:
            std::fill( sampling_weights.begin(), sampling_weights.end(), 1);
            break;
        case 1:
            prob = 1/static_cast<double>(n);
            sampling_weights = rmultinom_wrapper(prob, n, NU, kk);
            break;
        case 2:
            prob = static_cast<double>(NU)/static_cast<double>(n);
            for(unsigned int i = 0; i < n; i++){
                for(unsigned int k = 0; k < kk; k++){
                    if(R::runif(0,1) < prob ) sampling_weights(i, k) = 1;
                }
            }
            break;
        case 3:
            std::mt19937 randomizer(SEED+t);
            std::fill(vector_weights.begin(), vector_weights.end(), 0) ;
            std::shuffle(sublik_pool.begin(), sublik_pool.end(), randomizer);
            for(unsigned int draw = 0; draw < int(NU*kk); draw ++){
                vector_weights[sublik_pool[draw]] = 1;
                // int row = sublik_pool[draw] / kk;
                // int col = sublik_pool[draw] % kk;
                // sampling_weights(row, col) = 1;
            }

            std::copy(vector_weights.begin(), vector_weights.end(), sampling_weights.begin());
            break;
        }
        clock.tock("sampling_step");

        //initialize iteration quantities
        Eigen::VectorXd ngradient_t = Eigen::VectorXd::Zero(d);
        double ncl = 0;

        ///////////////////////////
        /* GRADIENT COMPUTATION  */
        ///////////////////////////
        // looping over units
        clock.tick("stochastic_gradient");
        for(unsigned int i = 0; i < n; i++){
            Eigen::VectorXd data_i = DATA.row(i);

            // looping over nodes
            for(unsigned int node = 0; node < kk; node++){
                unsigned int weight = sampling_weights(i, node);
                if(weight != 0){

                    Node.setup_(data_i, theta_t, CONSTRAINTS, p, node);
                    ngradient_t -= weight * Node.gradient_();
                    // ncl -= Node.ll_();

                }
            }

        }
        clock.tock("stochastic_gradient");

        ncl *= scale;
        ngradient_t *= scale;

        ///////////////////////////
        /*    PARAMETERS UPDATE  */
        ///////////////////////////
        clock.tick("update");
        double stepsize_t = STEPSIZE * pow(t, -.501);
        theta_t -= stepsize_t * ngradient_t;
        clock.tock("update");


        /////////////////////////////////
        /* STORE ITERATION QUANTITIES  */
        /////////////////////////////////
        path_theta.row(t ) = theta_t;
        path_grad.row(t-1) = ngradient_t;
        // path_ncl(t-1)      = ncl;

        // averaging after burnsize
        if(t <= BURN){
            path_av_theta.row(t) = path_theta.row(t);
        }else{
            path_av_theta.row(t) = ( (t - BURN - 1) * path_av_theta.row(t - 1) + path_theta.row(t) ) / (t - BURN);
        }
    }

    clock.tock("main");
    clock.stop("clock");





    // Rcpp::NumericMatrix sampling_weights(n,kk);
    // std::fill(sampling_weights.begin(), sampling_weights.end(), 0) ;
    // std::fill(vector_weights.begin(), vector_weights.end(), 0) ;
    // std::shuffle(sublik_pool.begin(), sublik_pool.end(), randomizer);
    // for(unsigned int draw = 0; draw < p; draw ++){
    //     vector_weights[sublik_pool[draw]] = 1;
    // }
    // std::copy(vector_weights.begin(), vector_weights.end(), sampling_weights.begin());


    // Rcpp::Rcout << "5/3 =(" << 5/3 << "," << 5%3<< ")\n" ;

    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("path_theta") = path_theta,
        Rcpp::Named("path_av_theta") = path_av_theta,
        Rcpp::Named("path_grad") = path_grad,
        // Rcpp::Named("path_ncl") = path_ncl,
        Rcpp::Named("scale") = scale,
        Rcpp::Named("n") = n,
        // Rcpp::Named("weights") = weights,
        // Rcpp::Named("sublik_pool") = sublik_pool,
        // Rcpp::Named("sampling_weights") = sampling_weights,


        Rcpp::Named("methodflag") = METHODFLAG
    );

    return output;

}
