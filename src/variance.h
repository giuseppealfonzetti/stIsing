#ifndef variance_H
#define variance_H

#include "binarynodeClass.h"

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd sampleH(
        Eigen::Map<Eigen::VectorXd> THETA,
        Eigen::Map<Eigen::MatrixXd> DATA,
        std::vector<bool> &CONSTRAINTS,
        const bool INVERTFLAG = false,
        const bool VERBOSEFLAG = false
){
    // Identify dimensions
    const unsigned int n = DATA.rows();
    const unsigned int d = THETA.size();
    const unsigned int p = DATA.cols();
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(d,d);

    binarynodeClass Node;
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
    for(unsigned int i = 0; i < n; i++){

        Eigen::VectorXd data_i = DATA.row(i);
        for(unsigned int node = 0; node < p; node++){

            Node.setup_(data_i, THETA, CONSTRAINTS, p, node);
            Eigen::VectorXd node_gradient = Node.gradient_();
            out += node_gradient * node_gradient.transpose();

        }


    }

    out /= n;

    if(INVERTFLAG){
        const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(d,d);
        Eigen::LLT<Eigen::MatrixXd> llt;
        llt.compute(out);
        out = llt.solve(Id);
    }
    return out;

}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd sampleJ(
        Eigen::Map<Eigen::VectorXd> THETA,
        Eigen::Map<Eigen::MatrixXd> DATA,
        std::vector<bool> &CONSTRAINTS,
        const bool VERBOSEFLAG = false
){
    // Identify dimensions
    const unsigned int n = DATA.rows();
    const unsigned int d = THETA.size();
    const unsigned int p = DATA.cols();
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(d,d);

    binarynodeClass Node;

    for(unsigned int i = 0; i < n; i++){

        Eigen::VectorXd unit_gradient = Eigen::VectorXd::Zero(d);
        Eigen::VectorXd data_i = DATA.row(i);
        for(unsigned int node = 0; node < p; node++){

            Node.setup_(data_i, THETA, CONSTRAINTS, p, node);
            unit_gradient += Node.gradient_();
        }

        J += unit_gradient * unit_gradient.transpose();



    }

    J /= n;


    return J;

}


//' @export
// [[Rcpp::export]]
Rcpp::List sampleVar(
        Eigen::Map<Eigen::VectorXd> THETA,
        Eigen::Map<Eigen::MatrixXd> DATA,
        std::vector<bool> &CONSTRAINTS,
        const unsigned int NU,
        const unsigned int METHOD,
        const unsigned int RANGE,
        const bool TOTFLAG,
        const bool PRINTFLAG
){

    // Identify dimensions
    const unsigned int d = THETA.size();
    const unsigned int n = DATA.rows();
    const unsigned int kk= DATA.cols();

    // Initialise variance matrices and std errors vectors
    Eigen::MatrixXd tot, jmat, sandwich, cond;
    Eigen::VectorXd stochastic_se, statistical_se, tot_se;


    Rcpp::List var;
    Rcpp::List se;

    // Compute inverse of negative hessian
    Eigen::MatrixXd inv_hmat = sampleH(THETA, DATA, CONSTRAINTS, true, PRINTFLAG);

    if(METHOD == 0){
        jmat = sampleJ(THETA, DATA, CONSTRAINTS, PRINTFLAG);
        sandwich = inv_hmat*jmat*inv_hmat;
        tot = sandwich/n;

        var = Rcpp::List::create( Rcpp::Named("var_tot") =  sandwich/static_cast<double>(n));
        se  = Rcpp::List::create( Rcpp::Named("se_tot" ) = (sandwich/static_cast<double>(n)).diagonal().array().sqrt() );

    }

    if(METHOD == 1){
        jmat = sampleJ(THETA, DATA, CONSTRAINTS, PRINTFLAG);
        sandwich = inv_hmat*jmat*inv_hmat;
        double st_scale = NU*RANGE;

        if(TOTFLAG) {

            var = Rcpp::List::create(
                Rcpp::Named("var_stoc") = sandwich/st_scale,
                Rcpp::Named("var_stat") = sandwich/n,
                Rcpp::Named("var_tot")  = sandwich.array()*( (1/st_scale) + (1/static_cast<double>(n)) )
            );

            se = Rcpp::List::create(
                Rcpp::Named("se_stoc") = (sandwich/st_scale).diagonal().array().sqrt(),
                Rcpp::Named("se_stat") = (sandwich/n).diagonal().array().sqrt(),
                Rcpp::Named("se_tot")  = (sandwich*( (1/st_scale) + (1/static_cast<double>(n)) )).diagonal().array().sqrt()
            );
        }else{
            var = Rcpp::List::create(Rcpp::Named("var_stoc") = sandwich/st_scale);
            se  = Rcpp::List::create(Rcpp::Named("se_stoc") = (sandwich/st_scale).diagonal().array().sqrt());
        }
    }

    if(METHOD == 2){

        double st_scale = NU*RANGE;

        if(TOTFLAG) {

            jmat = sampleJ(THETA, DATA, CONSTRAINTS, PRINTFLAG);
            sandwich = inv_hmat*jmat*inv_hmat;
            var = Rcpp::List::create(
                Rcpp::Named("var_stoc") = inv_hmat/st_scale,
                Rcpp::Named("var_stat") = sandwich/static_cast<double>(n),
                Rcpp::Named("var_tot")  = (sandwich/static_cast<double>(n))+(inv_hmat/st_scale)
            );

            se = Rcpp::List::create(
                Rcpp::Named("se_stoc") = (inv_hmat/st_scale).diagonal().array().sqrt(),
                Rcpp::Named("se_stat") = (sandwich/static_cast<double>(n)).diagonal().array().sqrt(),
                Rcpp::Named("se_tot")  = ((sandwich/static_cast<double>(n))+(inv_hmat/st_scale)).diagonal().array().sqrt()
            );
        }else{
            var = Rcpp::List::create(Rcpp::Named("var_stoc") = inv_hmat/st_scale);
            se  = Rcpp::List::create(Rcpp::Named("se_stoc") = (inv_hmat/st_scale).diagonal().array().sqrt());
        }
    }

    if(METHOD == 3){
        double st_scale = NU*RANGE;
        jmat = sampleJ(THETA, DATA, CONSTRAINTS, PRINTFLAG);
        double alpha_1 = static_cast<double>(n-1)/static_cast<double>(n) + (static_cast<double>(n-1)/static_cast<double>(n*kk-1))/(static_cast<double>(n));
        double alpha_2 = -(static_cast<double>(n-1)/static_cast<double>(n*kk-1))/(static_cast<double>(n));

        if(TOTFLAG) {

            jmat = sampleJ(THETA, DATA, CONSTRAINTS, PRINTFLAG);
            sandwich = inv_hmat*jmat*inv_hmat;
            var = Rcpp::List::create(
                Rcpp::Named("var_stoc") = (alpha_1*inv_hmat + alpha_2*sandwich)/st_scale,
                Rcpp::Named("var_stat") = sandwich/static_cast<double>(n),
                Rcpp::Named("var_tot")  = (sandwich/static_cast<double>(n))+((alpha_1*inv_hmat + alpha_2*sandwich)/st_scale)
            );

            se = Rcpp::List::create(
                Rcpp::Named("se_stoc") = ((alpha_1*inv_hmat + alpha_2*sandwich)/st_scale).diagonal().array().sqrt(),
                Rcpp::Named("se_stat") = (sandwich/static_cast<double>(n)).diagonal().array().sqrt(),
                Rcpp::Named("se_tot")  = ((sandwich/static_cast<double>(n))+((alpha_1*inv_hmat + alpha_2*sandwich)/st_scale)).diagonal().array().sqrt()
            );
        }else{
            var = Rcpp::List::create(Rcpp::Named("var_stoc") = (alpha_1*inv_hmat + alpha_2*sandwich)/st_scale);
            se  = Rcpp::List::create(Rcpp::Named("se_stoc") = ((alpha_1*inv_hmat + alpha_2*sandwich)/st_scale).diagonal().array().sqrt());
        }
    }


    Rcpp::List out = Rcpp::List::create(Rcpp::Named("var") = var, Rcpp::Named("se") = se);
    return out;
}
#endif
