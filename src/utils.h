#ifndef utils_H
#define utils_H

/* Function to print vector to rsession in a tidy way*/
void print_vec(Eigen::VectorXd vector){
    Rcpp::Rcout<<"(";
    for(int i = 0; i < vector.size(); i++){
        Rcpp::Rcout << vector(i);
        if(i != vector.size()-1) Rcpp::Rcout<<", ";
    }
    Rcpp::Rcout<<")\n";
}

//'@export
// [[Rcpp::export]]
Rcpp::NumericMatrix rmultinom_wrapper(const double prob, const unsigned int classes, const unsigned int batch, const unsigned int K) {

    Rcpp::NumericVector probs(classes, prob);
    Rcpp::IntegerVector outcome(classes);
    R::rmultinom(batch, probs.begin(), classes, outcome.begin());


    Rcpp::NumericMatrix out(classes, K);
    for(unsigned int j = 0; j < K; j++){
        out(Rcpp::_,j) = outcome;
    }

    return out;
}

#endif
