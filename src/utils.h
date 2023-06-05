#ifndef utils_H
#define utils_H

#include <random>
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


//'@export
// [[Rcpp::export]]
std::vector<int> hyper_sampling(const unsigned int K, const unsigned int N, const unsigned int SEED){
    std::mt19937 randomizer(SEED);
    std::vector<int> pool(N*K);
    std::iota (std::begin(pool), std::end(pool), 0);
    std::shuffle(pool.begin(), pool.end(), randomizer);
    return(pool);
}

//'@export
// [[Rcpp::export]]
std::vector<int> unit_sampling(const unsigned int N, const unsigned int SEED){
     std::mt19937 randomizer(SEED);
     std::vector<int> pool(N);
     std::iota (std::begin(pool), std::end(pool), 0);
     std::shuffle(pool.begin(), pool.end(), randomizer);
     return(pool);
}

//'@export
// [[Rcpp::export]]
 std::vector<int> components_given_unit(const unsigned int UNIT, const unsigned int K){
     std::vector<int> pool(K);
     std::iota (std::begin(pool), std::end(pool), UNIT*K);
     return(pool);
 }

//'@export
// [[Rcpp::export]]
std::vector<int> bernoulli_sampling(const unsigned int K, const unsigned int N, const double PROB){
     std::vector<int> pool;
     for( int iterator = 0; iterator < N*K; iterator++){
        if(R::runif(0,1) < PROB ) pool.push_back(iterator);
     }
     return(pool);
}

//'@export
// [[Rcpp::export]]
std::vector<int> index_to_component(const unsigned int P, const unsigned int N, const unsigned int INDEX){
    if(INDEX >= N*P) Rcpp::stop("'INDEX' must be less than N*P (count starts from 0)");
    int i = INDEX / P;
    int p = INDEX % P;
    std::vector<int> component{i,p};
    return(component);
}

#endif
