#ifndef binarynodesClass_H
#define binarynodesClass_H

class binarynodeClass{
    public:
        Eigen::VectorXd _data;
        Eigen::VectorXd _theta;
        std::vector<bool> _constraints;
        unsigned int _n_nodes;
        unsigned int _node_index;

        // intermediate quantities
        double _eta;
        double _d;



        void setup_(Eigen::VectorXd DATA, Eigen::VectorXd THETA, std::vector<bool> CONSTRAINTS, unsigned int N_NODES, unsigned int NODE_INDEX);
        double ll_(const bool VERBOSEFLAG = false);
        Eigen::VectorXd gradient_(const bool VERBOSEFLAG = false);

    private:
        double eta_(const bool VERBOSEFLAG = false);

};

void binarynodeClass::setup_(Eigen::VectorXd DATA, Eigen::VectorXd THETA, std::vector<bool> CONSTRAINTS, unsigned int N_NODES, unsigned int NODE_INDEX){
    _data        = DATA;
    _theta       = THETA;
    _constraints = CONSTRAINTS;
    _n_nodes     = N_NODES;
    _node_index  = NODE_INDEX;
    _eta         = eta_();
    _d           = THETA.size();
}

double binarynodeClass::eta_(const bool VERBOSEFLAG){

    // Initialize eta at node potential parameter
    double eta = _theta(_node_index);

    for(unsigned int j = 0; j < _n_nodes; j++){
        if(_node_index!=j){
            unsigned int beta_index;
            if(j < _node_index){
                if(j == 0){
                    beta_index = _node_index;
                    // Rcpp::Rcout << "j: " << j <<", beta_index:" << beta_index << ", beta:"<< theta(beta_index) << ", data:"<< data(j) <<"\n";

                }else{
                    unsigned int accumulator = 0;
                    for(unsigned int i = 2; i < j+2; i++) accumulator += i;
                    beta_index = _node_index + j * _n_nodes - accumulator;
                    // Rcpp::Rcout << "j: " << j << ", beta_index:" << beta_index << ", beta:"<< theta(beta_index) << ", data:"<< data(j) <<"\n";
                }
            }else
                if ( j > _node_index){
                    unsigned int accumulator = 0;
                    for(unsigned int i = 2; i < _node_index+2; i++) accumulator += i ;
                    beta_index = _node_index + 1 + _node_index * _n_nodes - accumulator + (j - _node_index - 1);
                    // Rcpp::Rcout << "j: " << j << ", beta_index:" << beta_index << ", beta:"<< theta(beta_index) << ", data:"<< data(j) <<"\n";

                }


                // shift by nodes params
                beta_index = _n_nodes + beta_index - 1;
                if(VERBOSEFLAG)Rcpp::Rcout << "j: " << j <<", beta_index:" << beta_index << ", beta:"<< _theta(beta_index) << ", data:"<< _data(j) <<"\n";

                eta += _theta(beta_index) * _data(j);
        }




    }

    return eta;
}

double binarynodeClass::ll_(const bool VERBOSEFLAG){

    double ll = log(exp(_data(_node_index) * _eta)/(1 + exp(_eta)));
    return ll;
}

Eigen::VectorXd binarynodeClass::gradient_(const bool VERBOSEFLAG){

    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(_d);

    gradient[_node_index] = _data[_node_index] - exp(_eta)/(1+exp(_eta));
    for(unsigned int j = 0; j < _n_nodes; j++){
        if(_node_index!=j){
            unsigned int beta_index;
            if(j < _node_index){
                if(j == 0){
                    beta_index = _node_index;
                    // Rcpp::Rcout << "j: " << j <<", beta_index:" << beta_index << ", beta:"<< theta(beta_index) << ", data:"<< data(j) <<"\n";

                }else{
                    unsigned int accumulator = 0;
                    for(unsigned int i = 2; i < j+2; i++) accumulator += i;
                    beta_index = _node_index + j * _n_nodes - accumulator;
                    // Rcpp::Rcout << "j: " << j << ", beta_index:" << beta_index << ", beta:"<< theta(beta_index) << ", data:"<< data(j) <<"\n";

                }
            }else
                if ( j > _node_index){
                    unsigned int accumulator = 0;
                    for(unsigned int i = 2; i < _node_index+2; i++) accumulator += i ;
                    beta_index = _node_index + 1 + _node_index * _n_nodes - accumulator + (j - _node_index - 1);
                    // Rcpp::Rcout << "j: " << j << ", beta_index:" << beta_index << ", beta:"<< theta(beta_index) << ", data:"<< data(j) <<"\n";

                }


                // shift by nodes params
                beta_index = _n_nodes + beta_index - 1;
                if(VERBOSEFLAG)Rcpp::Rcout << "j: " << j <<", beta_index:" << beta_index << ", beta:"<< _theta(beta_index) << ", data:"<< _data(j) <<"\n";

                if(_constraints[beta_index])gradient(beta_index) += _data(_node_index) * _data(j) - exp(_eta)/(1+exp(_eta))*_data(j);
        }



    }


    return gradient;
}
#endif

