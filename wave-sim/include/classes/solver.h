#pragma once

#include <vector>
#include <string>
#include <iostream>

#include <armadillo>

class Solver{
     
    public:
    arma::sp_mat global_matrix;
    std::vector<arma::vec> point_values;
    std::vector<double> times;

    void mode_analysis(arma::mat global_matrix);
    void unsteadry_analysis(arma::mat global_matrix,arma::vec init_condition,double delta_t, int iter);
};
