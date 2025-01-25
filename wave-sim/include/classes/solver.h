#pragma once

#include <vector>
#include <string>
#include <iostream>

#include <armadillo>

class Solver{
     
    public:

    struct EachModeResult
    {
        std::vector<double> times;
        std::vector<arma::vec> point_values;
    };

    struct SortedEigs{
        arma::cx_vec eign_values;
        arma::cx_mat eign_vectors;
    };

    struct ModeAnalysisResult
    {
        arma::vec freqency;
        int display_mode_num;

        std::vector<EachModeResult> modes;
    };

    struct UnsteadryAnalysisResult
    {
        arma::sp_mat global_matrix;
        std::vector<arma::vec> point_values;
        std::vector<double> times;
    };

    void mode_analysis_frequency(arma::sp_mat global_wave_matrix,arma::sp_mat gloabl_nodal_matrix);

    SortedEigs compute_eigs(arma::mat global_matrix);

    ModeAnalysisResult mode_analysis(arma::mat global_matrix,int mode_num);
    UnsteadryAnalysisResult unsteadry_analysis(arma::mat global_matrix,arma::vec init_condition,double delta_t, int iter);


};
