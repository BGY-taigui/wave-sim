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

    int mode_num(const arma::cx_vec& eign_val,const arma::cx_mat& eign_vec,double start_mode_freq,double end_mode_freq);

    std::vector<EachModeResult> output_time_series_mode_num(const arma::cx_vec& eign_val,const arma::cx_mat& eign_vec,int start_mode_num,int end_mode_num,int display_timestep_num);
    std::vector<EachModeResult> output_time_series_freq_range(const arma::cx_vec& eign_val,const arma::cx_mat& eign_vec,double start_mode_freq,double end_mode_freq,int display_timestep_num);

    ModeAnalysisResult mode_analysis(arma::mat global_matrix,int mode_num);
    UnsteadryAnalysisResult unsteadry_analysis(arma::mat global_matrix,arma::vec init_condition,double delta_t, int iter);


};
