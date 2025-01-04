#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <complex>
#include <numbers>

#include <armadillo>

#include <solver.h>

using namespace arma;

Solver::ModeAnalysisResult Solver::mode_analysis(mat global_matrix, int display_mode_num){
    
    Solver::ModeAnalysisResult result;

    result.display_mode_num = display_mode_num;

    arma::cx_vec eign_val;
    arma::cx_mat eign_vec;

    std::cout<<"computing eign value and vector"<<std::endl;

    arma::eig_gen(eign_val,eign_vec,global_matrix);

    arma::cx_vec angular_velocity(eign_val.size());
    arma::vec freqency(eign_val.size());

    for(int i=0;i<eign_val.size();i++){
        angular_velocity[i] = pow(eign_val[i],0.5);
        freqency[i] = pow(eign_val[i],0.5).imag()/(2 * M_PI);
    }

    std::cout<<"sorting modes"<<std::endl;
    arma::uvec sorted_index = arma::sort_index(freqency);

    result.freqency = freqency;

    //TODO これを引数にする もしくは自動にする
    int display_step_number = 100;

    double time_span = 10;
    double time_step = time_span/display_step_number;

    std::vector<double> times = {};

    for(int i=0;i<int(time_span/time_step);i++){
        times.push_back(time_step * i);
    }


    result.times = times;
    result.values = std::vector<std::vector<arma::vec>>(display_mode_num,std::vector<vec>(times.size(),vec(eign_val.size())));

    //result.modes = std::vector<EachModeResult>(display_mode_num);

    std::vector<vec> values={};

    // 振動数が少ない順にソートされたモードを指定された数だけ出力する
    for(int i =0; i<display_mode_num;i++){
        for(int j=0; j<times.size();j++){

            cx_vec result_complex = eign_vec.col(sorted_index[i]) * exp(angular_velocity[sorted_index[i]]*times[j]);

            vec result_real(result_complex.size());
            for (int k=0;k<result_real.size();k++){
                result_real(k) = result_complex(k).real();
            }

            result.values[i][j] = result_real;
        }
    }

    return result;
}



Solver::UnsteadryAnalysisResult Solver::unsteadry_analysis(mat global_matrix,vec init_condition,double delta_t, int iter){

    UnsteadryAnalysisResult result;

    result.point_values = {init_condition};
    result.times ={};

    for(int i=0;i<iter;i++){
        
        // u(i+1)=d2p/dt2 * DT**2 + 2u(i) -u(i-1)
        if(result.point_values.size()>1){
            result.point_values.push_back(
                delta_t*delta_t*global_matrix*result.point_values[i] + 2* result.point_values[i] - result.point_values[i-1]
            );
        }else{
            result.point_values.push_back(
                delta_t*delta_t*global_matrix*result.point_values[i] + result.point_values[i]
            );
        }

        result.times.push_back(delta_t*(i+1));

    } 


    return result;
}