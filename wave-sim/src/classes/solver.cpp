#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <complex>
#include <numbers>

#include <armadillo>

#include <solver.h>

using namespace arma;

void Solver::mode_analysis_frequency(sp_mat global_wave_matrix,sp_mat gloabl_nodal_matrix){

    double start_freqency = 0;
    double end_freqency = 0.1;
    double freqency_step = (end_freqency - start_freqency)/10;
    std::vector<double> det_values(0);

    std::cout<<"w m:"<<global_wave_matrix.max()<<"n m"<<gloabl_nodal_matrix.max()<<std::endl;
    std::cout<<det(mat(global_wave_matrix))<<std::endl;
    std::cout<<det(mat(gloabl_nodal_matrix))<<std::endl;
    std::cout<<mat(global_wave_matrix/global_wave_matrix.max()).max()<<std::endl;

    for(double freqency = start_freqency;freqency<end_freqency;freqency+=freqency_step){

        std::cout<<"computing det value. freq:"<<freqency;

        double alpha = 2 * M_PI * freqency;

        mat mat = arma::mat(global_wave_matrix+pow(alpha,2)*gloabl_nodal_matrix);

        arma::mat L, U, P;
        bool success = arma::lu(L, U, P, mat); // LU分解

        std::cout<<"LU success:"<<success;

        // 対数行列式の計算
        //double log_det_U = arma::sum(arma::log(U.diag())); // U の対角成分の対数の和

        std::complex<double> log_det_U(0.0, 0.0); // 対数行列式の初期値

        for (arma::uword i = 0; i < U.n_rows; ++i) {
            double diag_elem = U(i, i);

            if (diag_elem > 0) {
                log_det_U += std::log(diag_elem); // 正の対角成分はそのまま
            } else if (diag_elem < 0) {
                log_det_U += std::log(-diag_elem) + std::complex<double>(0.0, arma::datum::pi); // 負の対角成分
            } else {
                std::cerr << "Matrix is singular (zero diagonal element)." << std::endl;
            }
        }



        //float det_value = log10(arma::det(mat));
        //float det_value = arma::det(mat);

        std::cout<<" det:"<<log_det_U<<std::endl;

        //det_values.push_back(log_det_U.real());
    }

    std::cout<<"eig pair"<<std::endl;
    arma::eig_pair(arma::mat(global_wave_matrix),arma::mat(gloabl_nodal_matrix));

}



Solver::ModeAnalysisResult Solver::mode_analysis(mat global_matrix, int display_mode_num){
    
    Solver::ModeAnalysisResult result;

    result.display_mode_num = display_mode_num;

    arma::cx_vec eign_val;
    arma::cx_mat eign_vec;

    //std::cout<<"computing eign value and vector"<<std::endl;

    //arma::eig_gen(eign_val,eign_vec,global_matrix);

    std::cout<<"computing sp mat"<<std::endl;

    arma::sp_mat sp_gm(global_matrix);

    std::cout<<"computing eign limited value and vector"<<std::endl;

    arma::eigs_gen(eign_val,eign_vec,sp_gm,display_mode_num,"sm");


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
    int display_time_step_number = 100;

    result.modes = std::vector<EachModeResult>(display_mode_num);

    // 振動数が少ない順にソートされたモードを指定された数だけ出力する
    for(int i =0; i<display_mode_num;i++){

        std::vector<vec> each_mode_point_values(display_time_step_number);
        std::vector<double> each_mode_times(display_time_step_number);

        double time_span = 1/freqency(sorted_index[i]);
        double time_step = time_span/display_time_step_number;

        for(int j=0;j<display_time_step_number;j++){
            each_mode_times[j]= time_step * j;
        }

        for(int j=0; j<display_time_step_number;j++){

            cx_vec result_complex = eign_vec.col(sorted_index[i]) * exp(angular_velocity[sorted_index[i]]*each_mode_times[j]);

            vec result_real(result_complex.size());

            for (int k=0;k<result_real.size();k++){
                result_real(k) = result_complex(k).real();
            }

            each_mode_point_values[j] = result_real;
        }

        result.modes[i].point_values = each_mode_point_values;
        result.modes[i].times = each_mode_times;
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