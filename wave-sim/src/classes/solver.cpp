#include <vector>
#include <string>
#include <iostream>

#include <armadillo>

#include <solver.h>

using namespace arma;

void Solver::mode_analysis(mat global_matrix){

    arma::cx_vec eign_val;
    arma::cx_mat eign_vec;

    arma::eig_gen(eign_val,eign_vec,global_matrix);

    std::cout<<eign_val<<std::endl;

}

void Solver::unsteadry_analysis(mat global_matrix,vec init_condition,double delta_t, int iter){
    point_values = {init_condition};
    times ={};

    for(int i=0;i<iter;i++){
        
        // u(i+1)=d2p/dt2 * DT**2 + 2u(i) -u(i-1)
        if(point_values.size()>1){
            point_values.push_back(
                -delta_t*delta_t*global_matrix*point_values[i] + 2* point_values[i] - point_values[i-1]
            );
        }else{
            point_values.push_back(
                -delta_t*delta_t*global_matrix*point_values[i] + point_values[i]
            );
        }

        times.push_back(delta_t*(i+1));

    } 
}