
#include <iostream>
#include <vector>
#include <armadillo>
#include <shape-function.h>

#include <element-matrix.h>

using namespace arma;


//k[ij]の計算をする
double ElementMatrix::gauss_integral_m_two_spatial_derivative_term(int i,int j){

    double weight = 1;
    // xi plus minus
    double xi_pm = 0.577305;
    double eta_pm = 0.577305;
    double zeta_pm = 0.577305;

    //TODO本当に2次のガウス積分であってるか
    return 
        std::pow(weight,3)*integrand_spatial_derivative_term(xi_pm,eta_pm,zeta_pm,i,j)+
        std::pow(weight,3)*integrand_spatial_derivative_term(xi_pm,eta_pm,-zeta_pm,i,j)+ 
        std::pow(weight,3)*integrand_spatial_derivative_term(xi_pm,-eta_pm,zeta_pm,i,j)+
        std::pow(weight,3)*integrand_spatial_derivative_term(xi_pm,-eta_pm,-zeta_pm,i,j)+ 
        std::pow(weight,3)*integrand_spatial_derivative_term(-xi_pm,eta_pm,zeta_pm,i,j)+
        std::pow(weight,3)*integrand_spatial_derivative_term(-xi_pm,eta_pm,-zeta_pm,i,j)+ 
        std::pow(weight,3)*integrand_spatial_derivative_term(-xi_pm,-eta_pm,zeta_pm,i,j)+
        std::pow(weight,3)*integrand_spatial_derivative_term(-xi_pm,-eta_pm,-zeta_pm,i,j)
    ;
}


double ElementMatrix::gauss_integral_m_two_temporal_derivative_term(int i,int j){
    
    double weight = 1;
    // xi plus minus
    double xi_pm = 0.577305;
    double eta_pm = 0.577305;
    double zeta_pm = 0.577305;

    return 
        std::pow(weight,3)*integrand_temporal_derivative_term(xi_pm,eta_pm,zeta_pm,i,j)+
        std::pow(weight,3)*integrand_temporal_derivative_term(xi_pm,eta_pm,-zeta_pm,i,j)+ 
        std::pow(weight,3)*integrand_temporal_derivative_term(xi_pm,-eta_pm,zeta_pm,i,j)+
        std::pow(weight,3)*integrand_temporal_derivative_term(xi_pm,-eta_pm,-zeta_pm,i,j)+ 
        std::pow(weight,3)*integrand_temporal_derivative_term(-xi_pm,eta_pm,zeta_pm,i,j)+
        std::pow(weight,3)*integrand_temporal_derivative_term(-xi_pm,eta_pm,-zeta_pm,i,j)+ 
        std::pow(weight,3)*integrand_temporal_derivative_term(-xi_pm,-eta_pm,zeta_pm,i,j)+
        std::pow(weight,3)*integrand_temporal_derivative_term(-xi_pm,-eta_pm,-zeta_pm,i,j)
    ;
}

double ElementMatrix::integrand_spatial_derivative_term(double xi,double eta,double zeta,int i,int j){
        //逆行列計算の回数を減らすために、一度Nxyzの計算結果を保存する
        mat Nxyz=N.xyz(xi,eta,zeta);
        vec Nx = Nxyz.row(0).t();
        vec Ny = Nxyz.row(1).t();
        vec Nz = Nxyz.row(2).t();
        double jacobian_det = N.jacobian_det(xi,eta,zeta);

    return -sound_speed*sound_speed*(
        Nx(i)*Nx(j)+Ny(i)*Ny(j)+Nz(i)*Nz(j)
    )*jacobian_det;
}

double ElementMatrix::integrand_temporal_derivative_term(double xi,double eta,double zeta,int i,int j){
        //逆行列計算の回数を減らすために、一度Nxyzの計算結果を保存する
        vec N_func=N.N_func(xi,eta,zeta);
        double jacobian_det = N.jacobian_det(xi,eta,zeta);

    return N_func(i)*N_func(j)*jacobian_det;
}


void ElementMatrix::make_element_matrix(){

    for(int i=0;i<8;i++){
        for(int j=0;j<8;j++){
            wave_matrix(i,j) = gauss_integral_m_two_spatial_derivative_term(i,j);
            nodal_matrix(i,j) = gauss_integral_m_two_temporal_derivative_term(i,j);
        }
    }

}

ElementMatrix::ElementMatrix(mat init_points) : 
    N(init_points),
    wave_matrix(8,8,arma::fill::zeros),
    nodal_matrix(8,8,arma::fill::zeros) 
{
    points = init_points;
    make_element_matrix();

}
