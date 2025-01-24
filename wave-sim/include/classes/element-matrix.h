#pragma once

#include <iostream>
#include <vector>
#include <armadillo>
#include <shape-function.h>


class ElementMatrix{
    
    private: 
        double sound_speed = 1;
        arma::mat points;
        ShapeFunction N;
    
    public:
        arma::mat wave_matrix;
        arma::mat nodal_matrix;

        //k[ij]の計算をする
        double gauss_integral_m_two_spatial_derivative_term(int i,int j);
        double gauss_integral_m_two_temporal_derivative_term(int i,int j);
        double integrand_spatial_derivative_term(double xi,double eta,double zeta,int i,int j);
        double integrand_temporal_derivative_term(double xi,double eta,double zeta,int i,int j);
        void make_element_matrix();
        ElementMatrix(arma::mat init_points);
        ElementMatrix(arma::mat init_points,double init_sound_speed);
};