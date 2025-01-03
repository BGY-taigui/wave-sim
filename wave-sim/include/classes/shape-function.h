#pragma once

#include <iostream>
#include <vector>
#include <armadillo>

class ShapeFunction{
    private:
        arma::mat points;

    public:
        ShapeFunction(arma::mat points_init);
        arma::vec N_func(double xi,double eta,double zeta);
        arma::vec diff_xi(double xi,double eta,double zeta);
        arma::vec diff_eta(double xi,double eta,double zeta);
        arma::vec diff_zeta(double xi,double eta,double zeta);
        arma::mat diff_xietazeta(double xi,double eta,double zeta);
        arma::mat jacobian(double xi,double eta, double zeta);
        arma::mat jacobian_inv(double xi,double eta,double zeta);
        arma::mat xyz(double xi,double eta,double zeta);
        double jacobian_det(double xi,double eta,double zeta);
};
