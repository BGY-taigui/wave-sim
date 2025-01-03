#include <iostream>
#include <vector>
#include <armadillo>

#include <shape-function.h>

using namespace arma;

ShapeFunction::ShapeFunction(mat points_init):
    points(points_init)
    {}

arma::vec ShapeFunction::N_func(double xi,double eta,double zeta){
    vec result = {
        1.0/8*(1.0-xi)*(1.0-eta)*(1.0-zeta),
        1.0/8*(1.0+xi)*(1.0-eta)*(1.0-zeta),
        1.0/8*(1.0+xi)*(1.0+eta)*(1.0-zeta),
        1.0/8*(1.0-xi)*(1.0+eta)*(1.0-zeta),
        1.0/8*(1.0-xi)*(1.0-eta)*(1.0+zeta),
        1.0/8*(1.0+xi)*(1.0-eta)*(1.0+zeta),
        1.0/8*(1.0+xi)*(1.0+eta)*(1.0+zeta),
        1.0/8*(1.0-xi)*(1.0+eta)*(1.0+zeta),
    };
    return result;
}

vec ShapeFunction::diff_xi(double xi,double eta,double zeta){
    vec result = {
        1.0/8*-1.0*(1.0-eta)*(1.0-zeta),
        1.0/8*(1.0-eta)*(1.0-zeta),
        1.0/8*(1.0+eta)*(1.0-zeta),
        1.0/8*-1.0*(1.0+eta)*(1.0-zeta),
        1.0/8*-1.0*(1.0-eta)*(1.0+zeta),
        1.0/8*(1.0-eta)*(1.0+zeta),
        1.0/8*(1.0+eta)*(1.0+zeta),
        1.0/8*-1.0*(1.0+eta)*(1.0+zeta),
    };
    return result;
}

vec ShapeFunction::diff_eta(double xi,double eta,double zeta){
    vec result = {
        1.0/8*(1.0-xi)*-1.0*(1.0-zeta),
        1.0/8*(1.0+xi)*-1.0*(1.0-zeta),
        1.0/8*(1.0+xi)*(1.0-zeta),
        1.0/8*(1.0-xi)*(1.0-zeta),
        1.0/8*(1.0-xi)*-1.0*(1.0+zeta),
        1.0/8*(1.0+xi)*-1.0*(1.0+zeta),
        1.0/8*(1.0+xi)*(1.0+zeta),
        1.0/8*(1.0-xi)*(1.0+zeta),
    };
    return result;
}

vec ShapeFunction::diff_zeta(double xi,double eta,double zeta){
    vec result = {
        1.0/8*(1.0-xi)*(1.0-eta)*-1.0,
        1.0/8*(1.0+xi)*(1.0-eta)*-1.0,
        1.0/8*(1.0+xi)*(1.0+eta)*-1.0,
        1.0/8*(1.0-xi)*(1.0+eta)*-1.0,
        1.0/8*(1.0-xi)*(1.0-eta),
        1.0/8*(1.0+xi)*(1.0-eta),
        1.0/8*(1.0+xi)*(1.0+eta),
        1.0/8*(1.0-xi)*(1.0+eta),
    };
    return result;
}

        
mat ShapeFunction::diff_xietazeta(double xi,double eta,double zeta){
    mat result(3,8);
    result.row(0) = diff_xi(xi,eta,zeta).t();
    result.row(1) = diff_eta(xi,eta,zeta).t();
    result.row(2) = diff_zeta(xi,eta,zeta).t();

    return result;
}


//直交座標系で指定された位置のヤコビアンを計算する
//自然から直行座標への変換は線形変換ではないので、ヤコビアンは座標によって変わる
mat ShapeFunction::jacobian(double xi,double eta, double zeta){

    mat jacobian(3,3);

    jacobian.row(0) = {
        arma::dot(diff_xi(xi,eta,zeta),points.col(0)),
        arma::dot(diff_eta(xi,eta,zeta),points.col(0)),
        arma::dot(diff_zeta(xi,eta,zeta),points.col(0))
    };

    jacobian.row(1) = {
        arma::dot(diff_xi(xi,eta,zeta),points.col(1)),
        arma::dot(diff_eta(xi,eta,zeta),points.col(1)),
        arma::dot(diff_zeta(xi,eta,zeta),points.col(1))
    };

    jacobian.row(2) = {
        arma::dot(diff_xi(xi,eta,zeta),points.col(2)),
        arma::dot(diff_eta(xi,eta,zeta),points.col(2)),
        arma::dot(diff_zeta(xi,eta,zeta),points.col(2))
    };
    return jacobian;

}


mat ShapeFunction::jacobian_inv(double xi,double eta,double zeta){
    return arma::inv(jacobian(xi,eta,zeta));
}


//直交座標系での偏微分を自然座標系に変換して、結果を行列で返す
mat ShapeFunction::xyz(double xi,double eta,double zeta){
    return jacobian_inv(xi,eta,zeta)*diff_xietazeta(xi,eta,zeta);
}

double ShapeFunction::jacobian_det(double xi,double eta,double zeta){
    return arma::det(jacobian(xi,eta,zeta));
}