#include <iostream>
#include <string>
#include <vector>
#include <armadillo>

using namespace arma;

class ShapeFunction{
    public:
        mat points;

        ShapeFunction(mat points_init){
            points = points_init;
        }

        vec N_func(double xi,double eta,double zeta){
            vec result = {
                1/8*(1-xi)*(1-eta)*(1-zeta),
                1/8*(1+xi)*(1-eta)*(1-zeta),
                1/8*(1+xi)*(1+eta)*(1-zeta),
                1/8*(1-xi)*(1+eta)*(1-zeta),
                1/8*(1-xi)*(1-eta)*(1+zeta),
                1/8*(1+xi)*(1-eta)*(1+zeta),
                1/8*(1+xi)*(1+eta)*(1+zeta),
                1/8*(1-xi)*(1+eta)*(1+zeta),
            };
            return result;
        }

        vec diff_xi(double xi,double eta,double zeta){
            vec result = {
                1/8*-1*(1-eta)*(1-zeta),
                1/8*(1-eta)*(1-zeta),
                1/8*(1+eta)*(1-zeta),
                1/8*-1*(1+eta)*(1-zeta),
                1/8*-1*(1-eta)*(1+zeta),
                1/8*(1-eta)*(1+zeta),
                1/8*(1+eta)*(1+zeta),
                1/8*-1*(1+eta)*(1+zeta),
            };
            return result;
        }

        vec diff_eta(double xi,double eta,double zeta){
            vec result = {
                1/8*(1-xi)*-1*(1-zeta),
                1/8*(1+xi)*-1*(1-zeta),
                1/8*(1+xi)*(1-zeta),
                1/8*(1-xi)*(1-zeta),
                1/8*(1-xi)*-1*(1+zeta),
                1/8*(1+xi)*-1*(1+zeta),
                1/8*(1+xi)*(1+zeta),
                1/8*(1-xi)*(1+zeta),
            };
            return result;
        }

        vec diff_zeta(double xi,double eta,double zeta){
            vec result = {
                1/8*(1-xi)*(1-eta)*-1,
                1/8*(1+xi)*(1-eta)*-1,
                1/8*(1+xi)*(1+eta)*-1,
                1/8*(1-xi)*(1+eta)*-1,
                1/8*(1-xi)*(1-eta),
                1/8*(1+xi)*(1-eta),
                1/8*(1+xi)*(1+eta),
                1/8*(1-xi)*(1+eta),
            };
            return result;
        }

        mat diff_xietazeta(double xi,double eta,double zeta){
            mat result(8,3);
            result.row(0) = diff_xi(xi,eta,zeta);
            result.row(1) = diff_eta(xi,eta,zeta);
            result.row(2) = diff_zeta(xi,eta,zeta);

            return result;
        }

        mat jacobian(double xi,double eta, double zeta){

            mat jacobian(3,3);

            jacobian.row(0) = {
                arma::dot(diff_xi(xi,eta,zeta),points.col(0)),
                arma::dot(diff_xi(xi,eta,zeta),points.col(0)),
                arma::dot(diff_xi(xi,eta,zeta),points.col(0))
            };

            jacobian.row(1) = {
                arma::dot(diff_eta(xi,eta,zeta),points.col(1)),
                arma::dot(diff_eta(xi,eta,zeta),points.col(1)),
                arma::dot(diff_eta(xi,eta,zeta),points.col(1))
            };

            jacobian.row(2) = {
                arma::dot(diff_zeta(xi,eta,zeta),points.col(2)),
                arma::dot(diff_zeta(xi,eta,zeta),points.col(2)),
                arma::dot(diff_zeta(xi,eta,zeta),points.col(2))
            };
            //return jacobian;

        }

        mat jacobian_inv(double xi,double eta,double zeta){
            return arma::inv(jacobian(xi,eta,zeta));
        }

        mat xyz(double xi,double eta,double zeta){
            return jacobian_inv(xi,eta,zeta)*diff_xietazeta(xi,eta,zeta);
        }

        double calcurate_jacobian_det(double xi,double eta,double zeta){
            return arma::det(jacobian(xi,eta,zeta));
        }
};

class Cell{
    public:
    
        double sound_speed;

        mat points;
        mat element_matrix;

        ShapeFunction N;

        double gauss_integral(){

    
            return 0.0;
        }

        double weak_form_int_part(double xi,double eta,double zeta){

        }

        double calcurate_weak_form_integration_term(){

            return 0.0;
        }

        void make_element_matrix(){

        }

        Cell(mat init_points){
            points = init_points;

            N = ShapeFunction(points);

            calcurate_weak_form_integration_term();

            make_element_matrix();

        }
};


int main(){

    return 0;
}