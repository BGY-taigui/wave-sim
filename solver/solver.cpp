#include <iostream>
#include <string>
#include <vector>
#include <armadillo>

using namespace arma;

class ShapeFunction{
    public:
        //自然座標系での微分値
        vec x;
        vec y;
        vec z;

        mat jacobian;
        mat jacobian_inv;

        double jacobian_det;

        //直行座標系での微分値
        vec xi = {
            -1/8,1/8,1/8,-1/8,-1/8,1/8,1/8,-1/8
        };
        vec eta = {
            -1/8,-1/8,1/8,1/8,-1/8,-1/8,1/8,1/8
        };
        vec zeta = {
            -1/8,-1/8,-1/8,-1/8,1/8,1/8,1/8,1/8,
        };


    ShapeFunction(mat points){

        jacobian = mat(3,3,arma::fill::zeros);
        make_jacobian(points);
        calcurate_jacobian_inv();
        
        calcurate_Nxyz();

    }

    void make_jacobian(mat points){

        jacobian.row(0) = arma::dot(xi,points);
        jacobian.row(1) = arma::dot(eta,points);
        jacobian.row(2) = arma::dot(zeta,points);

    }

    void calcurate_jacobian_inv(){

        jacobian_inv = arma::inv(jacobian);

    }

    void calcurate_Nxyz(){

        mat N_xietazeta(3,3);
        mat Nxyz(3,3);

        N_xietazeta.row(0)=xi;
        N_xietazeta.row(1)=eta;
        N_xietazeta.row(2)=zeta;

        Nxyz = arma::dot(jacobian_inv,N_xietazeta);

        x = Nxyz.row(0);
        y = Nxyz.row(1);
        z = Nxyz.row(2);

    }

    void calcurate_jacobian_det(){
        jacobian_det = arma::det(jacobian);
    }
};

class Cell{
    public:

        mat points;
        mat element_matrix;

        double gauss_integral(){

    
            return 0.0;
        }

        double calcurate_weak_form_integration_term(){

            return 0.0;
        }

        void make_element_matrix(){

        }

        Cell(mat init_points){
            points = init_points;

            ShapeFunction N(points);

            calcurate_weak_form_integration_term();

            make_element_matrix();

        }
};


int main(){

    return 0;
}