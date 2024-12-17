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
                1.0/8*(1-xi)*(1-eta)*(1-zeta),
                1.0/8*(1+xi)*(1-eta)*(1-zeta),
                1.0/8*(1+xi)*(1+eta)*(1-zeta),
                1.0/8*(1-xi)*(1+eta)*(1-zeta),
                1.0/8*(1-xi)*(1-eta)*(1+zeta),
                1.0/8*(1+xi)*(1-eta)*(1+zeta),
                1.0/8*(1+xi)*(1+eta)*(1+zeta),
                1.0/8*(1-xi)*(1+eta)*(1+zeta),
            };
            return result;
        }

        //形状関数をxiで微分した値(直交座標系の位置によって値が変わる)
        vec diff_xi(double xi,double eta,double zeta){
            vec result = {
                1.0/8*-1*(1-eta)*(1-zeta),
                1.0/8*(1-eta)*(1-zeta),
                1.0/8*(1+eta)*(1-zeta),
                1.0/8*-1*(1+eta)*(1-zeta),
                1.0/8*-1*(1-eta)*(1+zeta),
                1.0/8*(1-eta)*(1+zeta),
                1.0/8*(1+eta)*(1+zeta),
                1.0/8*-1*(1+eta)*(1+zeta),
            };
            return result;
        }

        vec diff_eta(double xi,double eta,double zeta){
            vec result = {
                1.0/8*(1-xi)*-1*(1-zeta),
                1.0/8*(1+xi)*-1*(1-zeta),
                1.0/8*(1+xi)*(1-zeta),
                1.0/8*(1-xi)*(1-zeta),
                1.0/8*(1-xi)*-1*(1+zeta),
                1.0/8*(1+xi)*-1*(1+zeta),
                1.0/8*(1+xi)*(1+zeta),
                1.0/8*(1-xi)*(1+zeta),
            };
            return result;
        }

        vec diff_zeta(double xi,double eta,double zeta){
            vec result = {
                1.0/8*(1-xi)*(1-eta)*-1,
                1.0/8*(1+xi)*(1-eta)*-1,
                1.0/8*(1+xi)*(1+eta)*-1,
                1.0/8*(1-xi)*(1+eta)*-1,
                1.0/8*(1-xi)*(1-eta),
                1.0/8*(1+xi)*(1-eta),
                1.0/8*(1+xi)*(1+eta),
                1.0/8*(1-xi)*(1+eta),
            };
            return result;
        }

        mat diff_xietazeta(double xi,double eta,double zeta){
            mat result(3,8);
            result.row(0) = diff_xi(xi,eta,zeta).t();
            result.row(1) = diff_eta(xi,eta,zeta).t();
            result.row(2) = diff_zeta(xi,eta,zeta).t();

            return result;
        }

        //直交座標系で指定された位置のヤコビアンを計算する
        //自然から直行座標への変換は線形変換ではないので、ヤコビアンは座標によって変わる
        mat jacobian(double xi,double eta, double zeta){

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

        mat jacobian_inv(double xi,double eta,double zeta){
            return arma::inv(jacobian(xi,eta,zeta));
        }

        //直交座標系での偏微分を自然座標系に変換して、結果を行列で返す
        mat xyz(double xi,double eta,double zeta){
            return jacobian_inv(xi,eta,zeta)*diff_xietazeta(xi,eta,zeta);
        }

        double jacobian_det(double xi,double eta,double zeta){
            return arma::det(jacobian(xi,eta,zeta));
        }
};

class Cell{
    public:
    
        double sound_speed;

        mat points;
        mat element_matrix;

        ShapeFunction N;

        //k[ij]の計算をする
        double gauss_integral_m_two(double xi,double eta,double zeta,int i,int j){
            double weight = 1;
            // xi plus minus
            double xi_pm = 0.577305;
            double eta_pm = 0.577305;
            double zeta_pm = 0.577305;

            return 
                weight*integrand(xi_pm,eta_pm,zeta_pm,i,j)+weight*integrand(xi_pm,eta_pm,-zeta_pm,i,j)+ 
                weight*integrand(xi_pm,-eta_pm,zeta_pm,i,j)+weight*integrand(xi_pm,-eta_pm,-zeta_pm,i,j)+ 
                weight*integrand(-xi_pm,eta_pm,zeta_pm,i,j)+weight*integrand(-xi_pm,eta_pm,-zeta_pm,i,j)+ 
                weight*integrand(-xi_pm,-eta_pm,zeta_pm,i,j)+weight*integrand(-xi_pm,-eta_pm,-zeta_pm,i,j)
            ;
        }

        double integrand(double xi,double eta,double zeta,int i,int j){
                //逆行列計算の回数を減らすために、一度Nxyzの計算結果を保存する
                mat Nxyz=N.xyz(xi,eta,zeta);
                vec Nx = Nxyz.row(0);
                vec Ny = Nxyz.row(1);
                vec Nz = Nxyz.row(2);
                double jacobian = N.jacobian_det(xi,eta,zeta);

            return sound_speed*sound_speed*(
                Nx(i)*Nx(j)+Ny(i)*Ny(j)+Nz(i)*Nz(j)
            )*jacobian;
        }

        void make_element_matrix(){

        }

        Cell(mat init_points) : N(init_points) {
            points = init_points;


            make_element_matrix();

        }
};


class WaveMatrix{

};

int main(){

    mat points = {
        {0,0,0},
        {1,0,0},
        {1,1,0},
        {0,1,0},
        {0,0,1},
        {1,0,1},
        {1,1,1},
        {0,1,1},
    };

    ShapeFunction N(points);

    std::cout << "シミュレーションが終了しました。\n";

    return 0;
}