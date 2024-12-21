#include <iostream>
#include <string>
#include <vector>
#include <armadillo>

using namespace arma;

class ShapeFunction{
    private:
        mat points;

    public:
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

class ElementMatrix{
    
    private: 
        double sound_speed = 1;

        mat points;
        ShapeFunction N;
        mat k_matrix;
        mat n_matrix;
    
    public:

        //k[ij]の計算をする
        double gauss_integral_m_two_spatial_derivative_term(int i,int j){

            double weight = 1;
            // xi plus minus
            double xi_pm = 0.577305;
            double eta_pm = 0.577305;
            double zeta_pm = 0.577305;

            //TODO もしかしてWeight3乗しなきゃ行けない？
            return 
                weight*integrand_spatial_derivative_term(xi_pm,eta_pm,zeta_pm,i,j)+
                weight*integrand_spatial_derivative_term(xi_pm,eta_pm,-zeta_pm,i,j)+ 
                weight*integrand_spatial_derivative_term(xi_pm,-eta_pm,zeta_pm,i,j)+
                weight*integrand_spatial_derivative_term(xi_pm,-eta_pm,-zeta_pm,i,j)+ 
                weight*integrand_spatial_derivative_term(-xi_pm,eta_pm,zeta_pm,i,j)+
                weight*integrand_spatial_derivative_term(-xi_pm,eta_pm,-zeta_pm,i,j)+ 
                weight*integrand_spatial_derivative_term(-xi_pm,-eta_pm,zeta_pm,i,j)+
                weight*integrand_spatial_derivative_term(-xi_pm,-eta_pm,-zeta_pm,i,j)
            ;
        }

        double gauss_integral_m_two_temporal_derivative_term(int i,int j){
            
            double weight = 1;
            // xi plus minus
            double xi_pm = 0.577305;
            double eta_pm = 0.577305;
            double zeta_pm = 0.577305;

            //TODO もしかしてWeight3乗しなきゃ行けない？
            return 
                weight*integrand_temporal_derivative_term(xi_pm,eta_pm,zeta_pm,i,j)+
                weight*integrand_temporal_derivative_term(xi_pm,eta_pm,-zeta_pm,i,j)+ 
                weight*integrand_temporal_derivative_term(xi_pm,-eta_pm,zeta_pm,i,j)+
                weight*integrand_temporal_derivative_term(xi_pm,-eta_pm,-zeta_pm,i,j)+ 
                weight*integrand_temporal_derivative_term(-xi_pm,eta_pm,zeta_pm,i,j)+
                weight*integrand_temporal_derivative_term(-xi_pm,eta_pm,-zeta_pm,i,j)+ 
                weight*integrand_temporal_derivative_term(-xi_pm,-eta_pm,zeta_pm,i,j)+
                weight*integrand_temporal_derivative_term(-xi_pm,-eta_pm,-zeta_pm,i,j)
            ;
        }

        double integrand_spatial_derivative_term(double xi,double eta,double zeta,int i,int j){
                //逆行列計算の回数を減らすために、一度Nxyzの計算結果を保存する
                mat Nxyz=N.xyz(xi,eta,zeta);
                vec Nx = Nxyz.row(0).t();
                vec Ny = Nxyz.row(1).t();
                vec Nz = Nxyz.row(2).t();
                double jacobian_det = N.jacobian_det(xi,eta,zeta);

            return sound_speed*sound_speed*(
                Nx(i)*Nx(j)+Ny(i)*Ny(j)+Nz(i)*Nz(j)
            )*jacobian_det;
        }

        double integrand_temporal_derivative_term(double xi,double eta,double zeta,int i,int j){
                //逆行列計算の回数を減らすために、一度Nxyzの計算結果を保存する
                vec N_func=N.N_func(xi,eta,zeta);
                double jacobian_det = N.jacobian_det(xi,eta,zeta);

            return N_func(i)*N_func(j)*jacobian_det;
        }

        void make_element_matrix(){

            mat k(8,8);
            mat n(8,8);

            for(int i=0;i<8;i++){
                for(int j=0;j<8;j++){
                    k(i,j) = gauss_integral_m_two_spatial_derivative_term(i,j);
                }
            }
            for(int i=0;i<8;i++){
                for(int j=0;j<8;j++){
                    n(i,j) = gauss_integral_m_two_temporal_derivative_term(i,j);
                }
            }

            k_matrix = arma::inv(n)*k;

        }

        ElementMatrix(mat init_points) : 
            N(init_points),
            k_matrix(8,8,arma::fill::zeros) 
        {
            points = init_points;

            make_element_matrix();

            std::cout<<k_matrix;

        }
};


class WaveMatrix{

    double sound_speed;

    vec point_vector;
    mat k_matrix;

    void cell_connectivity(){

    }

    void voundary_condittion(){

    }

    void solve_inv(){

    }

};

class Solver{

};

class Point{
    vec vector;
    int id;
};

class Cell{
    std::vector<int> point_id;
    int id;
};

class MeshUtils{

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

    ElementMatrix cell(points);

    std::cout << "シミュレーションが終了しました。\n";

    return 0;
}