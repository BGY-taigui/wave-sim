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
    
    public:
        mat wave_matrix;
        mat nodal_matrix;

        //k[ij]の計算をする
        double gauss_integral_m_two_spatial_derivative_term(int i,int j){

            double weight = 1;
            // xi plus minus
            double xi_pm = 0.577305;
            double eta_pm = 0.577305;
            double zeta_pm = 0.577305;

            //TODO本当に2次のガウス積分であってるか
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

            for(int i=0;i<8;i++){
                for(int j=0;j<8;j++){
                    wave_matrix(i,j) = gauss_integral_m_two_spatial_derivative_term(i,j);
                }
            }
            for(int i=0;i<8;i++){
                for(int j=0;j<8;j++){
                    nodal_matrix(i,j) = gauss_integral_m_two_temporal_derivative_term(i,j);
                }
            }

            /*
            std::cout<<"print matrix"<<std::endl;
            std::cout<<wave_matrix<<std::endl;
            std::cout<<nodal_matrix<<std::endl;

            std::cout<<"total matrix"<<std::endl;
            std::cout<<arma::inv(nodal_matrix)*wave_matrix<<std::endl;
            */
        }

        ElementMatrix(mat init_points) : 
            N(init_points),
            wave_matrix(8,8,arma::fill::zeros),
            nodal_matrix(8,8,arma::fill::zeros) 
        {
            points = init_points;
            make_element_matrix();

        }
};


class Point{
    public:
        vec vector;
        int id;

        Point(vec init_vector,int init_id){
            vector = init_vector;
            id = init_id;
        }
};

class Cell{
    public:
        std::vector<int> point_ids;
        int id;

        Cell(std::vector<int> init_point_ids,int init_id){
            point_ids = init_point_ids;
            id = init_id;
        }
};

class GatherMatrix{
    public:
        sp_mat L;

        GatherMatrix(std::vector<int> element_point_ids,std::vector<int>& corresponding_point_ids):
            L(8,corresponding_point_ids.size()) 
        {
            for(int i=0;i<corresponding_point_ids.size();i++){
                for(int j=0;j<element_point_ids.size();j++){
                    if(element_point_ids[j]==corresponding_point_ids[i]){
                        L(j,i) = 1;
        }}}} 
};

class GlobalMatrix{

    public:

    double sound_speed;

    vec point_vector;
    std::vector<int> corresponding_point_ids;
    mat global_wave_matrix;
    sp_mat global_nodal_matrix;
    mat global_matrix;

    void boundary_condittion_zero_point(std::vector<Point> zero_value_points){

        for(auto& zero_value_point : zero_value_points){

        }

    }

    GlobalMatrix(std::vector<Cell>& mesh_cells,std::vector<Point>& mesh_points)
    {
        global_wave_matrix = sp_mat(mesh_points.size(),mesh_points.size());
        global_nodal_matrix = sp_mat(mesh_points.size(),mesh_points.size());

        corresponding_point_ids = std::vector<int>(mesh_points.size());
        for(int i=0;i<mesh_points.size();i++) {
            corresponding_point_ids[i]=mesh_points[i].id;
        }

        for(auto& cell : mesh_cells){
            //全体行列を一旦8行、列分0で拡大する

            mat points(8,3);
            //TODO ここら辺に、要素の数がズレていた時のエラー処理を実装する
            for(int i=0;i<cell.point_ids.size(); i++){
                //各要素のポイントの座標をまとめる
                points.row(i) = mesh_points[cell.point_ids[i]].vector.t();
            }

            ElementMatrix element_matrix(points);
            GatherMatrix gather_matrix(cell.point_ids,corresponding_point_ids);

            //要素行列を全体行列にマッピングする
            global_wave_matrix += gather_matrix.L.t() * element_matrix.wave_matrix * gather_matrix.L;
            global_nodal_matrix += gather_matrix.L.t() * element_matrix.nodal_matrix * gather_matrix.L;
        }

        std::cout<<"global single matrix"<<std::endl;
        //疎行列を使ってメモリ使用を効率化するなら、lapackじゃなくてsuperLUを使うように変更
        global_matrix = arma::spsolve(global_nodal_matrix,global_wave_matrix,"lapack");
        std::cout<<"end gathering"<<std::endl;

    }

};

class Solver{
     
    public:
    sp_mat global_matrix;

    void mode_analysis(){}
    void unsteadry_analysis(mat global_matrix,vec init_condition,double delta_t){
        int iter =3;
        std::vector<vec> point_value = {init_condition};

        for(int i=0;i<iter;i++){
            
            // u(i+1)=d2p/dt2 * DT**2 + 2u(i) -u(i-1)
            if(point_value.size()>1){
                std::cout<<"step:"<<i<<std::endl;
                point_value.push_back(
                    delta_t*delta_t*global_matrix*point_value[i] + 2* point_value[i] - point_value[i-1]
                );
            }else{

                std::cout<<"first step"<<std::endl;
                point_value.push_back(
                    delta_t*delta_t*global_matrix*point_value[i] + point_value[i]
                );
            }
            //std::cout<<(global_matrix*point_value[i]).t()<<std::endl;
            //std::cout<<point_value[i].t()<<std::endl;
            std::cout<<(global_matrix*point_value[i]).t()<<std::endl;
            std::cout<<point_value[i](0)<<","<<point_value[i](11)<<std::endl;
        } 

    }

};

class Rendere{};

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

    //ElementMatrix cell(points);
    /*
    std::vector<Point> mesh_points={
       Point(vec({0,0,0}),0),
       Point(vec({1,0,0}),1),
       Point(vec({1,1,0}),2),
       Point(vec({0,1,0}),3),
       Point(vec({0,0,1}),4),
       Point(vec({1,0,1}),5),
       Point(vec({1,1,1}),6),
       Point(vec({0,1,1}),7),

       Point(vec({0,0,2}),8),
       Point(vec({1,0,2}),9),
       Point(vec({1,1,2}),10),
       Point(vec({0,1,2}),11)
    }; 

    std::vector<Cell> mesh_cells={
        Cell(std::vector<int>{0,1,2,3,4,5,6,7},1),
        Cell(std::vector<int>{4,5,6,7,8,9,10,11},2)
    };

    GlobalMatrix global_matrix(mesh_cells,mesh_points);

    std::cout<<global_matrix.global_matrix<<std::endl;

    Solver solver;
    solver.unsteadry_analysis(global_matrix.global_matrix,vec({1,0,0,0,0,0,0,0,0,0,0,0}),0.1);

    std::cout << "シミュレーションが終了しました。\n";
    */
    std::vector<vec> point_value={vec({1,0,0,0,0,0,0,0})};
    for(int i=0;i<1;i++){
        point_value.push_back(
            arma::inv(ElementMatrix(points).nodal_matrix)*ElementMatrix(points).wave_matrix*point_value[i]
        );
    }

    std::cout<<point_value[1]<<std::endl;

    return 0;
}