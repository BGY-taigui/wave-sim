
#include <vector>
#include <armadillo>
#include <algorithm>

#include <mesh-class.h>
#include <element-matrix.h>
#include <gather-matrix.h>

#include <global-matrix.h>

using namespace arma;
    
void GlobalMatrix::define_condittion_zero_point(arma::vec center_vec, double radius,std::vector<Point>& mesh_points){

    std::vector<int> zero_value_point_ids;

    for(int i=0;i<corresponding_point_ids.size();i++){
        if (arma::norm(center_vec - mesh_points[corresponding_point_ids[i]].vector) < radius){
            zero_value_point_ids.push_back(corresponding_point_ids[i]);
        }
    }

    boundary_condittion_zero_point(zero_value_point_ids);
}

void GlobalMatrix::boundary_condittion_zero_point(std::vector<int> zero_value_point_ids){

    for(int zero_value_point_id : zero_value_point_ids){
        int zero_value_point_index = corresponding_point_id_indexes[zero_value_point_id];
        for(int i=0;i<global_wave_matrix.n_rows;i++){
            global_matrix(zero_value_point_index,i) = 0;
            global_matrix(i,zero_value_point_index) = 0;
        }
    }
}

mat GlobalMatrix::get_single_global_matrix(bool use_superlu){

    std::cout<<"cerating global single matrix ..."<<std::endl;
    //疎行列を使ってメモリ使用を効率化するなら、lapackじゃなくてsuperLUを使うように変更

    global_matrix = arma::spsolve(global_nodal_matrix,global_wave_matrix,"lapack");
    std::cout<<"end gathering"<<std::endl;

    /*
    if (use_superlu){
        global_matrix = arma::spsolve(global_nodal_matrix,global_wave_matrix,"superlu");
        std::cout<<"end gathering"<<std::endl;
    }else{
        global_matrix = arma::spsolve(global_nodal_matrix,global_wave_matrix,"lapack");
        std::cout<<"end gathering"<<std::endl;
    }
    */

    return global_matrix;
}

void GlobalMatrix::get_corresponding_point_id_indexes(){

    std::vector<int> indexes(corresponding_point_ids.size()+1);

    for(int row=0;row<corresponding_point_ids.size();row++){
        indexes[corresponding_point_ids[row]] = row;
    }

    corresponding_point_id_indexes = indexes;

};


GlobalMatrix::GlobalMatrix(std::vector<Cell>& mesh_cells,std::vector<Point>& mesh_points)
{
    global_wave_matrix = sp_mat(mesh_points.size(),mesh_points.size());
    global_nodal_matrix = sp_mat(mesh_points.size(),mesh_points.size());

    corresponding_point_ids = std::vector<int>(mesh_points.size());
    for(int i=0;i<mesh_points.size();i++) {
        corresponding_point_ids[i]=mesh_points[i].id;
    }

    get_corresponding_point_id_indexes();
    
    for(auto& cell : mesh_cells){
        //全体行列を一旦8行、列分0で拡大する

        std::cout<<"mapping cell ID:"<<cell.id<<std::endl;

        mat points(8,3);
        //TODO ここら辺に、要素の数がズレていた時のエラー処理を実装する
        for(int i=0;i<cell.point_ids.size(); i++){
            //各要素のポイントの座標をまとめる
            points.row(i) = mesh_points[cell.point_ids[i]].vector.t();
        }

        ElementMatrix element_matrix(points,cell.sound_speed);

        //要素行列を全体行列にマッピングする

        /*
        GatherMatrix gather_matrix(cell.point_ids,corresponding_point_ids);
        global_wave_matrix += gather_matrix.L.t() * element_matrix.wave_matrix * gather_matrix.L;
        global_nodal_matrix += gather_matrix.L.t() * element_matrix.nodal_matrix * gather_matrix.L;
        */


        for(int i=0;i<8;i++){
            for(int j=0;j<8;j++){
                int global_row = corresponding_point_id_indexes[cell.point_ids[i]];
                int global_col = corresponding_point_id_indexes[cell.point_ids[j]];

                global_wave_matrix(global_row,global_col) = 
                    global_wave_matrix(global_row,global_col) + 
                    element_matrix.wave_matrix(i,j);

                global_nodal_matrix(global_row,global_col) = 
                    global_nodal_matrix(global_row,global_col) + 
                    element_matrix.nodal_matrix(i,j);

            }
        }
    }
}


// 一応何もしないデフォルトコンストラクターを用意した
GlobalMatrix::GlobalMatrix(){}