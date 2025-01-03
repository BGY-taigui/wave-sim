
#include <vector>
#include <armadillo>

#include <mesh-class.h>
#include <element-matrix.h>
#include <gather-matrix.h>

#include <global-matrix.h>

using namespace arma;
    

void GlobalMatrix::boundary_condittion_zero_point(std::vector<Point> zero_value_points){

    for(auto& zero_value_point : zero_value_points){

    }

}


GlobalMatrix::GlobalMatrix(std::vector<Cell>& mesh_cells,std::vector<Point>& mesh_points)
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
