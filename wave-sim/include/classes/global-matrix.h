#pragma once

#include <vector>
#include <armadillo>


#include <mesh-class.h>
#include <element-matrix.h>
#include <gather-matrix.h>


class GlobalMatrix{

    public:

    double sound_speed;

    arma::vec point_vector;
    // index:row value:point_id
    std::vector<int> corresponding_point_ids;
    // index:point_id value:row
    std::vector<int> corresponding_point_id_indexes;
    arma::mat global_wave_matrix;
    arma::sp_mat global_nodal_matrix;
    arma::mat global_matrix;

    void boundary_condittion_zero_point(std::vector<int> zero_value_point_ids);

    void get_corresponding_point_id_indexes();

    arma::mat get_single_global_matrix(bool use_superlu);

    GlobalMatrix(std::vector<Cell>& mesh_cells,std::vector<Point>& mesh_points);

};