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
    std::vector<int> corresponding_point_ids;
    arma::mat global_wave_matrix;
    arma::sp_mat global_nodal_matrix;
    arma::mat global_matrix;

    void boundary_condittion_zero_point(std::vector<Point> zero_value_points);

    GlobalMatrix(std::vector<Cell>& mesh_cells,std::vector<Point>& mesh_points);

};