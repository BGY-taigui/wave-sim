#pragma once

#include <vector>
#include <armadillo>


class GatherMatrix{
    public:
        arma::sp_mat L;
        GatherMatrix(std::vector<int> element_point_ids,std::vector<int>& corresponding_point_ids);
};