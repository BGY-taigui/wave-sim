#pragma once

#include <armadillo>
#include <vector>

class Point{
    public:
        arma::vec vector;
        int id;
        Point(arma::vec init_vector,int init_id);
};

class Cell{
    public:
        std::vector<int> point_ids;
        int id;

        Cell(std::vector<int> init_point_ids,int init_id);
};