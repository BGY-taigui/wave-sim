
#include <armadillo>
#include <vector>

#include <mesh-class.h>

using namespace arma;

Point::Point(vec init_vector,int init_id):
        vector(init_vector),
        id(init_id)
    {}

Cell::Cell(std::vector<int> init_point_ids,int init_id):
    point_ids(init_point_ids),
    id(init_id)
    {}