#include <string>
#include <vector>

#include <armadillo>

#include <mesh-class.h>
#include <global-matrix.h>
#include <mesh-utils.h>
#include <solver.h>

using namespace arma;


int main(){

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
       Point(vec({0,1,2}),11),

       
       Point(vec({0,0,3}),12),
       Point(vec({1,0,3}),13),
       Point(vec({1,1,3}),14),
       Point(vec({0,1,3}),15),
       
       Point(vec({0,0,4}),16),
       Point(vec({1,0,4}),17),
       Point(vec({1,1,4}),18),
       Point(vec({0,1,4}),19)
    }; 

    std::vector<Cell> mesh_cells={
        Cell(std::vector<int>{0,1,2,3,4,5,6,7},1),
        Cell(std::vector<int>{4,5,6,7,8,9,10,11},2),
        Cell(std::vector<int>{8,9,10,11,12,13,14,15},3),
        Cell(std::vector<int>{12,13,14,15,16,17,18,19},4),
    };

    GlobalMatrix global_matrix(mesh_cells,mesh_points);

    std::cout<<global_matrix.global_matrix<<std::endl;

    Solver solver;
    solver.unsteadry_analysis(global_matrix.global_matrix,vec({1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}),0.01,1000);

    solver.mode_analysis(global_matrix.global_matrix);

    std::cout << "シミュレーションが終了しました。\n";

    MeshUtils mesh_utils;
    
    mesh_utils.write_mesh(mesh_points,mesh_cells,solver.point_values,solver.times);

    return 0;
}