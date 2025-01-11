#include <string>
#include <vector>

#include <armadillo>

#include <mesh-class.h>
#include <global-matrix.h>
#include <mesh-utils.h>
#include <solver.h>

using namespace arma;


int main(){


    //TODO メッシュ読み込みを実装する
    //TODO condition_importerを実装する

    MeshUtils meshutils;
    //MeshUtils::points_cells point_cells = meshutils.read_mesh("aa.vtk");

    //meshutils.read_vtm_file("aa.vtm");
    MeshUtils::points_cells point_cells = meshutils.read_cgns_file("aa.cgns");


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
       Point(vec({0,1,2}),11),

       
       Point(vec({0,0,2.7}),12),
       Point(vec({1,0,2.7}),13),
       Point(vec({1,1,3.3}),14),
       Point(vec({0,1,3.3}),15),
       
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
    MeshUtils::points_cells point_cells;
    point_cells.cells = mesh_cells;
    point_cells.points = mesh_points;

    */

    std::cout<<"creating Global Matrix ..."<<std::endl;
    GlobalMatrix global_matrix(point_cells.cells,point_cells.points);

    global_matrix.get_single_global_matrix(true);

    Solver solver;
    /*
    Solver::UnsteadryAnalysisResult unsteadry_analysis_result = solver.unsteadry_analysis(global_matrix.global_matrix,vec({1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}),0.01,1000);
    
    */

   solver.mode_analysis_frequency(sp_mat(global_matrix.global_wave_matrix),global_matrix.global_nodal_matrix);

    std::cout<<"Mode Analysing ..."<<std::endl;
    Solver::ModeAnalysisResult mode_analysis_result= solver.mode_analysis(global_matrix.global_matrix,10);

    MeshUtils mesh_utils;
    
    /*
    mesh_utils.write_mesh(mesh_points,mesh_cells,unsteadry_analysis_result.point_values,unsteadry_analysis_result.times,"output_datas");
    */

    for(int i=0;i<mode_analysis_result.display_mode_num;i++){
        mesh_utils.write_mesh(point_cells.points,point_cells.cells,mode_analysis_result.modes[i].point_values,mode_analysis_result.modes[i].times,"output_mode"+std::to_string(i+1));
    }


    return 0;

}