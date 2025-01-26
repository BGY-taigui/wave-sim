#include <string>
#include <vector>
#include <unistd.h>
#include <getopt.h>

#include <armadillo>

#include <mesh-class.h>
#include <global-matrix.h>
#include <mesh-utils.h>
#include <solver.h>
#include <project-operator.h>

using namespace arma;


int main(int argc , char* argv[]){


    /* TODO

        新しく実装が必要なクラス
        プロジェクトオペレーター
        レンダラー
        MeshUtilsのメッシュのデータ構造解析機能
        境界条件のファイル解析機能

        1. wave-eq newproject projectname
        2. wave-eq read-mesh -f cgnsfilename.cgns -p projectfile.pjtORprojectname
        3. wave-eq loadbc -f bcfilename.csv -p projectfile.pjtORprojectname
        4. wave-eq run-modeanalysis -p projectfile.pjtORprojectname
        5. wave-eq render -p projectfile.pjtORprojectname
        
         . wave-eq describe -f meshfilename.cgnsORvtkORvtm
         . wave-eq describe -p projectfile.pjtORprojectname

        コマンドオプション
            f filename
            m mode analysis
            u unsteadry analysis


        コマンドユーティリティとして実装したいこと
        　・プロジェクトファイルを最初に作成する
        　・プロジェクトファイルにファイルなどを書き込んでいく
        　・大規模計算を行うときは、サブコマンドとプロジェクトファイルの指定で行う
        　・大規模計算の結果は、専用のファイルに出力し、プロジェクトファイルに内容を記録しておく
        　・プロジェクトファイルの内容を出力して表示する機能を実装
        　・メッシュファイルの内容を出力して表示する機能を実装
    */

    std::string mesh_filename = "aa,cgns";

    bool do_mode_analysis = false;
    bool do_rendering = false;

    std::string filename;
    std::string projectnameorfilename;


    std::string subcommand =argv[1];

    int opt;

    if(subcommand == "newproject"){
        while ((opt = getopt(argc - 1, argv + 1, "p:")) != -1) {
            switch (opt) {
                case 'p':
                    projectnameorfilename = optarg;
                    break;
                default:
                    std::cerr << "Usage: wave-eq newproject -p <project>" << std::endl;
                    return 1;
            }
        }

        ProjectOperator project_operator;
        project_operator.NewProject(projectnameorfilename);
        

    }else if(subcommand == "read-mesh"){
        while ((opt = getopt(argc - 1, argv + 1, "f:p:")) != -1) {
            switch (opt) {
                case 'f':
                    filename = optarg;
                    break;
                case 'p':
                    projectnameorfilename = optarg;
                    break;
                default:
                    std::cerr << "Usage: wave-eq read-mesh -f <filename> -p <project>" << std::endl;
                    return 1;
            }
        }
        ProjectOperator project_operator(projectnameorfilename);
        project_operator.ReadMeshFile(filename);
    }else if(subcommand == "mode-analysis"){
        while ((opt = getopt(argc - 1, argv + 1, "f:p:")) != -1) {
            switch (opt) {
                case 'f':
                    filename = optarg;
                    break;
                case 'p':
                    projectnameorfilename = optarg;
                    break;
                default:
                    std::cerr << "Usage: wave-eq read-mesh -f <filename> -p <project>" << std::endl;
                    return 1;
            }
        }
        ProjectOperator project_operator(projectnameorfilename);
        project_operator.RunModeAnalysis();
    }else if(subcommand == "describe"){
        while ((opt = getopt(argc - 1, argv + 1, "f:")) != -1) {
            switch (opt) {
                case 'f':
                    filename = optarg;
                    break;
                default:
                    std::cerr << "Usage: wave-eq read-mesh -f <filename> -p <project>" << std::endl;
                    return 1;
            }
        }
        ProjectOperator project_operator;
        project_operator.Describe(filename);
    }else if(subcommand == "cleanall"){
        while ((opt = getopt(argc - 1, argv + 1, "p:")) != -1) {
            switch (opt) {
                case 'p':
                    projectnameorfilename = optarg;
                    break;
                default:
                    std::cerr << "Usage: wave-eq read-mesh -f <filename> -p <project>" << std::endl;
                    return 1;
            }
        }
        ProjectOperator project_operator(projectnameorfilename);
        project_operator.CleanAll();

    }else if(subcommand == "output-timeseries"){
        while ((opt = getopt(argc - 1, argv + 1, "p:")) != -1) {
            switch (opt) {
                case 'p':
                    projectnameorfilename = optarg;
                    break;
                default:
                    std::cerr << "Usage: wave-eq read-mesh -f <filename> -p <project>" << std::endl;
                    return 1;
            }
        }
        ProjectOperator project_operator(projectnameorfilename);
        project_operator.OutputTimeSeriesVTK();
    }else if(subcommand == "copy-project"){
        std::string origin_projectnamerofilename;
        while ((opt = getopt(argc - 1, argv + 1, "o:p:")) != -1) {
            switch (opt) {
                case 'o':
                    origin_projectnamerofilename = optarg;
                    break;
                case 'p':
                    projectnameorfilename = optarg;
                    break;
                default:
                    std::cerr << "Usage: wave-eq read-mesh -f <filename> -p <project>" << std::endl;
                    return 1;
            }
        }
        ProjectOperator project_operator(projectnameorfilename);
        project_operator.CopyProject(origin_projectnamerofilename);
    }


    /*

    //TODO condition_importerを実装する

    MeshUtils meshutils;
    //MeshUtils::points_cells point_cells = meshutils.read_mesh("aa.vtk");

    //meshutils.read_vtm_file("aa.vtm");
    MeshUtils::points_cells point_cells = meshutils.read_cgns_file(mesh_filename);

    for(int i=0;i<point_cells.cells.size();i++){
        point_cells.cells[i].sound_speed = 1466;
    }
    */

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

    /*
    std::cout<<"creating Global Matrix ..."<<std::endl;
    GlobalMatrix global_matrix(point_cells.cells,point_cells.points);

    global_matrix.get_single_global_matrix(true);

    //global_matrix.define_condittion_zero_point(vec({0,0,0}),0.3,point_cells.points);
    global_matrix.define_condittion_zero_point(
        vec({66.756378,-22.345131,-19.052559}),
        0.3,point_cells.points
        );
 
 
    Solver solver;
    */
    /*
    Solver::UnsteadryAnalysisResult unsteadry_analysis_result = solver.unsteadry_analysis(global_matrix.global_matrix,vec({1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}),0.01,1000);
    
    */

   //solver.mode_analysis_frequency(sp_mat(global_matrix.global_wave_matrix),global_matrix.global_nodal_matrix);

    /*
    std::cout<<"Mode Analysing ..."<<std::endl;
    Solver::ModeAnalysisResult mode_analysis_result= solver.mode_analysis(global_matrix.global_matrix,100);

    MeshUtils mesh_utils;
    */
    /*
    mesh_utils.write_mesh(mesh_points,mesh_cells,unsteadry_analysis_result.point_values,unsteadry_analysis_result.times,"output_datas");
    */

    /*
    for(int i=0;i<mode_analysis_result.display_mode_num;i++){
        mesh_utils.write_mesh(point_cells.points,point_cells.cells,mode_analysis_result.modes[i].point_values,mode_analysis_result.modes[i].times,"output_mode"+std::to_string(i+1));
    }
    */

    return 0;

}