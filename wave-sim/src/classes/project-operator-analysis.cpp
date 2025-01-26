#include <project-operator.h>


void ProjectOperator::RunModeAnalysis(){

    //TODO ModeAnalysisの設定を読み込む処理を追加

    if(global_matrix.global_matrix.n_cols==0){
        global_matrix.get_single_global_matrix(true);

        std::string global_matrix_filename = project_json["project_name"].get<std::string>() + "_global_matrix.bin";
        SaveMatrixFile(global_matrix.global_matrix,global_matrix_filename,"global_matrix_file");

        UpdateProjectFile();
    }  

    Solver solver;
    Solver::SortedEigs sorted_eigs = solver.compute_eigs(global_matrix.global_matrix);
    std::string sorted_eign_vec_filename = project_json["project_name"].get<std::string>() + "_sorted_eign_vec.bin";
    std::string sorted_eign_val_filename = project_json["project_name"].get<std::string>() + "_sorted_eign_val.bin";
    SaveCxVectorFile(sorted_eigs.eign_values,sorted_eign_val_filename,"sorted_eign_values");
    SaveCxMatrixFile(sorted_eigs.eign_vectors,sorted_eign_vec_filename,"sorted_eign_vectors");

    UpdateProjectFile();
}

void ProjectOperator::OutputTimeSeriesVTK(){

    if(!project_json.contains("output_dir")){
        std::cerr<<"Define output_dir ("+project_json["project_name"].get<std::string>()+".json)"<<std::endl;
    }

    MeshUtils::points_cells points_cells = LoadPointsCells(project_json["mesh_filename"]);

    //mesh_utils.write_mesh(points_cells.points,points_cells.cells,mode_analysis_result.modes[i].point_values,mode_analysis_result.modes[i].times,output_dir);
}