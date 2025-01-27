#include <project-operator.h>

void ProjectOperator::ReadBoundaryCondition(){

    if(!project_json.contains("boundary_condition")){

        project_json["boundary_condition"] = {};
        UpdateProjectFile();
        std::cout<<"Nothing Boundary condition to load "<<std::endl;

        return ;
    }

    global_matrix.get_single_global_matrix(true);
    std::string global_matrix_filename = project_json["project_name"].get<std::string>() + "_global_matrix.bin";
    SaveMatrixFile(global_matrix.global_matrix,global_matrix_filename,"global_matrix_file");

    UpdateProjectFile();


    MeshUtils::points_cells points_cells = LoadPointsCells(
        project_json["mesh_filename"].get<std::string>()
    );


    for (auto boundary_condition : project_json["boundary_condition"]){
        std::string boundary_condition_type = boundary_condition["type"].get<std::string>();


        if(boundary_condition_type == "zero_sphere"){

            std::cout<<"Boundary Type : zero_sphere";

            arma::vec vector(3);
            vector(0) = boundary_condition["vector"][0].get<double>();
            vector(1) = boundary_condition["vector"][1].get<double>();
            vector(2) = boundary_condition["vector"][2].get<double>();
            double radius = boundary_condition["radius"].get<double>();

            std::cout<<" Vector x:"<<vector(0)<<" y:"<<vector(1)<<" z:"<<vector(2)<<std::endl;

            global_matrix.define_condittion_zero_point(
                vector,radius,
                points_cells.points
            );
        }
    }

    SaveMatrixFile(
        global_matrix.global_matrix,
        project_json["files"]["global_matrix_file"],
        "global_matrix"
    );


    UpdateProjectFile();
}

void ProjectOperator::RunModeAnalysis(){

    if(global_matrix.global_matrix.n_cols==0){
        global_matrix.get_single_global_matrix(true);

        std::string global_matrix_filename = project_json["project_name"].get<std::string>() + "_global_matrix.bin";
        SaveMatrixFile(global_matrix.global_matrix,global_matrix_filename,"global_matrix_file");

        UpdateProjectFile();
    }  

    Solver solver;
    //TODOファイル名を_file.binにする
    Solver::SortedEigs sorted_eigs = solver.compute_eigs(global_matrix.global_matrix);
    std::string sorted_eign_vec_filename = project_json["project_name"].get<std::string>() + "_sorted_eign_vec.bin";
    std::string sorted_eign_val_filename = project_json["project_name"].get<std::string>() + "_sorted_eign_val.bin";
    SaveCxVectorFile(sorted_eigs.eign_values,sorted_eign_val_filename,"sorted_eign_values");
    SaveCxMatrixFile(sorted_eigs.eign_vectors,sorted_eign_vec_filename,"sorted_eign_vectors");

    UpdateProjectFile();
}




void ProjectOperator::OutputTimeSeriesVTK(){

    if(!project_json["files"].contains("global_matrix_file")){
        std::cerr<<"Calcurate Global Matrix first"<<std::endl;
        return ;
    }
    //TODO 固有値、ベクトルファイルがないときもreturnする処理

    bool auto_conf = false;

    if(!project_json["analysis_confs"].contains("mode_analysis_type")){
        project_json["analysis_confs"]["mode_analysis_type"]="mode_num";
        auto_conf = true;
    }

    if(!project_json["analysis_confs"].contains("mode_analysis_num_range")){
        project_json["analysis_confs"]["mode_analysis_num_range"]={0,10};
        auto_conf = true;
    }

    if(!project_json["analysis_confs"].contains("mode_analysis_freq_range")){
        project_json["analysis_confs"]["mode_analysis_freq_range"]={0,10};
        auto_conf = true;
    }
    if(!project_json["analysis_confs"].contains("mode_analysis_display_timestep_num")){
        project_json["analysis_confs"]["mode_analysis_display_timestep_num"]=100;
        auto_conf = true;
    }
    
    if(auto_conf){UpdateProjectFile();}

    if(!project_json.contains("output_dir")){
        std::cerr<<"Define output_dir ("+project_json["project_name"].get<std::string>()+".json)"<<std::endl;
        return ;
    }

    arma::cx_vec eign_val = LoadCxVectorFile(project_json["files"]["sorted_eign_values"]);

    arma::cx_mat eign_vec = LoadCxMatrixFile(project_json["files"]["sorted_eign_vectors"]);

    MeshUtils::points_cells points_cells = LoadPointsCells(project_json["mesh_filename"]);

    Solver solver;
    std::vector<Solver::EachModeResult> result;

    if(project_json["analysis_confs"]["mode_analysis_type"].get<std::string>() == "mode_num"){

        std::cout<<"mode analysis num range"<<std::endl;
        result = solver.output_time_series_mode_num(
            eign_val,eign_vec,
            project_json["analysis_confs"]["mode_analysis_num_range"][0].get<int>(),
            project_json["analysis_confs"]["mode_analysis_num_range"][1].get<int>(),
            project_json["analysis_confs"]["mode_analysis_display_timestep_num"].get<int>()
        );

    }else if(project_json["analysis_confs"]["mode_analysis_type"] == "freq_range"){

        std::cout<<"mode analysis freq range"<<std::endl;
        result = solver.output_time_series_freq_range(
            eign_val,eign_vec,
            project_json["analysis_confs"]["mode_analysis_freq_range"][0].get<double>(),
            project_json["analysis_confs"]["mode_analysis_freq_range"][1].get<double>(),
            project_json["analysis_confs"]["mode_analysis_display_timestep_num"].get<int>()
        );


    }else{
        std::cerr<<"Invalid Mode Analysis Type :"+project_json["analysis_confs"]["mode_analysis_type"].get<std::string>()<<std::endl;
    }

    MeshUtils mesh_utils;

    int mode_counter_display =0;
    for(auto each_mode_result : result){
        mode_counter_display++;
        std::cout<<"Saving Mode :"<<mode_counter_display<<"/"<<result.size()<<" ";

        mesh_utils.write_mesh(
            points_cells.points,points_cells.cells,
            each_mode_result.point_values,
            each_mode_result.times,
            project_json["output_dir"].get<std::string>()+"/mode"+std::to_string(mode_counter_display)
            );
    }


    UpdateProjectFile();
}


void ProjectOperator::ShowRelevantModes(){

    if(!project_json["files"].contains("global_matrix_file")){
        std::cerr<<"Calcurate Global Matrix first"<<std::endl;
        return ;
    }
    //TODO 固有値、ベクトルファイルがないときもreturnする処理

    bool auto_conf = false;

    if(!project_json["analysis_confs"].contains("mode_analysis_type")){
        project_json["analysis_confs"]["mode_analysis_type"]="mode_num";
        auto_conf = true;
    }

    if(!project_json["analysis_confs"].contains("mode_analysis_num_range")){
        project_json["analysis_confs"]["mode_analysis_num_range"]={0,10};
        auto_conf = true;
    }

    if(!project_json["analysis_confs"].contains("mode_analysis_freq_range")){
        project_json["analysis_confs"]["mode_analysis_freq_range"]={0,10};
        auto_conf = true;
    }

    if(!project_json["analysis_confs"].contains("mode_analysis_display_timestep_num")){
        project_json["analysis_confs"]["mode_analysis_display_timestep_num"]=100;
        auto_conf = true;
    }
    
    if(auto_conf){UpdateProjectFile();}

    arma::cx_vec eign_val = LoadCxVectorFile(project_json["files"]["sorted_eign_values"]);

    arma::cx_mat eign_vec = LoadCxMatrixFile(project_json["files"]["sorted_eign_vectors"]);


    Solver solver;

    if(project_json["analysis_confs"]["mode_analysis_type"].get<std::string>() == "mode_num"){

        std::cout<<"mode analysis num range"<<std::endl;
        solver.output_time_series_mode_num(
            eign_val,eign_vec,
            project_json["analysis_confs"]["mode_analysis_num_range"][0].get<int>(),
            project_json["analysis_confs"]["mode_analysis_num_range"][1].get<int>(),
            project_json["analysis_confs"]["mode_analysis_display_timestep_num"].get<int>()
        );

    }else if(project_json["analysis_confs"]["mode_analysis_type"] == "freq_range"){

        std::cout<<"mode analysis freq range"<<std::endl;
        solver.output_time_series_freq_range(
            eign_val,eign_vec,
            project_json["analysis_confs"]["mode_analysis_freq_range"][0].get<double>(),
            project_json["analysis_confs"]["mode_analysis_freq_range"][1].get<double>(),
            project_json["analysis_confs"]["mode_analysis_display_timestep_num"].get<int>()
        );


    }else{
        std::cerr<<"Invalid Mode Analysis Type :"+project_json["analysis_confs"]["mode_analysis_type"].get<std::string>()<<std::endl;
    }

}