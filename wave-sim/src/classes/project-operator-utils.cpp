#include <project-operator.h>

void ProjectOperator::NewProject(std::string project_name){
    project_json["project_name"] = project_name;
    project_json["properties"]["sound_speed"] = 1;
    //TODO 球じゃなくて平面で指定できるようにする。
    project_json["properties"]["zero_points"] = {};

    UpdateProjectFile();
}

ProjectOperator::ProjectOperator(){}

ProjectOperator::ProjectOperator(std::string projectnameorfile){
    std::string project_filename;

    if(projectnameorfile.ends_with(".json")){
        project_filename = projectnameorfile;
    }else{
        project_filename = projectnameorfile + ".json";
    }

    ReadProjectFile(project_filename);    

    if(project_json["files"].contains("wave_matrix_file")){
        global_matrix.global_wave_matrix = LoadMatrixFile(
            project_json["files"]["wave_matrix_file"]
        );
    }
    if(project_json["files"].contains("nodal_matrix_file")){
        global_matrix.global_nodal_matrix = LoadMatrixFile(
            project_json["files"]["nodal_matrix_file"]
        );
    }
    if(project_json["files"].contains("global_matrix_file")){
        global_matrix.global_matrix = LoadMatrixFile(
            project_json["files"]["global_matrix_file"]
        );
    }
    if(project_json["files"].contains("corresponding_point_ids_file")){
        global_matrix.corresponding_point_ids = LoadIntVectorFile(
            project_json["files"]["corresponding_point_ids_file"]
        );
    }
    if(project_json["files"].contains("corresponding_point_id_indexes_file")){
        global_matrix.corresponding_point_id_indexes = LoadIntVectorFile(
            project_json["files"]["corresponding_point_id_indexes_file"]
        );
    }
    if(project_json["files"].contains("point_vector_file")){
        global_matrix.point_vector = LoadVectorFile(
            project_json["files"]["point_vector_file"]
        );
    }

}

void ProjectOperator::ReadProjectFile(std::string project_filename){

    std::ifstream stream(project_filename);

    if (stream.is_open()){
        // json::acceptがフォーマットチェック時にpositionを進めてしまうので、先頭に戻す
        stream.seekg(0, std::ios::beg);
        project_json =  vtknlohmann::json::parse(stream);

    }else{
        std::cerr<<"Fail to open Json file:"+project_filename << std::endl;
        if (!vtknlohmann::json::accept(stream)){
            std::cerr<<"Invalid Joson stle"<<std::endl;
        }
    }
}

void ProjectOperator::UpdateProjectFile(){

    std::string project_filename = project_json["project_name"].get<std::string>() + ".json";

    std::ofstream outputFile(project_filename);

    if(outputFile.is_open()){
        outputFile << project_json.dump(4);
        outputFile.close();

        std::cout<<"Project file saved as:"+ project_filename<<std::endl;
    }else{
        std::cerr<<"Unable to save project file:"+project_filename<<std::endl;
    }
}

vtknlohmann::json ProjectOperator::GetProject(){
    return project_json;
}


void ProjectOperator::ReadMeshFile(std::string mesh_filename){

    MeshUtils meshutils;
    MeshUtils::points_cells points_cells;

    project_json["mesh_filename"] = mesh_filename;

    if(mesh_filename.ends_with(".vtk")){
        std::cout<<"mesh file type  vtk"<<std::endl;
        points_cells = meshutils.read_vtk_file(mesh_filename);
    }else if(mesh_filename.ends_with(".vtm")){
        std::cout<<"mesh file type  vtm"<<std::endl;
        points_cells = meshutils.read_vtm_file(mesh_filename);
    }else if(mesh_filename.ends_with(".cgns")){
        std::cout<<"mesh file type  cgns"<<std::endl;
        points_cells = meshutils.read_cgns_file(mesh_filename);
    }

    std::cout<<"complite read mesh file"<<std::endl;
    std::cout<<project_json["properties"]["sound_speed"]<<std::endl;

    for(int i=0;i<points_cells.cells.size();i++){
        points_cells.cells[i].sound_speed = project_json["properties"]["sound_speed"].get<double>();
    }

    std::cout<<"Creating Global Matrix..."<<std::endl;

    GlobalMatrix global_matrix = GlobalMatrix(points_cells.cells,points_cells.points);


    std::string wave_matrix_filename = project_json["project_name"].get<std::string>() + "_wave_matrix.bin";
    SaveMatrixFile(global_matrix.global_wave_matrix,wave_matrix_filename,"wave_matrix_file");

    std::string nodal_matrix_filename = project_json["project_name"].get<std::string>() + "_nodal_matrix.bin";
    arma::mat global_nodal_matrix(global_matrix.global_nodal_matrix);
    SaveMatrixFile(global_nodal_matrix,nodal_matrix_filename,"nodal_matrix_file");

    std::string point_vector_filename = project_json["project_name"].get<std::string>() + "_pv.bin";
    SaveVectorFile(global_matrix.point_vector,point_vector_filename,"point_vector_file");

    std::string corresponding_point_ids_filename = project_json["project_name"].get<std::string>() + "_cpi.bin";
    SaveIntVectorFile(global_matrix.corresponding_point_ids,corresponding_point_ids_filename,"corresponding_point_ids_file");

    std::string corresponding_point_id_indexes_filename = project_json["project_name"].get<std::string>() + "_cpii.bin";
    SaveIntVectorFile(global_matrix.corresponding_point_id_indexes,corresponding_point_id_indexes_filename,"corresponding_point_id_indexes_file");

    UpdateProjectFile();
}

void ProjectOperator::SaveMatrixFile(arma::mat& matrix,std::string matrix_filename,std::string description){

    std::cout<<"Saving matrix "+description+" : "+matrix_filename+" ...";
    if(matrix.save(matrix_filename,arma::arma_binary)){
        project_json["files"][description] = matrix_filename;
        std::cout<<"SUCCESS"<<std::endl;
    }else{
        std::cerr<<"Fail to save matrix as:" + matrix_filename<<std::endl;
    }
}

void ProjectOperator::SaveCxMatrixFile(arma::cx_mat& matrix,std::string matrix_filename,std::string description){

    std::cout<<"Saving matrix "+description+" : "+matrix_filename+" ...";
    if(matrix.save(matrix_filename,arma::arma_binary)){
        project_json["files"][description] = matrix_filename;
        std::cout<<"SUCCESS"<<std::endl;
    }else{
        std::cerr<<"Fail to save matrix as:" + matrix_filename<<std::endl;
    }
}

void ProjectOperator::SaveVectorFile(arma::vec & vector,std::string vector_filename, std::string description){

    std::cout<<"Saving vector "+description+" : "+vector_filename+" ...";
    if(vector.save(vector_filename,arma::arma_binary)){
        project_json["files"][description] = vector_filename;
        std::cout<<"SUCCESS"<<std::endl;
    }else{
        std::cerr<<"Fail to save vector as:" + vector_filename<<std::endl;
    }
}

void ProjectOperator::SaveCxVectorFile(arma::cx_vec & vector,std::string vector_filename, std::string description){

    std::cout<<"Saving vector "+description+" : "+vector_filename+" ...";
    if(vector.save(vector_filename,arma::arma_binary)){
        project_json["files"][description] = vector_filename;
        std::cout<<"SUCCESS"<<std::endl;
    }else{
        std::cerr<<"Fail to save vector as:" + vector_filename<<std::endl;
    }
}

void ProjectOperator::SaveIntVectorFile(std::vector<int> vector,std::string vector_filename,std::string description){
    std::ofstream file(vector_filename, std::ios::binary);
    if (!file) {
        throw std::ios_base::failure("Unable to open file:"+vector_filename);
    }

    size_t size =vector.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));

    std::cout<<"Saving vector "+description+":"+vector_filename+" ...";
    if(file.write(reinterpret_cast<const char*>(vector.data()), size * sizeof(int))){
        project_json["files"][description] = vector_filename;
        std::cout<<"SUCCESS"<<std::endl;
    }else{
        std::cerr<<"Fail to save int vecotr as:" + vector_filename<<std::endl;
    }

    file.close();
}

arma::mat ProjectOperator::LoadMatrixFile(std::string matrix_filename){
    
    arma::mat matrix;
    std::cout<<"Loading matrix from:"+matrix_filename+" ...";
    if(matrix.load(matrix_filename)){
        std::cout<<"SUCCESS"<<std::endl;
    }else{
        std::cerr<<"Fail to load matrix"<<std::endl;    
    }
    return matrix;
}

arma::cx_mat ProjectOperator::LoadCxMatrixFile(std::string matrix_filename){
    
    arma::cx_mat matrix;
    std::cout<<"Loading matrix from:"+matrix_filename+" ...";
    if(matrix.load(matrix_filename)){
        std::cout<<"SUCCESS"<<std::endl;
    }else{
        std::cerr<<"Fail to load matrix"<<std::endl;    
    }
    return matrix;
}

arma::vec ProjectOperator::LoadVectorFile(std::string vector_filename){
    
    arma::vec vector;
    std::cout<<"Loading vector from:"+vector_filename + " ...";
    if(vector.load(vector_filename)){
        std::cout<<"SUCCESS"<<std::endl;
    }else{
        std::cerr<<"Fail to load vector"<<std::endl;    
    }
    return vector;
}

arma::cx_vec ProjectOperator::LoadCxVectorFile(std::string vector_filename){
    
    arma::cx_vec vector;
    std::cout<<"Loading vector from:"+vector_filename + " ...";
    if(vector.load(vector_filename)){
        std::cout<<"SUCCESS"<<std::endl;
    }else{
        std::cerr<<"Fail to load vector"<<std::endl;    
    }
    return vector;
}

std::vector<int> ProjectOperator::LoadIntVectorFile(std::string vector_filename){
   std::ifstream file(vector_filename, std::ios::binary);
    if (!file) {
        throw std::ios_base::failure("Unable to open "+vector_filename);
    }

    size_t size;
    file.read(reinterpret_cast<char*>(&size), sizeof(size));


    std::vector<int> data(size);
    std::cout<<"Loading Int Vector from:"+vector_filename+" ...";
    if(file.read(reinterpret_cast<char*>(data.data()), size * sizeof(int))){
        std::cout<<"SUCCESS"<<std::endl;
    }
    file.close();

    return data;
}

MeshUtils::points_cells ProjectOperator::LoadPointsCells(std::string mesh_filename){
    MeshUtils meshutils;
    MeshUtils::points_cells points_cells;

    project_json["mesh_filename"] = mesh_filename;

    if(mesh_filename.ends_with(".vtk")){
        std::cout<<"mesh file type  vtk"<<std::endl;
        points_cells = meshutils.read_vtk_file(mesh_filename);
    }else if(mesh_filename.ends_with(".vtm")){
        std::cout<<"mesh file type  vtm"<<std::endl;
        points_cells = meshutils.read_vtm_file(mesh_filename);
    }else if(mesh_filename.ends_with(".cgns")){
        std::cout<<"mesh file type  cgns"<<std::endl;
        points_cells = meshutils.read_cgns_file(mesh_filename);
    }

    return points_cells;
}


void ProjectOperator::Describe(std::string mesh_filename){
    MeshUtils mesh_utils;

    if(mesh_filename.ends_with(".vtk")){
        std::cout<<"mesh file type  vtk"<<std::endl;
        mesh_utils.describe_vtk_file(mesh_filename);
    }else if(mesh_filename.ends_with(".vtm")){
        std::cout<<"mesh file type  vtm"<<std::endl;
        mesh_utils.describe_vtm_file(mesh_filename);
    }else if(mesh_filename.ends_with(".cgns")){
        std::cout<<"mesh file type  cgns"<<std::endl;
        mesh_utils.describe_cgns_file(mesh_filename);
    }
}

void ProjectOperator::CleanAll(){
    if(project_json.contains("files")){

        std::vector<std::string> removed_keys;

        for(auto [key,file] :project_json["files"].items()){
            std::string filename = file.get<std::string>();
            std::cout<<"Removing :"+filename+" ...";
            if(std::filesystem::remove(filename)){
                removed_keys.push_back(key);
                std::cout<<"SUCCESS"<<std::endl;
            }else{
                std::cerr<<"Unabel to remove file"<<std::endl;
            }
        }

        for(auto removed_key : removed_keys){
            project_json["files"].erase(project_json["files"].find(removed_key));
        }

    }else{
        std::cout<<"Nothing to Clean"<<std::endl;
    }

    UpdateProjectFile();

}

void ProjectOperator::CopyProject(std::string origin_projectnamerofile){

    std::cout<<"LOADING ORIGIN PROJECT"<<std::endl;
    ProjectOperator origin_project(origin_projectnamerofile);

    vtknlohmann::json origin_project_json = origin_project.GetProject();

    project_json["properties"] = origin_project_json["properties"];

    std::string origin_project_name = origin_project_json["project_name"];
    std::string project_name = project_json["project_name"];

    for(auto [key,filename] : origin_project_json["files"].items()){
        std::string filename_str = filename;
        size_t pos = filename_str.find(origin_project_json["project_name"]);
        std::string origin_filename = filename_str;
        std::string new_filename = filename_str.replace(
            pos,
            origin_project_name.length(),
            project_name
        );

        std::cout<<"Copying File :"+origin_filename+" to "+new_filename;
        if(std::filesystem::copy_file(
            origin_filename, new_filename, 
            std::filesystem::copy_options::overwrite_existing)
        ){
            project_json["files"][key]= new_filename;
            std::cout<<"SUCCESS"<<std::endl;

        }else{
            std::cerr<<"Unable to copy file"<<std::endl;
        }

    }

    UpdateProjectFile();
}

void ProjectOperator::ShowHelp(){
    std::cout<<"This is help documentation"<<std::endl;
}