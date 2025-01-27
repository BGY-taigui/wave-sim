#include <iostream>
#include <string>
#include <filesystem>

#include <armadillo>

#include <vtkUnstructuredGrid.h>
#include <vtkNlohmann/json.hpp>

#include <global-matrix.h>
#include <mesh-utils.h>
#include <solver.h>

class ProjectOperator{
    private:

        vtknlohmann::json project_json;    
        GlobalMatrix global_matrix;
        

    public:

        void ReadProjectFile(std::string project_filename);
        vtknlohmann::json GetProject();
        void UpdateProjectFile();
        void NewProject(std::string project_name);

        void ReadMeshFile(std::string mesh_filename);
        
        //TODO constつける
        void SaveMatrixFile(arma::mat& matrix,std::string matrix_filename,std::string description);
        void SaveCxMatrixFile(arma::cx_mat& matrix,std::string matrix_filename,std::string description);
        void SaveVectorFile(arma::vec& vector,std::string vector_filename,std::string description);
        void SaveCxVectorFile(arma::cx_vec& vector,std::string vector_filename,std::string description);
        void SaveIntVectorFile(std::vector<int>,std::string vector_filename,std::string description);

        MeshUtils::points_cells LoadPointsCells(std::string mesh_filename);

        arma::mat LoadMatrixFile(std::string matrix_filename);
        arma::cx_mat LoadCxMatrixFile(std::string matrix_filename);
        arma::vec LoadVectorFile(std::string vector_filename);
        arma::cx_vec LoadCxVectorFile(std::string vector_filename);
        std::vector<int> LoadIntVectorFile(std::string vectro_filename);

        void OutputTimeSeriesVTK();
        void ShowRelevantModes();

        void ReadBoundaryCondition();
        void RunModeAnalysis();
        void RunUnsteadryAnalysis();
        void Rendering();

        void Describe(std::string mesh_filename);
        void CleanAll();
        void CopyProject(std::string origin_projectnameorfile);

        void ShowHelp();

        ProjectOperator();
        ProjectOperator(std::string projectnameorfile);
};

