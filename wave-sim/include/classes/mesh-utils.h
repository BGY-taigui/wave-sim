#pragma once

#include <string>
#include <vector>
#include <armadillo>
#include <filesystem>

#include <mesh-class.h>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkQuad.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkXMLMultiBlockDataReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>


class MeshUtils{

    public:

    struct points_cells{
        std::vector<Point> points;
        std::vector<Cell> cells;
    };

    points_cells read_vtk_file(std::string vtk_filename);
    points_cells read_vtm_file(std::string vtm_filename);
    points_cells read_cgns_file(std::string cgns_filename);

    void describe_vtk_file(std::string vtk_filename);
    void describe_vtm_file(std::string vtm_filename);
    void describe_cgns_file(std::string cgns_filename);

    void make_paraview_data_file(std::string output_dir,std::vector<double> times,std::vector<std::string> datafile_names);
    void write_mesh(std::vector<Point> points,std::vector<Cell> cells,std::vector<arma::vec> values,std::vector<double> times,std::string output_dir);

};