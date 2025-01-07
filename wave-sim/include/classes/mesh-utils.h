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


class MeshUtils{

    public:

    struct points_cells{
        std::vector<Point> points;
        std::vector<Cell> cells;
    };

    points_cells read_mesh(std::string vtk_filename);
    points_cells resd_vtm_file(std::string vtm_filename);
    void make_paraview_data_file(std::string output_dir,std::vector<double> times,std::vector<std::string> datafile_names);
    void write_mesh(std::vector<Point> points,std::vector<Cell> cells,std::vector<arma::vec> values,std::vector<double> times,std::string output_dir);

};