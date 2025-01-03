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


class MeshUtils{

    public:

    void read_mesh(){}

    void make_paraview_data_file(std::string output_dir,std::vector<double> times,std::vector<std::string> datafile_names);
    void write_mesh(std::vector<Point> points,std::vector<Cell> cells,std::vector<arma::vec> values,std::vector<double> times);

};