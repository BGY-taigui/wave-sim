#include <string>
#include <vector>
#include <armadillo>
#include <filesystem>
#include <mesh-class.h>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkQuad.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkCGNSReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>


#include <mesh-utils.h>

using namespace arma;

MeshUtils::points_cells MeshUtils::read_cgns_file(std::string cgns_filename){
    points_cells points_cells;

    std::cout<<"-------------------"<<std::endl;

     // CGNSファイルのリーダー
    auto reader = vtkSmartPointer<vtkCGNSReader>::New();
    reader->SetFileName(cgns_filename.c_str());
    reader->Update();

    // マルチブロックデータセットを取得
    auto multiBlockDataSet = reader->GetOutput();
    if (!multiBlockDataSet) {
        std::cerr << "Error: Unable to read the CGNS file." << std::endl;
    }

    std::cout << "Number of blocks: " << multiBlockDataSet->GetNumberOfBlocks() << std::endl;

    // 最初のブロックを処理（例として）
    auto block = multiBlockDataSet->GetBlock(0);
    if (!block) {
        std::cerr << "Error: Block 0 is empty." << std::endl;
    }

    for (unsigned int i = 0; i < multiBlockDataSet->GetNumberOfBlocks(); ++i) {
        auto block = multiBlockDataSet->GetBlock(i);
        if (!block) {
            std::cout << "Block " << i << " is empty." << std::endl;
            continue;
        }

        // ブロックの型を出力
        std::cout << "Block " << i << " class name: " << block->GetClassName() << std::endl;
        
        if (auto nestedMultiBlock = vtkMultiBlockDataSet::SafeDownCast(block)) {
            std::cout<<"this is vtkm dataset "<<std::endl;

            for (unsigned int j = 0; j < nestedMultiBlock->GetNumberOfBlocks(); ++j) {
                auto nestedBlock = nestedMultiBlock->GetBlock(j);
                if (!nestedBlock) {
                    std::cout << "Nested block " << j << " is empty." << std::endl;
                    continue;
                }

                auto vtk_dataset = vtkDataSet::SafeDownCast(nestedBlock);

                // ネストされたブロックの型を出力
                std::cout << "Nested block " << j << " class name: " << nestedBlock->GetClassName() << std::endl;

                //pointとCellの数を表示

                std::cout << "Number of points: " << vtk_dataset->GetNumberOfPoints() << std::endl;
                std::cout << "Number of cells: " << vtk_dataset->GetNumberOfCells() << std::endl;

                
                // 各点の座標を表示
                for (vtkIdType i = 0; i < vtk_dataset->GetNumberOfPoints(); ++i) {
                    double p[3];
                    vtk_dataset->GetPoint(i, p);
                    //TODO これ消す
                    //std::cout<<p[0]<<","<<p[1]<<","<<p[2]<<std::endl;

                    Point point(vec(p,sizeof(p)/sizeof(p[0])),i);
                    points_cells.points.push_back(point);

                }

                // 各セル（要素）の情報を表示

                int cell_id_count =0;
                for (vtkIdType i = 0; i < vtk_dataset->GetNumberOfCells(); ++i) {
                    vtkCell* cell = vtk_dataset->GetCell(i);
                    if(cell->GetCellType() == 12){
                        
                        std::vector<int> point_ids(8);
                        for(int num=0;num<cell->GetNumberOfPoints();num++){
                            point_ids[num] = cell->GetPointId(num);
                        }

                        cell_id_count +=1;

                        Cell cell(point_ids,cell_id_count);

                        points_cells.cells.push_back(cell);

                    }

                }

            }
        }

    }


    return points_cells;
}


MeshUtils::points_cells MeshUtils::read_vtm_file(std::string vtm_filename){

    points_cells points_cells;

    // VTMファイルのリーダー
    auto reader = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
    reader->SetFileName(vtm_filename.c_str());
    reader->Update();


    std::cout<<"-------------------"<<std::endl;


    // マルチブロックデータセットを取得
    //auto multiBlockDataSet = reader->GetOutput();

    auto multiBlockDataSet = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());


    if (!multiBlockDataSet) {
        std::cerr << "Error: Unable to read the VTM file." << std::endl;
    }else{
        std::cout << "Successfully loaded :"<<vtm_filename<<std::endl;
    }

    std::cout << "Number of blocks: " << multiBlockDataSet->GetNumberOfBlocks() << std::endl;

    //とりあえずブロックは一つの場合だけ処理する
    //TODO 複数ブロックにも対応させるか?
    if(multiBlockDataSet->GetNumberOfBlocks() != 1){
        std::cout << "WARNING::VTM file contain more than 1 blocks" << std::endl;
        std::cout << "Processing Only One Block" << std::endl;
    }

    // 各ブロックを処理
    for (unsigned int i = 0; i < multiBlockDataSet->GetNumberOfBlocks(); ++i) {
        auto block = multiBlockDataSet->GetBlock(i);
        if (block) {
            std::cout << "Block " << i << ": " << block->GetClassName() << std::endl;
        } else {
            std::cout << "Block " << i << " is empty." << std::endl;
        }
    }

    auto block = multiBlockDataSet->GetBlock(0);
    if (block) {

      // ブロックをvtkDataSetとして扱う
        auto dataSet = vtkDataSet::SafeDownCast(block);
        if (!dataSet) {
            std::cerr << "Block is not a vtkDataSet." << std::endl;
        }else{
            //TODO　全体処理をここに入れる
        }

        // Pointsを取得
        auto points = dataSet->GetPoints();
        if (points) {
            std::cout << "Number of points: " << points->GetNumberOfPoints() << std::endl;
            for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
                double p[3];
                points->GetPoint(i, p);
                std::cout << "Point " << i << ": (" << p[0] << ", " << p[1] << ", " << p[2] << ")" << std::endl;

                Point point(vec(p,sizeof(p)/sizeof(p[0])),i); 
                points_cells.points.push_back(point);
            }
        }

        // Cellsを取得
        vtkIdType numberOfCells = dataSet->GetNumberOfCells();
        std::cout << "Number of cells: " << numberOfCells << std::endl;
        for (vtkIdType i = 0; i < numberOfCells; ++i) {
            auto cell = dataSet->GetCell(i);
            std::cout << "Cell " << i << " has " << cell->GetNumberOfPoints() << " points." << std::endl;

            if(cell->GetCellType() == 12){
                
                std::vector<int> point_ids(8);
                for(int num=0;num<cell->GetNumberOfPoints();num++){
                    point_ids[num] = cell->GetPointId(num);
                }

            Cell cell(point_ids,i);

            points_cells.cells.push_back(cell);

        }


        }
    } else {
        std::cout << "Block is empty." << std::endl;
    }

    return points_cells;

}

MeshUtils::points_cells MeshUtils::read_mesh(std::string vtk_filename){

    MeshUtils::points_cells points_cells;
    points_cells.points = {};
    points_cells.cells = {};
    

    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(vtk_filename.c_str());
    std::cout<<"loading vtk file"<<std::endl;
    reader->Update();

    vtkUnstructuredGrid* mesh = reader->GetOutput();

    // メッシュの基本情報を表示
    std::cout << "Successfully loaded: " << vtk_filename << std::endl;
    std::cout << "Number of points: " << mesh->GetNumberOfPoints() << std::endl;
    std::cout << "Number of cells: " << mesh->GetNumberOfCells() << std::endl;



    // 各点の座標を表示
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
        double p[3];
        mesh->GetPoint(i, p);
        //TODO これ消す
        //std::cout<<p[0]<<","<<p[1]<<","<<p[2]<<std::endl;

        Point point(vec(p,sizeof(p)/sizeof(p[0])),i);
        points_cells.points.push_back(point);

    }

    // 各セル（要素）の情報を表示

    int cell_id_count =0;
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i) {
        vtkCell* cell = mesh->GetCell(i);
        if(cell->GetCellType() == 12){
            
            std::vector<int> point_ids(8);
            for(int num=0;num<cell->GetNumberOfPoints();num++){
                point_ids[num] = cell->GetPointId(num);
            }

            cell_id_count +=1;

            Cell cell(point_ids,cell_id_count);

            points_cells.cells.push_back(cell);

        }


        /*
        std::cout << "  Cell " << i << " (type " << cell->GetCellType() << "): ";
        for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); ++j) {
            std::cout << cell->GetPointId(j) << " ";
        }
        std::cout << std::endl;
*/
    }


    return points_cells;

}


void MeshUtils::make_paraview_data_file(std::string output_dir,std::vector<double> times,std::vector<std::string> datafile_names){    // PVDファイルのパス
    std::string pvd_path = output_dir + "/output_datas.pvd";

    // ファイルを開く
    std::ofstream pvd_file(pvd_path);
    if (!pvd_file.is_open()) {
        std::cerr << "Error: Unable to create PVD file: " << pvd_path << std::endl;
        return;
    }

    // PVDファイルのヘッダー部分
    pvd_file << R"(<?xml version="1.0"?>)" << "\n";
    pvd_file << R"(<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="none">)" << "\n";
    pvd_file << "  <Collection>\n";

    // 各タイムステップに対応するDataSetを記述
    for (size_t i = 0; i < times.size(); ++i) {
        std::string filename = datafile_names[i];
        pvd_file << "    <DataSet timestep=\"" << times[i] 
                << "\" part=\"0\" file=\"" << filename << "\"/>\n";
    }

    // PVDファイルのフッター部分
    pvd_file << "  </Collection>\n";
    pvd_file << "</VTKFile>\n";

    pvd_file.close();

    std::cout << "PVD file created: " << pvd_path << std::endl;
}


void MeshUtils::write_mesh(std::vector<Point> points,std::vector<Cell> cells,std::vector<vec> values,std::vector<double> times,std::string output_dir){

    if (!std::__fs::filesystem::exists(output_dir)) {
        std::__fs::filesystem::create_directories(output_dir);
    }


    // vtkPointsオブジェクトを作成して座標を追加
    auto vtk_points = vtkSmartPointer<vtkPoints>::New();
    for (const auto& point : points) {
        vtk_points->InsertNextPoint(point.vector(0),point.vector(1),point.vector(2));
    }

    // vtkCellArrayオブジェクトを作成してセル情報を追加
    auto vtk_cells = vtkSmartPointer<vtkCellArray>::New();
    for (const auto& cell : cells) {
        vtk_cells->InsertNextCell(cell.point_ids.size());
        for (int pointID : cell.point_ids) {
            vtk_cells->InsertCellPoint(pointID);
        }
    }

    // vtkIntArrayを使ってPoint IDを追加
    auto point_id_array = vtkSmartPointer<vtkIntArray>::New();
    point_id_array->SetName("PointIDs");
    for (auto point : points) {
        point_id_array->InsertNextValue(point.id);
    }

    // vtkIntArrayを使ってCell IDを追加
    auto cell_id_array = vtkSmartPointer<vtkIntArray>::New();
    cell_id_array->SetName("CellIDs");
    for (auto cell : cells) {
        cell_id_array->InsertNextValue(cell.id);
    }

    for(int i=0;i<times.size();i++){

        auto data_array = vtkSmartPointer<vtkDoubleArray>::New();
        data_array->SetName("Pressure (Pa)");
        data_array->SetNumberOfComponents(1); 
        data_array->SetNumberOfTuples(vtk_points->GetNumberOfPoints());

        for(int point_num=0;point_num<points.size();point_num++){
            data_array->SetValue(points[point_num].id, values[i](point_num));
        }

        // vtkUnstructuredGridオブジェクトを作成
        auto unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        unstructured_grid->SetPoints(vtk_points);
        unstructured_grid->SetCells(VTK_HEXAHEDRON, vtk_cells); // VTK_QUADは四角形セルのタイプ
        unstructured_grid->GetPointData()->AddArray(point_id_array);
        unstructured_grid->GetCellData()->AddArray(cell_id_array);

        unstructured_grid->GetPointData()->SetScalars(data_array);

        // ファイルに書き出し

        std::string filename_at_the_time = output_dir +"/" + std::to_string(times[i]) + ".vtu";

        auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(filename_at_the_time.c_str()); // 出力ファイル名
        writer->SetInputData(unstructured_grid);
        writer->Write();
    }

std::vector<std::string> output_filenames;
for(auto time : times){
    output_filenames.push_back(std::to_string(time)+".vtu");
}
make_paraview_data_file(output_dir,times,output_filenames);

}
