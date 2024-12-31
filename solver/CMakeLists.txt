cmake_minimum_required(VERSION 3.12)
project(SolverProject)

# C++標準の設定
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 実行ファイルの名前とソースコード
add_executable(solver solver.cpp)

# VTKライブラリの設定を見つける
find_package(VTK REQUIRED)


# VTKヘッダーとライブラリをリンク
target_link_libraries(solver PRIVATE ${VTK_LIBRARIES})
target_include_directories(solver PRIVATE ${VTK_INCLUDE_DIRS})

# Armadilloライブラリをリンク
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(solver PRIVATE ${ARMADILLO_LIBRARIES})