## ライブラリをビルド

### openblas

git clone OPENBLAS LIBRARY REPOSITORY
mkdir build
cd build 
cmake -DBUILD_SHARED_LIBS=ON ..
make
sudo make install

### lapack

wget LAPACK LIBRARY REPOSITORY
mkdir build
cd build
cmake -DCBLAS:BOOL=ON -DLAPACKE:BOOL=ON -DBUILD_TESTING=ON -DLAPACKE_WITH_TMG:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=ON  ../
make
make test
sudo make install

### arpack

git clone ARPACK-NG LIBRARY REPOSITORY
mkdir build
cd build
cmake -D EXAMPLES=ON -D ICB=ON -D EIGEN=ON ..
make
make all check
sudo make install

### superLU
git clone SUPERLU LIBRARY REPOSITORY
// armadillo library require superlu version 5
git checkout v5.3.0
mkdir build
cd build
cmake ..
make 
sudo make install

### armadillo
wget ARMADILLO LIBRARY REPOSITORY
mkdir build
cmake -DALLOW_OPENBLAS_MACOS=ON ..   
make
make test
sudo make install