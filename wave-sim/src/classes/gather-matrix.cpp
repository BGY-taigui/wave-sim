#include <vector>
#include <armadillo>

#include <gather-matrix.h>

using namespace arma;


GatherMatrix::GatherMatrix(std::vector<int> element_point_ids,std::vector<int>& corresponding_point_ids):
            L(8,corresponding_point_ids.size()) 
        {
            for(int i=0;i<corresponding_point_ids.size();i++){
                for(int j=0;j<element_point_ids.size();j++){
                    if(element_point_ids[j]==corresponding_point_ids[i]){
                        L(j,i) = 1;
        }}}} 