#include <iostream>
#include <string>
#include <vector>

class Cell{
    public:
        std::vector<double> points;
        std::vector<std::vector<double>> jacobian;
        std::vector<std::vector<double>> element_matrix;

        double gauss_integral(){

    
            return 0.0;
        }
        double calcurate_weak_form_integration_term(){

            return 0.0;
        }

        void calcurate_jacobian(){

        }

        void make_element_matrix(){

        }

        Cell(std::vector<double> init_points){
            points = init_points;

            calcurate_jacobian();

            calcurate_weak_form_integration_term();

            make_element_matrix();

        }
};


void main(){

}