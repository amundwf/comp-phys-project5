#include <armadillo>
#include <iostream>
#include "utils.hpp"

using namespace std;
using namespace arma;

int main(){
    int n_x=100; int N_sum=12000;
    //double tFinal=1e-6; double tStep=tFinal*0.1/5;
    double tFinal=1.4e-7; double tStep=5e-10;
    analytical_solution_1D(n_x, tFinal, tStep, N_sum);

    /*
    int dim; string lithosphereOrNot;
    
    cout << "Please enter if you want to run in 1D or 2D diffusion equation (int)" << endl;
    cin >> dim;

    // If statement for the dimension of the problem..
    if (dim == 1){
        diffusion1D();
    }
    
    if (dim == 2){
        cout << "Do you want to run the function for the Lithosphere? [Y/N]" << endl;
        cin >> lithosphereOrNot;

        if (lithosphereOrNot == "y"){
            diffusion2DLithosphere();
        }
        else{
            diffusion2D();
        }
    }
    */

    return 0;
}
