#include <armadillo>
#include <iostream>
#include "utils.hpp"

using namespace std;
using namespace arma;

int main(){
    int dim;
    
    cout << "Please enter if you want to run in 1D or 2D diffusion equation (int)" << endl;
    cin >> dim;

    if (dim == 1){
        diffusion1D();
    }
    
    if (dim == 2){
        diffusion2D();
    }

    return 0;
}
