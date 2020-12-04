#include <armadillo>
#include <iostream>
#include "utils.hpp"

using namespace std;
using namespace arma;

int main(){
    int dim; string lithosphereOrNot;
    
    cout << "Please enter if you want to run in 1D or 2D diffusion equation (int)" << endl;
    cin >> dim;

    if (dim == 1){
        diffusion1D();
    }
    
    if (dim == 2){
        cout << "Do you want to run the function with physical constant? [Y/N]" << endl;
        cin >> lithosphereOrNot;

        if (lithosphereOrNot == "y"){
            diffusion2DLithosphere();
        }
        else{
            diffusion2D();
        }

        
    }

    return 0;
}
