#include <armadillo>
#include <iostream>
#include "utils.hpp"

using namespace std;
using namespace arma;

int main(){
    /*
    run_5c();

    */
    int dim; string lithosphereOrNot; bool enrichment; string EnrichmentOrNot;
    
    cout << "Please enter if you want to run in 1D or 2D diffusion equation (int)" << endl;
    cin >> dim;

    // If statement for the dimension of the problem..
    if (dim == 1){
        diffusion1D();
    }
    
    if (dim == 2){
        cout << "Do you want to run the Lithosphere problem [Y] or Unit Test [N]? [Y/N]" << endl;
        cin >> lithosphereOrNot;

        if (lithosphereOrNot == "y"){
            
            cout << "Do you want to run after the mantle has been enriched? [Y/N]" << endl;
            cin >> EnrichmentOrNot;

            if (EnrichmentOrNot == "y"){
                enrichment = true;
            }
            else{
                enrichment = false;
            }
            lithosphere(enrichment);
        }
        else{
            diffusion2D();
        }
    }
    

    return 0;
}
