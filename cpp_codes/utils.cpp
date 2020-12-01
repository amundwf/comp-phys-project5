#include "utils.hpp"
#include <iostream>
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

void writeGeneralMatrixToCSV(mat results, field<string> columnLabels, string filename, string directory){
    // columnLabels contains the labels of the columns, e.g. "t", "x", "y", "z" or "t", "L".
    // It should have the same amount of elements as 'results' has columns.
    ofstream ofile;
    string filePath = directory + filename;

    // Save matrix in CSV format with the column labels in the header:
    //results.save(csv_name("results.csv", header));
    results.save(csv_name(filePath, columnLabels));
}

void writeGeneralMatrixToCSV_noLabels(mat results, string filename, string directory){
    // columnLabels contains the labels of the columns, e.g. "t", "x", "y", "z" or "t", "L".
    // It should have the same amount of elements as 'results' has columns.
    ofstream ofile;
    string filePath = directory + filename;

    // Save matrix in CSV format with the column labels in the header:
    //results.save(csv_name("results.csv", header));
    results.save(filePath, csv_ascii);
}


vec ThomasAlgorithm(int n, vec u, double a, double b, bool verbose) {
    /* Returns solution to the Thomas algorith algorithm for a tri-
    diagonal matrix.

    Inputs:
        alpha: double, will give diagonal elements a, b 
        n: integer, number of grid points
        x_0: double, point 0
        x_n: double, point n+1
        h: double, step size in spacial dimension.
    Outputs:
        v: double, solution of size n+1 
    */
    
    vec G = u;
    vec B = vec(n+1);
    vec unew = vec(n+1);
    B[0] = b;

    if (verbose==true){
            cout << "\nForward substitution... "<< endl;
        }
    // Forward substitution:
    for (int i=1; i < n; i++) { // Start at index 2, or third row. up to n-1.
        double factorDiv = B[i-1];
        double factorMult = a;
        double factor = factorMult/factorDiv;

        B[i] = b - a*factor;
        G[i] = G[i] - G[i-1]*factor;

        if (verbose==true){
            cout << "B[i] is: "<< B[i] << endl;
            cout << "G[i] is: "<< G[i] << endl;
        }
    }
    // Need the last element of B[n] for backward sub.
    B[n] = b;

    // Just in case set the boundary condition maunally.
    unew[n] = 1.0; unew[0] = 0.0;
    G[n] = 1.0; G[0] = 0.0;

    if (verbose==true){
            cout << "\nBackward substitution... "<< endl;
        }
    // Backward substitution:
    for (int i=n-1; i > 0; i--) { // Start at row n-1 (index n-2), end at first row (row 1, index 0).
        double factorDiv = B[i+1];
        double factorMult = a;
        double factor = factorMult/factorDiv;
        // All upper diagonal elements gets eliminated.
        // Final element in the row:
        G[i] = G[i] - G[i+1]*factor; //in very first run i.e t=1 then G[i+1] = G[n]

        if (verbose==true){
            cout << "G[i] is: "<< G[i] << endl;
        }
    }
    if (verbose==true){
            cout << "\nScaling... "<< endl;
        }
    // Normalize the diagonal (divide all row i by b[i] for all rows) in order to get the 
    // solution for v:
    for (int i=1; i < n; i++) {
        //std::cout << "g[i]: " << g[i] << std::endl;
        //std::cout << "b[i]: " << b[i] << std::endl;
        unew[i] = G[i]/B[i];

        if (verbose==true){
            cout << "unew[i+1]: "<< unew[i+1] << endl;
        }
    }
    
    
    return unew;
}

void implicitScheme(int n, int tFinal, double tStep){
    // First we set initialise the new and old vectors
    // Set boundary conditions.
    // We have n+1 grid points so from 0 to L. starting at x_0 to x_n, which are both boundaries.
    vec u = vec(n+1, fill::zeros);
    vec unew = vec(n+1, fill::zeros);

    // Set the boundary conditions.
    double u_0 = 0.0;
    double u_n = 1.0;

    u(0) = unew(0) = u_0;
    u(n) = unew(n) = u_n;

    // Evaluate Delta x.
    double xStep = (u(n)- u(0)) / n;

    // Find Delta t and the number of time steps to get to tFinal.
    int tSteps = int(tFinal/ tStep);

    // Evaluate alpha , i.e Delta t / (Delta x * Delta x). 
    double alpha = tStep / (xStep*xStep);

    cout << "n is: " << n <<endl;
    cout << "xStep is: " << xStep <<endl;
    cout << "tFinal is: " << tFinal <<endl;
    cout << "tStep is: " << tStep <<endl;
    cout << "tSteps is: " << tSteps <<endl;
    cout << "alpha is: " << alpha << endl;
    
    // Set up the table to store solution at all time steps.
    mat results = mat(tSteps+1, n+1);

    // Add initial results to matrix.
    results(0, span(0,n)) = u.t();

    cout << "Initial u vec added to matrix" <<endl;
    // Time integration

    double b = 1 + 2*alpha;
    double a = -alpha;
    double hh = xStep*xStep;
    u = u*hh;
    cout << "Running Thomas algo for A^-1." << endl;
    
    bool verbose = true;
    for (int t = 1; t < 4; t++) {
        cout << "\nt is: " << t << endl;
        vec unew = ThomasAlgorithm(n, u, a, b, verbose);
        //  note that the boundaries are not changed.
        results(t, span(0,n)) = unew.t();
        u = unew;
    }  
    verbose = false;
    for (int t = 4; t <= tSteps; t++) {
        
        vec unew = ThomasAlgorithm(n, u, a, b, verbose);
        //  note that the boundaries are not changed.
        results(t, span(0,n)) = unew.t();
        u = unew;
    }  
    cout << "saving..." << endl;
    // Save the results.
        string directory = "../results/5c/";
        string filename = "Implicit_N=" + to_string(n) + "tSteps=" + to_string(tSteps) + ".csv";
        writeGeneralMatrixToCSV_noLabels(results, filename, directory); 
}

void explicitScheme(int n, int tFinal){
    // First we set initialise the new and old vectors
    // Set boundary conditions.
    // We have n+1 grid points so from 0 to L. starting at x_0 to x_n, which are both boundaries.
    vec u = vec(n+1, fill::zeros);
    vec unew = vec(n+1, fill::zeros);
    u(0) = unew(0) = 0.0;
    u(n) = unew(n) = 1.0;

    // Evaluate Delta x.
    double xStep = (u(n)- u(0)) / n;

    // Find Delta t and the number of time steps to get to tFinal.
    // Stability criteria, constrains t step. 
    double tStep = xStep*xStep/2;
    // Round upwards.
    int tSteps = ceil(tFinal/ tStep);  

    // Evaluate alpha , i.e Delta t / (Delta x * Delta x). 
    double alpha = tStep / (xStep*xStep);

    cout << "n is: " << n <<endl;
    cout << "xStep is: " << xStep <<endl;
    cout << "tFinal is: " << tFinal <<endl;
    cout << "tStep is: " << tStep <<endl;
    cout << "tSteps is: " << tSteps <<endl;
    
    // Set up the table to store solution at all time steps.
    mat results = mat(tSteps+1, n+1);

    // Add initial results to matrix.
    results(0, span(0,n)) = u.t();

    cout << "Initial u vec added to matrix" <<endl;
    // Time integration
    cout << "Running for all t..." << endl;
    for (int t = 1; t <= tSteps; t++) {
        for (int i = 1; i < n; i++) {
            // Discretized diff eq
            unew(i) = alpha * u(i-1) + (1 - 2*alpha) * u(i) + alpha * u(i+1);
        }
        //  note that the boundaries are not changed.
        results(t, span(0,n)) = unew.t();
        u = unew;
    }  
    cout << "saving..." << endl;
    // Save the results.
        string directory = "../results/5c/";
        string filename = "Explicit_N=" + to_string(n) + "tSteps=" + to_string(tSteps) + ".csv";
        writeGeneralMatrixToCSV_noLabels(results, filename, directory); 
}

void crankNicolsonScheme(int n, int tFinal, double tStep){
    // First we set initialise the new and old vectors
    // Set boundary conditions.
    // We have n+1 grid points so from 0 to L. starting at x_0 to x_n, which are both boundaries.
    vec u = vec(n+1, fill::zeros);
    vec r = vec(n+1, fill::zeros);
    vec unew = vec(n+1, fill::zeros);
    u(0) = unew(0) = 0.0;
    u(n) = unew(n) = 1.0;

    // Evaluate Delta x.
    double xStep = (u(n)- u(0)) / n;
    // Find Delta t and the number of time steps to get to tFinal.

    int tSteps = int(tFinal/ tStep);

    // Evaluate alpha , i.e Delta t / (Delta x * Delta x). and the diagonals. 
    double alpha = tStep / (xStep*xStep);
    double a = - alpha;
    double b = 2 + 2*alpha;
    double hh = xStep*xStep;

    // Need to scale by xStep*xStep.
    u = u*hh;

    cout << "n is: " << n <<endl;
    cout << "xStep is: " << xStep <<endl;
    cout << "tFinal is: " << tFinal <<endl;
    cout << "tStep is: " << tStep <<endl;
    cout << "tSteps is: " << tSteps <<endl;
    
    // Set up the table to store solution at all time steps.
    mat results = mat(tSteps+1, n+1);

    // Add initial results to matrix.
    results(0, span(0,n)) = u.t();

    bool verbose = false;
    // Time integration
    cout << "Running for all t..." << endl;
    for (int t = 1; t <= tSteps; t++) {
        // Calculate the right hand side for CN. 
        for (int i = 1; i < n; i++) {
        	r(i) = alpha*u(i-1) + (2 - 2*alpha)*u(i) + alpha*u(i+1);
        }
        r(0) = 0.0;
        r(n) = 1.0;

        // Run Thomas algo
        vec unew = ThomasAlgorithm(n, r, a, b, verbose);

        //  note that the boundaries are not changed.
        results(t, span(0,n)) = unew.t();
        u = unew;
    }  
    cout << "saving results..." << endl;
    // Save the results.
    string directory = "../results/5c/";
    string filename = "CrankNicolson_N=" + to_string(n) + "tSteps=" + to_string(tSteps) + ".csv";
    writeGeneralMatrixToCSV_noLabels(results, filename, directory); 
}

void diffusion1D(){
    int n = 10;
    int tFinal = 10; 
    double tStep = 0.005;
    explicitScheme(n, tFinal);
    implicitScheme(n, tFinal, tStep);
    crankNicolsonScheme(n, tFinal, tStep);
}