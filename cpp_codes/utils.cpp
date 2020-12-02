#include "utils.hpp"
#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip> // Needed for set precision. 

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
        a: double, the off diagonal element.
        b: double, the diagonal element.
        n: integer, number of grid points. There are actually n+1 grid points. 
        u: arma::vec, the input vector(n+1).
        verbose: bool, if you want to print out lots of stuff. 
    Outputs:
        unew: arma::vec, solution size(n+1)
    */
    
    // I could use u but G is more similar to the previous code.
    vec G = u;
    vec B = vec(n+1);
    vec unew = vec(n+1);
    // Set the first element of B.
    B[0] = b;
    // Need the last element of B[n] for backward sub.
    B[n] = b;

    // Printing if verbose.
    if (verbose==true){
        cout << "\nForward substitution... "<< endl;
    }
    // Forward substitution:
    // Start at index 2, up to n-1. Don't touch n. 
    for (int i=1; i < n; i++) { 
        double factorDiv = B[i-1];
        double factorMult = a;
        double factor = factorMult/factorDiv;

        B[i] = b - a*factor;
        G[i] = G[i] - G[i-1]*factor;

        // Printing if verbose.
        if (verbose==true){
            cout << "B[i] is: "<< B[i] << endl;
            cout << "G[i] is: "<< G[i] << endl;
        }
    }
    
    // Just in case set the boundary condition maunally.
    unew[n] = 1.0; unew[0] = 0.0;
    G[n] = 1.0; G[0] = 0.0;

    // Printing if verbose.
    if (verbose==true){
        cout << "\nBackward substitution... "<< endl;
    }
    // Backward substitution:
    // Start at index n-1, end at index 1. Dont touch index 0 or n. 
    for (int i=n-1; i > 0; i--) { 
        double factorDiv = B[i+1];
        double factorMult = a;
        double factor = factorMult/factorDiv;
        // All upper diagonal elements gets eliminated.
        
        // In very first run i.e t=1 then G[i+1] = G[n], so we acces n but dont change it. 
        G[i] = G[i] - G[i+1]*factor;

        // Printing if verbose.
        if (verbose==true){
            cout << "G[i] is: "<< G[i] << endl;
        }
    }
    // Normalize the diagonal (divide all row i by b[i] for all rows) in order to get the 
    // solution for v:
    for (int i=1; i < n; i++) {
        unew[i] = G[i]/B[i];

        if (verbose==true){
            cout << "\nScaling... "<< endl;
            cout << "unew: "<< unew[i+1] << endl;
        }
    }
    // Return the solution arma::vec unew.
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
    double xStep = (u(n) - u(0)) / n;

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
        string directory = "../results/1D_diffusion/";
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
    double xStep = (u(n) - u(0)) / n;

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
        string directory = "../results/1D_diffusion/";
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

        // Add data to results. 
        results(t, span(0,n)) = unew.t();
        // Reset u for the next time step.
        u = unew;
    }  
    cout << "saving results..." << endl;
    // Save the results.
    string directory = "../results/1D_diffusion/";
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

void diffusion2D(){
    int Npoints = 40;
    int Tpoints = 40;
    double ExactSolution;
    double dx = 1.0/(Npoints-1);
    double dt = 0.25*dx*dx;
    double tFinal = Tpoints * dt;
    double tolerance = 1.0e-14;
    mat A = zeros<mat>(Npoints,Npoints);
    mat q = zeros<mat>(Npoints,Npoints);
    cube results = cube(Npoints, Npoints, Tpoints);

    // setting up an additional source term. 
    // This must to the initial state by add heat at some spots.
    // For t=0
    for(int i = 0; i < Npoints; i++){
        for(int j = 0; j < Npoints; j++){
            q(i,j) = -2.0*M_PI*M_PI*sin(M_PI*dx*i)*sin(M_PI*dx*j);
        }
    }
    // Store initial conditions. 
    results(span::all, span::all, span(0)) = q;

    // Loop over time.
    for( int t = 1; t < Tpoints; t++){

        int itcount = JacobiSolver(Npoints,dx,dt,A,q,tolerance);

        // Store A in cube results.
        results( span::all, span::all, span(t)) = A;

        // Testing against exact solution
        double sum = 0.0;
        for(int i = 0; i < Npoints; i++){
            for(int j=0;j < Npoints; j++){
                ExactSolution = -sin(M_PI*dx*i)*sin(M_PI*dx*j)*exp(-2*M_PI*M_PI*t);
                
                sum += fabs((A(i,j) - ExactSolution));
            }
        }

        cout << setprecision(5) << setiosflags(ios::scientific);
        cout << "Jacobi method with error " << sum/Npoints << " in " << itcount << " iterations" << endl;
    }
    ofstream ofile;
    string directory = "../results/2D_diffusion/";
    string filename =  "Tpoints=" + to_string(Tpoints)+ "_Npoints=" + to_string(Npoints) + ".txt";
    string filePath = directory + filename;
    results.save(filePath, raw_ascii);
}

// Function for setting up the iterative Jacobi solver
int JacobiSolver(int N, double dx, double dt, mat &A, mat &q, double abstol){
    // A only has the boundary conditions, else it is zero. 
    // q is the previous time step solution. 

    int MaxIterations = 100000;
    // Aold is the inital guess. 
    mat Aold = zeros<mat>(N, N);
    
    double alpha = dt/(dx*dx);
    
    // This is the initial guess. 
    for(int i=1;  i < N-1; i++){
        for(int j=1; j < N-1; j++){
            Aold(i,j) = 1.0;
        }
    }
    
    // Boundary Conditions -- all zeros
    for(int i=0; i < N; i++){
    A(0,i) = 0.0; // Top of matrix
    A(N-1,i) = 0.0; // Bottom of matrix.
    A(i,0) = 0.0; // Left side.
    A(i,N-1) = 1.0; // Right side.
    }
    // Start the iterative solver
    for(int k = 0; k < MaxIterations; k++){
        for(int i=1; i < N-1; i++){
            for(int j=1; j < N-1; j++){
                A(i,j) = dt*q(i,j) + Aold(i,j) +
                alpha*(Aold(i+1,j) + Aold(i,j+1) - 4.0*Aold(i,j) + 
                Aold(i-1,j) + Aold(i,j-1));
            }
        }

        // Sum the error at each location.
        // And make Aold = A. 
        double sum = 0.0;
        for(int i = 0; i < N;i++){
            for(int j = 0; j < N;j++){
                sum += (Aold(i,j)-A(i,j))*(Aold(i,j)-A(i,j));
                Aold(i,j) = A(i,j);
            }
        }

        // Does the error reach the stopping criteria
        if(sqrt (sum) <abstol){
            return k;
        }
    }
    // A should go to the next time step as q. 
    for(int i = 0; i < N;i++){
            for(int j = 0; j < N;j++){
                A(i,j) = q(i,j);
            }
        }
    cerr << "Jacobi: Maximum Number of Interations Reached Without Convergence\n";
    return MaxIterations;
}