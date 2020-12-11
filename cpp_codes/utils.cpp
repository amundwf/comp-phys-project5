#include "utils.hpp"
#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip> // Needed for set precision.
#include <math.h> // M_PI
#include "omp.h"
#define NUM_THREADS 8

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
    // Saves the arma::mat. ALthough not as a csv but as ascii. 
    // I could be better to change the name of this function to not make
    // it confusing. 
    ofstream ofile;
    string filePath = directory + filename;

    // Save matrix in CSV format without a header, in ascii format.
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
    vec B = vec(n);
    vec unew = vec(n);
    // Set the first element of B.
    // Not sure why we cant start at index 0.
    B[1] = b;

    // Printing if verbose.
    if (verbose==true){
        cout << "\nForward substitution... "<< endl;
    }
    // Forward substitution:
    // Start at index 1, up to n-2. Don't write over n-1 i.e the last element. 
    // Changed from 1 to 2 since no index 0. 
    for (int i=2; i < n-1; i++) { 
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
    unew[n-1] = 1.0; unew[0] = 0.0;
    
    // Printing if verbose.
    if (verbose==true){
        cout << "\nBackward substitution... "<< endl;
    }
    // Backward substitution:
    // Start at index n-2, end at index 1. Dont touch index 0 or n-1. 
    for (int i=n-2; i > 0; i--) { 
        double factorDiv = B[i];
        double factorMult = a;
        double factor = factorMult/factorDiv;
        // All upper diagonal elements gets eliminated.
        
        // In very first run i.e t=1 then G[i+1] = G[n], so we acces n but dont change it. 
        unew[i] = G[i]/factorDiv - unew[i+1]*factor;

        // Printing if verbose.
        if (verbose==true){
            cout << "unew[i] is: "<< unew[i] << endl;
        }
    }

    // Return the solution arma::vec unew.
    return unew;
}

void analytical_solution_1D(int n_x, double x_start, double x_end, double tFinal, double tStep, int N_sum){
    /* A funtion to calculate the analytical solution to the one dimensional problem.

    Inputs:
    n_x: integer, the number of x points.
    x_start: double, the start of x.
    x_end: double, the end of x.
    tFinal: double, the final time.
    tStep: double, dt or the time step.
    N_sum: integer, Number of terms to include in the sum of the analytical solution.
    The higher N_sum is, the better the approximation will be.

    Outputs: Void

    The function saves all relevant data to file.
    */

    // Set the number of nodes.
    omp_set_num_threads(NUM_THREADS);

    double L = x_end-x_start;
    // The constant k pre-calculates pi/L:
    double k = M_PI/L;
    double k2 = k*k;

    // Make the list of x values. n_x values between x=0 and x=1.
    double xStep = (x_end-x_start)/double(n_x-1);

    vec xList = zeros(n_x);
    for (int i=0; i<=n_x-1; i++){
        xList(i) = i*xStep;
    }
    // Make the list of t values.
    int tPoints = ceil(tFinal/tStep);
    vec tList = zeros(tPoints);
    for (int i=0; i<=tPoints-1; i++){
        tList(i) = i*tStep;
    }
    // Now calculate the values of u(x,t) and put them in a 2D array:
    mat v_xt_array = zeros(n_x, tPoints+1);
    mat u_xt_array = zeros(n_x, tPoints);

    cout << "\nRunning analytical solution..." << endl; 

    for (int i=0; i<=n_x-1; i++){ // For all x
        double x = xList(i);
        for (int j=0; j<=tPoints-1; j++){ // For all t
            double t = tList(j);
            double v_xt = 0;

            // Calculate the sum for v(x,t) (as in u(x,t)=v(x,t)-f(x)):
            double sum_element;
            int n;
            # pragma omp parallel for default(shared) private (n, sum_element) reduction(+:v_xt)
            for (int n=1; n<=N_sum; n++){
                sum_element = (1/n)*pow(-1,n)*sin(k*n*x)*exp(-k2*(n*n)*t);
                v_xt += sum_element;
                // The term '+ x' is from f(x)=-x/L with L=1, as in u(x,t)=v(x,t)-f(x).
            }

            v_xt *= 2/M_PI; // Multiply with 2/pi to get the correct value
            // Add u(x,t) to u_xt_array:
            v_xt_array(i,j) = v_xt;
            u_xt_array(i,j) = v_xt + (x/L);
        }
    }

    // Save the results:
    string directory = "../results/1D_diffusion/";
    writeGeneralMatrixToCSV_noLabels(v_xt_array.t(), "v_xt.csv", directory);
    writeGeneralMatrixToCSV_noLabels(u_xt_array.t(), "analytical_1D.csv", directory);
    // ^ u_xt_array is transposed to get the plots right.
    // Also save xList and tList for plotting:
    writeGeneralMatrixToCSV_noLabels(xList, "analytical_1D_xList.csv", directory);
    writeGeneralMatrixToCSV_noLabels(tList, "analytical_1D_tList.csv", directory);
}

void explicitScheme(int n, double x_start, double x_end, double tFinal, bool verbose ){
    /* A function which runs the explicit scheme (forward euler) for the 
    one dimensional case.

    Inputs:
    n: integer, the number of x points.
    x_start: double, the start of x.
    x_end: double, the end of x.
    tFinal: double, the final time.
    verbose: bool, if you want to print more to terminal.

    Outputs: Void

    The function saves all relevant data to file.
    */

    // First we set initialise the new and old vectors
    // Set boundary conditions.
    // We have n+1 grid points so from 0 to L. starting at x_0 to x_n, which are both boundaries.
    vec u = vec(n, fill::zeros);
    vec unew = vec(n, fill::zeros);
    u(0) = unew(0) = x_start;
    u(n-1) = unew(n-1) = x_end;

    // Evaluate Delta x.
    double xStep = (u(n-1) - u(0)) / (n-1);

    // Find Delta t and the number of time steps to get to tFinal.
    // Stability criteria, constrains t step. 
    double tStep = xStep*xStep/2;
    // Round upwards.
    int tPoints = ceil(tFinal/ tStep);  

    // Evaluate alpha , i.e Delta t / (Delta x * Delta x). 
    double alpha = tStep / (xStep*xStep);

    cout << "\nRunning Explicit Scheme ..." << endl;
    cout << "N is:      " << n <<endl;
    cout << "xStep is:  " << xStep <<endl;
    cout << "tFinal is: " << tFinal <<endl;
    cout << "tStep is:  " << tStep <<endl;
    cout << "tPoints is: " << tPoints <<endl;
    
    // Set up the table to store solution at all time steps.
    mat results = mat(tPoints+1, n);

    // Add initial results to matrix.
    results(0, span(0,n-1)) = u.t();

    // Time integration
    for (int t = 1; t <= tPoints; t++) {
        for (int i = 1; i < n-1; i++) {
            // Discretized diff eq
            unew(i) = alpha * u(i-1) + (1 - 2*alpha) * u(i) + alpha * u(i+1);
        }
        //  note that the boundaries are not changed.
        results(t, span(0,n-1)) = unew.t();
        u = unew;
    }  
    
    // Save the results.
    string directory = "../results/1D_diffusion/";
    string filename = "Explicit_N=" + to_string(n) + "_tPoints=" + to_string(tPoints) + ".csv";
    writeGeneralMatrixToCSV_noLabels(results, filename, directory);
}

void explicitScheme_v2(int n, double x_start, double x_end, double tFinal, double tStep){
    // This does the same as explicitScheme(), but now with tStep as an input.
    // n: Number of x points.
    double xStep = (x_end-x_start)/double(n-1);

    // Boundary conditions for u(x,t):
    vec u = vec(n, fill::zeros);
    vec unew = vec(n, fill::zeros);
    u(0) = unew(0) = 0.0;
    u(n-1) = unew(n-1) = 1.0;

    int tPoints = ceil(tFinal/tStep); // Number of time steps.

    // Evaluate alpha , i.e Delta t / (Delta x * Delta x).
    double alpha = tStep / (xStep*xStep);

    cout << "\nRunning Explicit Scheme ..." << endl;
    cout << "Nx is:      " << n <<endl;
    cout << "xStep is:  " << xStep <<endl;
    cout << "tFinal is: " << tFinal <<endl;
    cout << "tStep is:  " << tStep <<endl;
    cout << "tPoints is: " << tPoints <<endl;

    // Set up the table to store solution at all time steps.
    mat results = mat(tPoints+1, n);

    // Add initial results to matrix.
    results(0, span(0,n-1)) = u.t();

    // Time integration
    for (int t = 1; t <= tPoints; t++) {
        for (int i = 1; i < n-1; i++) {
            // Discretized diff eq
            unew(i) = alpha * u(i-1) + (1 - 2*alpha) * u(i) + alpha * u(i+1);
        }
        //  note that the boundaries are not changed.
        results(t, span(0,n-1)) = unew.t();
        u = unew;
    }

    cout << "results.n_cols: " << results.n_cols << endl;
    cout << "results.n_rows: " << results.n_rows << endl;
    results.print("results explicit:");

    // Save the results.
    string directory = "../results/1D_diffusion/";
    string filename = "explicit_1D.csv";
    writeGeneralMatrixToCSV_noLabels(results, filename, directory);
}

void implicitScheme(int n, double x_start, double x_end, double tFinal, double tStep, bool verbose){
    /* A function which runs the implicit scheme (backward euler) for the 
    one dimensional case.

    Inputs:
    n: integer, the number of x points.
    x_start: double, the start of x.
    x_end: double, the end of x.
    tFinal: double, the final time.
    tStep: double, dt the time step.
    verbose: bool, if you want to print more to terminal.
    
    Outputs: Void

    The function saves all relevant data to file.
    */

    // First we set initialise the new and old vectors
    // Set boundary conditions.
    // We have n+1 grid points so from 0 to L. starting at x_0 to x_n, which are both boundaries.
    vec u = vec(n, fill::zeros);
    vec unew = vec(n, fill::zeros);

    // Set the boundary conditions.
    double u_0 = x_start;
    double u_n = x_end;

    u(0) = unew(0) = u_0;
    u(n-1) = unew(n-1) = u_n;

    // Evaluate Delta x.
    double xStep = (u(n-1) - u(0)) / (n-1);

    // Find Delta t and the number of time steps to get to tFinal.
    int tPoints = ceil(tFinal/ tStep);

    // Evaluate alpha , i.e Delta t / (Delta x * Delta x). 
    double alpha = tStep / (xStep*xStep);

    cout << "\nRunning Implicit Scheme ..." << endl;
    cout << "N is:       " << n <<endl;
    cout << "xStep is:   " << xStep <<endl;
    cout << "tFinal is:  " << tFinal <<endl;
    cout << "tStep is:   " << tStep <<endl;
    cout << "tPoints is: " << tPoints <<endl;
    cout << "alpha is:   " << alpha << endl;
    
    // Set up the table to store solution at all time steps.
    mat results = mat(tPoints+1, n);

    // Add initial results to matrix.
    results(0, span(0,n-1)) = u.t();

    double b = 1 + 2*alpha;
    double a = -alpha;
    double hh = xStep*xStep;

    // Multiply by hh. 
    u = u*hh;
    
    for (int t = 1; t < 4; t++) {
        if (verbose==true){
            cout << "\nt is: " << t << endl;
        }
        vec unew = ThomasAlgorithm(n, u, a, b, verbose);
        //  note that the boundaries are not changed.
        results(t, span(0,n-1)) = unew.t();
        u = unew;
    }  
    verbose = false;
    for (int t = 4; t <= tPoints; t++) {
        
        vec unew = ThomasAlgorithm(n, u, a, b, verbose);
        //  note that the boundaries are not changed.
        results(t, span(0,n-1)) = unew.t();
        u = unew;
    }  
    
    // Save the results.
        string directory = "../results/1D_diffusion/";
        string filename = "Implicit_N=" + to_string(n) + "_tPoints=" + to_string(tPoints) + ".csv";
        writeGeneralMatrixToCSV_noLabels(results, filename, directory); 
}

void implicitScheme_v2(int n, double x_start, double x_end, double tFinal, double tStep){
    // Same as implicitScheme(), but saving the results file with a different name.
    vec u = vec(n, fill::zeros);
    vec unew = vec(n, fill::zeros);

    double xStep = (x_end-x_start)/double(n-1);

    // Set the boundary conditions.
    double u_0 = 0.0;
    double u_n = 1.0;
    u(0) = unew(0) = u_0;
    u(n-1) = unew(n-1) = u_n;

    // Find Delta t and the number of time steps to get to tFinal.
    int tPoints = ceil(tFinal/ tStep);

    // Evaluate alpha , i.e Delta t / (Delta x * Delta x). 
    double alpha = tStep / (xStep*xStep);

    cout << "\nRunning Implicit Scheme ..." << endl;
    cout << "N is:       " << n <<endl;
    cout << "xStep is:   " << xStep <<endl;
    cout << "tFinal is:  " << tFinal <<endl;
    cout << "tStep is:   " << tStep <<endl;
    cout << "tPoints is: " << tPoints <<endl;
    cout << "alpha is:   " << alpha << endl;
    
    // Set up the table to store solution at all time steps.
    mat results = mat(tPoints+1, n);

    // Add initial results to matrix.
    results(0, span(0,n-1)) = u.t();

    double b = 1 + 2*alpha;
    double a = -alpha;
    double hh = xStep*xStep;

    // Multiply by hh. 
    u = u*hh;

    bool verbose = false;
    for (int t = 1; t < 4; t++) {
        //cout << "\nt is: " << t << endl;

        vec unew = ThomasAlgorithm(n, u, a, b, verbose);
        //  note that the boundaries are not changed.

        results(t, span(0,n-1)) = unew.t();
        u = unew;
    }  

    for (int t = 4; t <= tPoints; t++) {
        vec unew = ThomasAlgorithm(n, u, a, b, verbose);
        // Note that the boundaries are not changed.
        results(t, span(0,n-1)) = unew.t();
        u = unew;
    }
    
    cout << "results.n_cols: " << results.n_cols << endl;
    cout << "results.n_rows: " << results.n_rows << endl;
    results.print("results implicit:");
    results(0,span::all).print("row idx = 0:");
    results(tPoints-1,span::all).print("row idx = tPoints-1:");
    results(tPoints,span::all).print("row idx = tPoints:");

    // Save the results.
    string directory = "../results/1D_diffusion/";
    string filename = "implicit_1D.csv";
    writeGeneralMatrixToCSV_noLabels(results, filename, directory);
}

void crankNicolsonScheme(int n, double x_start, double x_end, double tFinal, double tStep, bool verbose){
    /* A function which runs the Crank Nicolson for the 
    one dimensional case, which is a mix between forward and 
    backward schemes.

    Inputs:
    n: integer, the number of x points.
    x_start: double, the start of x.
    x_end: double, the end of x.
    tFinal: double, the final time.
    tStep: double, dt the time step.
    verbose: bool, if you want to print more to terminal.
    
    Outputs: Void

    The function saves all relevant data to file.
    */

    // First we set initialise the new and old vectors
    // Set boundary conditions.
    // We have n+1 grid points so from 0 to L. starting at x_0 to x_n, which are both boundaries.
    vec u = vec(n, fill::zeros);
    vec r = vec(n, fill::zeros);
    vec unew = vec(n+1, fill::zeros);
    u(0) = unew(0) = x_start;
    u(n-1) = unew(n-1) = x_end;

    // Evaluate Delta x.
    double xStep = (u(n-1)- u(0)) / (n-1);
    // Find Delta t and the number of time steps to get to tFinal.

    int tPoints = ceil(tFinal/ tStep);

    // Evaluate alpha , i.e Delta t / (Delta x * Delta x). and the diagonals. 
    double alpha = tStep / (xStep*xStep);
    double a = - alpha;
    double b = 2 + 2*alpha;
    double hh = xStep*xStep;

    // Need to scale by xStep*xStep.
    u = u*hh;

    cout << "\nRunning Crank Nicolson Scheme ..." << endl;
    cout << "Npoints is:   " << n <<endl;
    cout << "xStep is:     " << xStep <<endl;
    cout << "tFinal is:    " << tFinal <<endl;
    cout << "tStep is:     " << tStep <<endl;
    cout << "tPoints is:   " << tPoints <<endl;
    
    // Set up the table to store solution at all time steps.
    mat results = mat(tPoints+1, n);

    // Add initial results to matrix.
    results(0, span(0,n-1)) = u.t();

    verbose = false;
    // Time integration
    for (int t = 1; t <= tPoints; t++) {
        // Calculate the right hand side for CN. 
        for (int i = 1; i < n-1; i++) {
        	r(i) = alpha*u(i-1) + (2 - 2*alpha)*u(i) + alpha*u(i+1);
        }
        r(0) = 0.0;
        r(n-1) = 1.0;

        // Run Thomas algo
        vec unew = ThomasAlgorithm(n, r, a, b, verbose);

        // Add data to results. 
        results(t, span(0,n-1)) = unew.t();
        // Reset u for the next time step.
        u = unew;
    }  
    
    // Save the results.
    string directory = "../results/1D_diffusion/";
    string filename = "CrankNicolson_N=" + to_string(n) + "_tPoints=" + to_string(tPoints) + ".csv";
    writeGeneralMatrixToCSV_noLabels(results, filename, directory); 
}

void crankNicolsonScheme_v2(int n, double x_start, double x_end, double tFinal, double tStep){
    // This function does the same as crankNicolsonScheme(), but changes the saved
    // results file name.

    // We have n+1 grid points so from 0 to L. starting at x_0 to x_n, which are both boundaries.
    vec u = vec(n, fill::zeros);
    vec r = vec(n, fill::zeros);
    vec unew = vec(n+1, fill::zeros);
    u(0) = unew(0) = 0.0;
    u(n-1) = unew(n-1) = 1.0;

    // Evaluate Delta x.
    double xStep = (u(n-1)- u(0)) / (n-1);
    // Find Delta t and the number of time steps to get to tFinal.

    int tPoints = ceil(tFinal/tStep);

    // Evaluate alpha , i.e Delta t / (Delta x * Delta x). and the diagonals. 
    double alpha = tStep/(xStep*xStep);
    double a = - alpha;
    double b = 2 + 2*alpha;
    double hh = xStep*xStep;

    // Need to scale by xStep*xStep.
    u = u*hh;

    cout << "\nRunning Crank Nicolson scheme ..." << endl;
    cout << "Npoints is:   " << n <<endl;
    cout << "xStep is:     " << xStep <<endl;
    cout << "tFinal is:    " << tFinal <<endl;
    cout << "tStep is:     " << tStep <<endl;
    cout << "tPoints is:   " << tPoints <<endl;
    
    // Set up the table to store solution at all time steps.
    mat results = mat(tPoints+1, n);

    // Add initial results to matrix.
    results(0, span(0,n-1)) = u.t();

    bool verbose = false;
    // Time integration
    for (int t = 1; t <= tPoints; t++) {
        // Calculate the right hand side for CN. 
        for (int i = 1; i < n-1; i++) {
        	r(i) = alpha*u(i-1) + (2 - 2*alpha)*u(i) + alpha*u(i+1);
        }
        r(0) = 0.0;
        r(n-1) = 1.0;

        // Run Thomas algo
        vec unew = ThomasAlgorithm(n, r, a, b, verbose);

        // Add data to results. 
        results(t, span(0,n-1)) = unew.t();
        // Reset u for the next time step.
        u = unew;
    }  
    
    // Save the results.
    string directory = "../results/1D_diffusion/";
    string filename = "crankNicolson_1D.csv";
    writeGeneralMatrixToCSV_noLabels(results, filename, directory); 
}

void diffusion1D(){
    /* Runs all one dimensional schemes including the analytical solver.

    Inputs: None

    Outputs: Void

    All necessary data is written to file.

    */

    double tFinal; int Npoints; double dt; string verboseOrNot; bool verbose; double x_start; double x_end;
    string dtExplicitOrNot;

    cout << "Please enter x_start (double)" << endl;
    cin >> x_start;

    cout << "Please enter x_end (double)" << endl;
    cin >> x_end;

    cout << "Please enter Npoints (int)" << endl;
    cin >> Npoints;

    double dx = (x_end-x_start)/(Npoints-1);
    cout << "dx is: " << dx << endl;

    cout << "\nPlease enter tFinal (double)" << endl;
    cin >> tFinal;
    
    cout << "Do you want to use dt dictated by the explicit scheme? [Y/N]" << endl;
    cin >> dtExplicitOrNot;

    if (dtExplicitOrNot == "y"){
        // Set dt to the same as explicit
        dt = (dx*dx)/2;

    } else{
        cout << "Please enter dt for Implicit and Crank Nicolson schemes(double)" << endl;
        // Ask for dt instead.
        cin >> dt;
    }

    cout << "Do you want to print info to terminal? (Y/N)?";
    cin >> verboseOrNot;

    if (verboseOrNot == "y"){
        verbose = true;
    }
    else{
        verbose = false;
    }

    // The number which the analytic solution should go to.
    int N_sum = 10000;

    explicitScheme(Npoints, x_start, x_end, tFinal, verbose);
    implicitScheme(Npoints, x_start, x_end, tFinal, dt, verbose);
    crankNicolsonScheme(Npoints, x_start, x_end, tFinal, dt, verbose);
    analytical_solution_1D(Npoints, x_start, x_end, tFinal, dt, N_sum);
}



// These functions below do not run with physical constants.

void diffusion2D(){
    /* Runs the two dimensional (dimensionless) problem including the analytical solver.

    Inputs: None

    Outputs: Void

    All necessary data is written to file.

    */

    omp_set_num_threads(NUM_THREADS);
    cout << "The number of processors available = " << omp_get_num_procs ( ) << endl;

    double tFinal; int Npoints; double dt;
    
    cout << "Please enter Npoints (int)..." << endl;
    cin >> Npoints;

    double dx = 1.0/(Npoints-1);
    cout << "dx is: " << dx << endl;

    cout << "Please enter tFinal (can be double)..." << endl;
    cin >> tFinal;
    
    cout << "Please enter dt (double)..." << endl;
    cin >> dt;

    int Tpoints = int(tFinal / dt);
    double tolerance = 1.0e-14;
    double totalRelativeError = 0.0;
    mat A = zeros<mat>(Npoints,Npoints);
    mat A_prev = zeros<mat>(Npoints,Npoints);
    cube results = cube(Npoints, Npoints, Tpoints);
    mat A_analytic = zeros<mat>(Npoints,Npoints);
    cube resultsAnalytic = cube(Npoints, Npoints, Tpoints);

    // setting up an additional source term. for analytic comparison.
    // For t=0
    for(int i = 1; i < Npoints-1; i++){
        for(int j = 1; j < Npoints-1; j++){
            A(i,j) = sin(M_PI*dx*i)*sin(M_PI*dx*j);
        }
    }
    
    // Boundary Conditions -- all zeros
    for(int i=0; i < Npoints; i++){
        A(0,i) = 0.0; // Top of matrix
        A(Npoints-1, i) = 0.0; // Bottom of matrix.
        A(i,0) = 0.0; // Left side.
        A(i, Npoints-1) = 0.0; // Right side.
    }

    // Store initial conditions. 
    results(span::all, span::all, span(0)) = A;
    resultsAnalytic(span::all, span::all, span(0)) = A;

    // Loop over time.
    for( int t = 1; t < Tpoints; t++){
        double time = dt*t;
        A_prev = A;
        //cout << A << endl;
        int itcount = JacobiSolver(Npoints,dx,dt,A,A_prev,tolerance);

        // Store A in cube results.
        results( span::all, span::all, span(t)) = A;

        // Testing against exact solution
        double sum = 0.0;
        for(int i=0; i < Npoints; i++){
            for(int j=0; j < Npoints; j++){
                // I removed the mius sign in exact.
                A_analytic(i,j) = sin(M_PI*dx*i)*sin(M_PI*dx*j)*exp(-2*M_PI*M_PI*time);
                sum += fabs( (A(i,j) - A_analytic(i,j))) ;
            }
        }
        // Store analytic result.
        resultsAnalytic( span::all, span::all, span(t)) = A_analytic;
        //totalRelativeError += sum;
        cout << setprecision(5) << setiosflags(ios::scientific);
        cout << "Jacobi method with error " << sum/Npoints << " in " << itcount << " iterations" << endl;
    }
    // End time loop.

    // Relative error per point
    //cout << "The relative error per point is: " << totalRelativeError/Npoints << endl;

    ofstream ofile;
    string directory = "../results/2D_diffusion/";
    string filename =  "Tpoints=" + to_string(Tpoints)+ "_Npoints=" + to_string(Npoints) + ".txt";
    string filePath = directory + filename;
    results.save(filePath, raw_ascii);

    string filenameAnalytic = "Tpoints=" + to_string(Tpoints)+ "_Npoints=" + to_string(Npoints) + "_Analytic.txt";
    string filePathAnalytic = directory + filenameAnalytic;
    resultsAnalytic.save(filePathAnalytic, raw_ascii);
}

// Function for setting up the iterative Jacobi solver
int JacobiSolver(int N, double dx, double dt, mat &A, mat &A_prev, double abstol){
    /* Function for the iterative Jacobi solver. The function returns
    the iteration it converges at, or the maxiteration without convergence.

    Inputs:
    N: int, the number of x and y points. 
    dx: double, the x and y step.
    dt: double, the time step.
    A: arma::mat, the solution. At t=0, A = 0.0. After that it is changed each loop.
    A_prev: arma::mat, the solution from the previous time step. At t=0 A_prev = 0.0. 
    abstol: double, the stopping tolerance for the Jacobi solver.

    Outputs:
    k: int, the iteration at which the solver reaches convergence or if not it returns Maxiteration. 
    */

    int MaxIterations = 100000;
    double alpha = dt/(dx*dx);

    // Aold is the inital guess which starts as 1.0
    // everywhere. As the iterations go on it is changed.
    mat Aold = zeros<mat>(N, N); 
    for(int i=1;  i < N-1; i++){
        for(int j=1; j < N-1; j++){
            Aold(i,j) = 1.0;
        }
    }

    // I think the boundary conditions need to be given to Aold as well each time iteration.
    // Boundary Conditions set each time step to make sure.
    for(int i=0; i < N; i++){
        Aold(0,i) = 0.0; // Top of matrix
        Aold(N-1, i) = 0.0; // Bottom of matrix.
        Aold(i,0) = 0.0; // Left side.
        Aold(i, N-1) = 0.0; // Right side.
    }

    
    // Start the iterative solver
    for(int k=0; k < MaxIterations; k++){
        double sum = 0.0;
        // Declare i and j from omp.
        int i, j;
        // Start parallel task.
        # pragma omp parallel default(shared) private(i,j) reduction(+:sum)
        {
            # pragma omp for
            for(int i=1; i < N-1; i++){
                for(int j=1; j < N-1; j++){
                    A(i,j) = (1/(1 + 4*alpha))*( A_prev(i,j) + alpha*( Aold(i+1,j) + Aold(i,j+1) + 
                    Aold(i-1,j) + Aold(i,j-1) ) );
                }
            }

            // Sum the error at each location.
            // And make Aold = A for the next iteration. 
            
            // In the example code this went from 0 to N, which i dont want.
            for(int i = 1; i < N-1;i++){
                for(int j = 1; j < N-1;j++){
                    sum += fabs( Aold(i,j) - A(i,j) );
                    Aold(i,j) = A(i,j);
                }
            }
        } // End parallel task.

        // Does the error reach the stopping criteria
        if( (sum /= (N*N)) < abstol){
            return k;
        }
    }
    cerr << "Jacobi: Maximum Number of Interations Reached Without Convergence\n";
    return MaxIterations;
}


// Functions to run before the enrichment. 

void lithosphere(bool enrichment){
    /* Runs the two dimensional diffusion solver for the lithosphere problem.
    It does not take into account any heat production from enrichment of 
    radioactive elements.  

    Inputs: None

    Outputs: Void

    All necessary data is written to file.

    */

    double tFinal; double dt; double maxDepth; double dx;

    // Time is in Gy (since 1Gy ago) so how long do you want to run in Gy
    cout << "Starting from 1 Giga year ago, please enter the final time in Giga years (double)" << endl;
    cin >> tFinal;

    // What about the time step..
    cout << "Please enter a time step (dt) in Gy (double)" << endl;
    cin >> dt;
    cout << "You have chosen dt= " << dt*1e9 << " years" << endl;

    // How deep in km
    cout << "Please enter max depth in km" << endl;
    cin >> maxDepth;

    // What about the distance step..
    cout << "Please enter the distance step (dx) in km" << endl;
    cin >> dx;

    int Tpoints = int(tFinal / dt);
    int Npoints = int(maxDepth / dx);

    cout << "You have chosen Tpoints= " << Tpoints << " and Npoints= " << Npoints << endl; 

    omp_set_num_threads(NUM_THREADS);
    cout << "The number of processors available = " << omp_get_num_procs ( ) << endl;

    // Calc alpha.
    double alpha = dt/(dx*dx);

    // density of lithosphere 3.510 Kg/m3.
    double rho = 3.510;
    // Thermal conductivity k, 2.5 W/m/K.
    double k = 2.5;
    // Specific heat capacity cp, 1000 J/Kg/K.
    double cp = 1000;

    // Now these constant will be changed to be in Joules, 
    // km, Giga years (Gy), Kelvin and kilograms (Kg).
    /*
    rho = rho * 1e9; // Kg/km3. 
    k = k * 3.15e19;

    // cp is unchanged, units are fine.
    double beta = 1/(rho*cp);
    cout << "The constant beta is: " << beta << endl;

    double gamma = k*beta*alpha;
    cout << "The constant gamma is: " << gamma << endl;
    */
    
    cout << "alpha is: " << alpha << endl;
    double eta = 3.15e7/rho;
    cout << "The constant eta is: " << eta << endl;
    double gamma = k*eta*alpha;
    cout << "The constant gamma is: " << gamma << endl;


    double tolerance = 1.0e-14;
    mat A = zeros<mat>(Npoints,Npoints);
    mat A_prev = zeros<mat>(Npoints,Npoints);
    cube results = cube(Npoints, Npoints, Tpoints);

    A += 281.15;

    // Boundary Conditions in Kelvin. 
    for(int i=0; i < Npoints; i++){
        A(0,i) = 281.15; // Top of matrix
        A(Npoints-1, i) = 1573.15; // Bottom of matrix.
        A(i,0) = 281.15; // Left side.
        A(i, Npoints-1) = 281.15; // Right side.
    }

    // Store initial conditions. 
    results(span::all, span::all, span(0)) = A;

    // Loop over time.
    for( int t = 1; t < Tpoints; t++){
        double time = dt*t; // time should be in Gy. 
        A_prev = A;

        int itcount = JacobiSolverLithosphere(Npoints,dx,dt,t,alpha,A,A_prev,tolerance,gamma,eta,enrichment);

        // Store A in cube results.
        results( span::all, span::all, span(t)) = A;

        cout << "Jacobi method with error in " << itcount << " iterations" << endl;
    }
    // End time loop.
    ofstream ofile;
    string directory;
    if (enrichment == true){
        directory = "../results/2D_diffusion_after_enrichment/";
    } else{
        directory = "../results/2D_diffusion_before_enrichment/";
    }
    
    string filename =  "Tpoints=" + to_string(Tpoints)+ "_Npoints=" + to_string(Npoints) + "Lithosphere.txt";
    string filePath = directory + filename;
    results.save(filePath, raw_ascii);
}

int JacobiSolverLithosphere(int N, double dx, double dt, double t, double alpha, mat &A, mat &A_prev, double abstol, double gamma, double eta, bool enrichment ){
    /* Function for the iterative Jacobi solver. The function returns
    the iteration it converges at, or the maxiteration without convergence.

    Inputs:
    N: int, the number of x and y points. 
    dx: double, the x and y step.
    dt: double, the time step.
    A: arma::mat, the solution. At t=0, A = 0.0. After that it is changed each loop.
    A_prev: arma::mat, the solution from the previous time step. At t=0 A_prev = 0.0. 
    abstol: double, the stopping tolerance for the Jacobi solver.
    eta: double, constants from the diffusion equation merged into one. 

    Outputs:
    k: int, the iteration at which the solver reaches convergence or if not it returns Maxiteration. 
    */

    int MaxIterations = 100000;
    
    // Aold is the inital guess which starts as 1.0
    // everywhere. As the iterations go on it is changed.
    mat Aold = zeros<mat>(N, N); 
    for(int i=1;  i < N-1; i++){
        for(int j=1; j < N-1; j++){
            Aold(i,j) = 1.0;
        }
    }
    
    // Boundary Conditions set each time step to make sure.
    // Boundary conditions in kelvin. 8 degree celcius at top to 1300 degrees celcius at bottom.
    for(int i=0; i < N; i++){
        Aold(0,i) = 281.15; // Top of matrix
        Aold(N-1, i) = 1573.15 ; // Bottom of matrix.
        Aold(i,0) = 281.15; // Left side.
        Aold(i, N-1) = 281.15; // Right side.
    }
    

    // Calculate heat production in mantle.
    double Qt = Qtime(t*dt);

    // Start the iterative solver
    for(int k=0; k < MaxIterations; k++){
        double sum = 0.0;
        // Declare i and j from omp.
        int i, j;
        // Start parallel task.
        # pragma omp parallel default(shared) private(i,j) reduction(+:sum)
        {
            # pragma omp for
            for(int i=1; i < N-1; i++){
                // Calculate the heat production as a function of depth.
                double Qi = Qdepth(i, dx);

                for(int j=1; j < N-1; j++){
                    // If we are in the mantle (below 40 km) and enrichment is true, add heat production from enichment. 
                    if ( i*dx > 40 && enrichment ){
                        double Qtotal = Qi + Qt;                  
                        A(i,j) = (1/(1 + 4*gamma))*( A_prev(i,j) + eta*dt*Qtotal + gamma*(Aold(i+1,j) + Aold(i,j+1) + 
                                Aold(i-1,j) + Aold(i,j-1)) );
                    }else{
                        A(i,j) = (1/(1 + 4*gamma))*( A_prev(i,j) + beta*dt*Qi + gamma*(Aold(i+1,j) + Aold(i,j+1) + 
                                Aold(i-1,j) + Aold(i,j-1)) );
                    }  
                }
            }

            // Sum the error at each location.
            // And make Aold = A for the next iteration.
            for(int i = 1; i < N-1;i++){
                for(int j = 1; j < N-1;j++){
                    sum += fabs( Aold(i,j) - A(i,j) );
                    Aold(i,j) = A(i,j);
                }
            }
        } // End parallel task.

        // Does the error reach the stopping criteria
        if( (sum /= (N*N)) < abstol){
            return k;
        }
    }
    cerr << "Jacobi: Maximum Number of Interations Reached Without Convergence\n";
    return MaxIterations;
}

double Qtime(double time){
    // Returns Qt, the heat production after time in Gy.

    // At t=0 this is the heat production formula.
    // where U is Uranium, Th is Thorium and K
    // is Potassium.
    // +0.5 = Q(t=0) = 0.4*U(t=0) + 0.4*Th(t=0) 0.2*K(t=0)

    // Half-life equation. starting with n0 = 1. 
    double nU = pow(0.5, (time/4.47)); 
    double nTh = pow(0.5, (time/14.0));
    double nK = pow(0.5, (time/1.25));

    // Calc the sum with different ratios in the mantle. 
    double sum = 0.4*nU + 0.4*nTh + 0.2*nK;
    // Multiply this fraction (sum) by the heat production at t=0. 
    // At t=0 is was 0.5 micro Watts / m^3, now it is in
    // Joules / Gy / km^3
    double Qt = 0.5*sum;

    return Qt;
}

double Qdepth(int i, double dx){
    // Return the heat production as a function of depth, Qi or Qdepth. 

    // depth in km.
    double depth = dx*i;

    // Q units were micro Watts /m^3 now they are
    // Joules / Gy / km^3

    // If less than or equal to 20 km.
    if (depth <= 20){
        return 1.40;
    }

    // If between 20 and 40 km.
    if (depth > 20 && depth <= 40 ){
        return 0.35;
    }
    // Else return Qdepth if greater than 40 km.
    return 0.05;
}

void run_5c(){
    /* Runs the one dimensional schemes without inputs from the terminal.
    */

    double tFinal = 1; // Final time
    int Nx = 10; // Number of x points between 0 and L=1. In 5c: 10 or 100)
    double x_start = 0; double x_end = 1;
    double dx = (x_end-x_start)/double(Nx-1);
    double dt = 0.5*(dx*dx); // Explicit scheme stability condition.

    int N_sum = 10000; // The number of terms to include in the sum of the
    // analytical solution.

    // Get the results (saved in .csv files):

    analytical_solution_1D(Nx, x_start, x_end, tFinal, dt, N_sum);
    explicitScheme_v2(Nx, x_start, x_end, tFinal, dt);
    implicitScheme_v2(Nx, x_start, x_end, tFinal, dt);
    //implicitScheme(Nx, tFinal, dt, true);
    crankNicolsonScheme_v2(Nx, x_start, x_end, tFinal, dt);
}
