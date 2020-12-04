#include "utils.hpp"
#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip> // Needed for set precision. 
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
    vec B = vec(n);
    vec unew = vec(n);
    // Set the first element of B.
    B[0] = b;
    // Need the last element of B[n] for backward sub.
    B[n-1] = b;

    // Printing if verbose.
    if (verbose==true){
        cout << "\nForward substitution... "<< endl;
    }
    // Forward substitution:
    // Start at index 1, up to n-2. Don't touch n-1 i.e the last element. 
    for (int i=1; i < n-1; i++) { 
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
    G[n-1] = 1.0; G[0] = 0.0;

    // Printing if verbose.
    if (verbose==true){
        cout << "\nBackward substitution... "<< endl;
    }
    // Backward substitution:
    // Start at index n-2, end at index 1. Dont touch index 0 or n-1. 
    for (int i=n-2; i > 0; i--) { 
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
    for (int i=1; i < n-1; i++) {
        unew[i] = G[i]/B[i];

        if (verbose==true){
            cout << "\nScaling... "<< endl;
            cout << "unew: "<< unew[i+1] << endl;
        }
    }
    // Return the solution arma::vec unew.
    return unew;
}

void implicitScheme(int n, int tFinal, double tStep, bool verbose){
    // First we set initialise the new and old vectors
    // Set boundary conditions.
    // We have n+1 grid points so from 0 to L. starting at x_0 to x_n, which are both boundaries.
    vec u = vec(n, fill::zeros);
    vec unew = vec(n, fill::zeros);

    // Set the boundary conditions.
    double u_0 = 0.0;
    double u_n = 1.0;

    u(0) = unew(0) = u_0;
    u(n-1) = unew(n-1) = u_n;

    // Evaluate Delta x.
    double xStep = (u(n-1) - u(0)) / (n-1);

    // Find Delta t and the number of time steps to get to tFinal.
    int tPoints = int(tFinal/ tStep);

    // Evaluate alpha , i.e Delta t / (Delta x * Delta x). 
    double alpha = tStep / (xStep*xStep);

    cout << "\nRunning Explicit Scheme ..." << endl;
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

void explicitScheme(int n, int tFinal, bool verbose ){
    // First we set initialise the new and old vectors
    // Set boundary conditions.
    // We have n+1 grid points so from 0 to L. starting at x_0 to x_n, which are both boundaries.
    vec u = vec(n, fill::zeros);
    vec unew = vec(n, fill::zeros);
    u(0) = unew(0) = 0.0;
    u(n-1) = unew(n-1) = 1.0;

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

void crankNicolsonScheme(int n, int tFinal, double tStep, bool verbose){
    // First we set initialise the new and old vectors
    // Set boundary conditions.
    // We have n+1 grid points so from 0 to L. starting at x_0 to x_n, which are both boundaries.
    vec u = vec(n, fill::zeros);
    vec r = vec(n, fill::zeros);
    vec unew = vec(n+1, fill::zeros);
    u(0) = unew(0) = 0.0;
    u(n-1) = unew(n-1) = 1.0;

    // Evaluate Delta x.
    double xStep = (u(n-1)- u(0)) / (n-1);
    // Find Delta t and the number of time steps to get to tFinal.

    int tPoints = int(tFinal/ tStep);

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

void diffusion1D(){

    double tFinal; int Npoints; double dt; string verboseOrNot; bool verbose;
    
    cout << "Please enter Npoints (int)..." << endl;
    cin >> Npoints;

    double dx = 1.0/(Npoints-1);
    cout << "dx is: " << dx << endl;

    cout << "Please enter tFinal (can be double)..." << endl;
    cin >> tFinal;
    
    cout << "Please enter dt for Implicit and CN (double)..." << endl;
    cin >> dt;

    cout << "Do you want verbose? (Y/N)?";
    cin >> verboseOrNot;

    if (verboseOrNot == "y"){
        verbose = true;
    }
    else{
        verbose = false;
    }

    explicitScheme(Npoints, tFinal, verbose);
    implicitScheme(Npoints, tFinal, dt, verbose);
    crankNicolsonScheme(Npoints, tFinal, dt, verbose);
}

void diffusion2D(){
    omp_set_num_threads(NUM_THREADS);
    cout << "The number of processors available = " << omp_get_num_procs ( ) << endl;

    double tFinal; int Npoints; double dt; double maxDepth;

    // This was the previous version before using physical constants. 
    /*
    cout << "Please enter Npoints (int)..." << endl;
    cin >> Npoints;

    double dx = 1.0/(Npoints-1);
    cout << "dx is: " << dx << endl;

    cout << "Please enter tFinal (can be double)..." << endl;
    cin >> tFinal;
    
    cout << "Please enter dt (double)..." << endl;
    cin >> dt;
    */

    // Time is in Gy (since 1Gy ago) so how long do you want to run in Gy
    cout << "Starting from 1 Giga year ago, how long do you want to run for in Gy (double)" << endl;
    cin >> tFinal;

    // What about the time step..
    cout << "What time step in Gy (double)" << endl;
    cin << dt;
    cout << "You have chosen dt= " << dt*1e9 << "years" << endl;

    // How deep in km
    cout << "Please enter max depth in km. Note that there are changes to Q between 0 and 120 km" << endl;
    cin >> maxDepth;

    // What about the distance step..
    cout << "Please enter the distance step (dx) in km" << endl;
    cin >> dx;

    int Tpoints = int(tFinal / dt);
    int Npoints = int(maxDepth / dx);

    cout << "You have chosen Tpoints= " << Tpoints << " and Npoints= " << Npoints >> endl; 

    // density of lithosphere 3.510 Kg/m3.
    double rho = 3.510;
    // Thermal conductivity k, 2.5 W/m/C.
    double k = 2.5;
    // Specific heat capacity cp, 1000 J/Kg/K.
    double cp = 1000;

    // Now these constant will be changed to be in Joules, 
    // km, Giga years (Gy), Kelvin and kilograms (Kg).
    rho = rho * 1e9;
    k = k * 3.15e19;
    // cp is unchanged, units are fine.
    double beta = 1/(rho*cp);

    double ExactSolution;
    double tolerance = 1.0e-14;
    mat A = zeros<mat>(Npoints,Npoints);
    mat A_prev = zeros<mat>(Npoints,Npoints);
    cube results = cube(Npoints, Npoints, Tpoints);

    // setting up an additional source term. 
    // This must to the initial state by add heat at some spots.
    // For t=0
    /*
    for(int i = 1; i < Npoints-1; i++){
        for(int j = 1; j < Npoints-1; j++){
            A_prev(i,j) = -2.0*M_PI*M_PI*sin(M_PI*dx*i)*sin(M_PI*dx*j);
        }
    }
    */
    // Boundary Conditions -- all zeros
    for(int i=0; i < Npoints; i++){
        A(0,i) = 0.0; // Top of matrix
        A(Npoints-1, i) = 1.0; // Bottom of matrix.
        A(i,0) = 0.0; // Left side.
        A(i, Npoints-1) = 0.0; // Right side.
    }

    // Store initial conditions. 
    results(span::all, span::all, span(0)) = A;

    // Loop over time.
    for( int t = 1; t < Tpoints; t++){
        double time = dt*t; // time should be in Gy. 
        A_prev = A;
        
        // Calculate Q, the heat production at time t.
        // Qtotal = Q(t) + Q(j). The depth part is added in JacobiSolver.
        double Qt = heatProduction(time);

        int itcount = JacobiSolver(Npoints,dx,dt,A,A_prev,tolerance,Qt,k,beta);

        // Store A in cube results.
        results( span::all, span::all, span(t)) = A;

        // Testing against exact solution
        double sum = 0.0;
        for(int i=0; i < Npoints; i++){
            for(int j=0; j < Npoints; j++){
                ExactSolution = -sin(M_PI*dx*i)*sin(M_PI*dx*j)*exp(-2*M_PI*M_PI*time);
                sum += fabs((A(i,j) - ExactSolution));
            }
        }

        cout << setprecision(5) << setiosflags(ios::scientific);
        cout << "Jacobi method with error " << sum/Npoints << " in " << itcount << " iterations" << endl;
    }
    // End time loop.
    ofstream ofile;
    string directory = "../results/2D_diffusion/";
    string filename =  "Tpoints=" + to_string(Tpoints)+ "_Npoints=" + to_string(Npoints) + ".txt";
    string filePath = directory + filename;
    results.save(filePath, raw_ascii);
}

// Function for setting up the iterative Jacobi solver
int JacobiSolver(int N, double dx, double dt, mat &A, mat &A_prev, double abstol, double Qt, double k, double beta){
    /* Return the iteration at which the Jacobi solver completes.

    Inputs:
    A, matrix. A at t=0 is A(Npoints,Npoints, fill::zeros). After that it is changed each loop.
    A_prev, matrix. A_prev at t=0 is a initial state and at t>0 it is equal to A from the 
    previous step.

    */

    // Check if A is all zeros 
    //cout << A_prev.is_zero() << endl;

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
        Aold(N-1, i) = 1.0; // Bottom of matrix.
        Aold(i,0) = 0.0; // Left side.
        Aold(i, N-1) = 0.0; // Right side.
    }

    
    // Start the iterative solver
    for(int k=0; k < MaxIterations; k++){
        double sum = 0.0;
        // Declare i and j from omp.
        int i, j;
        // Start parallel task.
        # pragma omp parallel default(shared) private(i,j, Qtotal) reduction(+:sum)
        {
            # pragma omp for
            for(int i=1; i < N-1; i++){
                for(int j=1; j < N-1; j++){
                    // Both Qt and Qdepth are in the correct units. 
                    // Joules / Gy / km^3
                    double Qtotal = Qdepth(j, dx) + Qt;

                    // Without physical constants. 
                    //A(i,j) = (1/(1 + 4*alpha))*( A_prev(i,j) + alpha*( Aold(i+1,j) + Aold(i,j+1) + 
                    //Aold(i-1,j) + Aold(i,j-1) ) );

                    // With physical constants. 
                    A(i,j) = (1/(1 + 4*alpha*beta*k))*( A_prev(i,j) + beta*(dt*Qtotal + k*alpha*( Aold(i+1,j) + Aold(i,j+1) + 
                    Aold(i-1,j) + Aold(i,j-1) ) ) );
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

double heatProduction(double time){
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
    double Qt = 0.5*3.15e19*sum;

    return Qt;
}

double Qdepth(int j, double dx){
    // Return the heat production as a function of depth, Qj or Qdepth. 

    // depth in km.
    double depth = dx*j;

    // Q units were micro Watts /m^3 now they are
    // Joules / Gy / km^3
    if (depth <= 40 && depth>= 20){
        return 0.35 * 3.15e19;
    }
    if (depth > 40){
        return 0.05 * 3.15e19;
    }
    if (depth < 20){
        return 1.40 * 3.15e19;
    }
}