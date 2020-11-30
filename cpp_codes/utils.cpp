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

/*
double* Thomas_algo(int n, double x_0, double x_np1, int a, int b) {
     Returns solution to the special algorithm.

    This function calculates v for the general case, i.e a 
    tri-diagonal matrix, where as opposed to the general case
    a and c are equal and all elements along the diagonal and 
    off-diagonals are identical. 
    Inputs:
        a, b: integers, diagonal elements
        n: interger, number of points
        x_0: double, point 0
        x_np1: double, point n+1
    Outputs:
        v: double, solution of size n+1 

    // Calculate h. Using long double for increased precision.
    double h = (x_np1-x_0)/(n+1);
    double hh = h*h;
    double *fList = new double[n+2]; double *g = new double[n+2]; 
    double *xList = new double[n+2]; double *v = new double[n+2];
    double *b_sub = new double[n+2]; double *g_sub = new double[n+2];

    // Create xList:
    xList[0] = x_0; // End points.
    xList[n+1] = x_np1;

    fList[0] = x_0; // End points.
    fList[n+1] = x_np1;
    for (int i=1; i<=n; i++) { xList[i] = x_0 + i*h; }    // Step size h between points.

    for (int i=0; i<=n+1; i++) {
        fList[i] = f(xList[i]);
        g[i] = hh * fList[i];
    }

    //clock_t start, finish; // declare start and final time
    //start = clock();

    // Forward substitution
    b_sub[1] = b;
    for (int i=2; i<=n+1; i++){
        b_sub[i] = (i+1)/i;
        g_sub[i] = g[i] + (g_sub[i-1]/b_sub[i-1]);
        } 
    v[n] = g_sub[n] / b_sub[n];

    // Backward substitution
    for (int i=n-1; i>0; i--){
        v[i] = (g_sub[i] + v[i+1])/b_sub[i];
        }
            
    //finish = clock();

    delete [] xList; delete [] fList; delete [] g;
    delete [] b_sub; delete [] g_sub;

    return v;//, ((finish - start)/CLOCKS_PER_SEC);
}
*/

void explicitScheme(int n){
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
    int tFinal = 10;
    // Stability criteria, constrains t step. 
    double tStep = xStep*xStep/2;
    int tSteps = int(tFinal/ tStep);

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
    }  
    cout << "saving..." << endl;
    // Save the results.
        string directory = "../results/5c/";
        string filename = "N=" + to_string(n) + "tSteps=" + to_string(tSteps) + ".csv";
        writeGeneralMatrixToCSV_noLabels(results, filename, directory); 
}

void task_5c(){
    int n = 100;
    explicitScheme(n);
}