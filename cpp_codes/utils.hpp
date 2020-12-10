#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <armadillo>

// Saving matrices to file.

void writeGeneralMatrixToCSV_noLabels(arma::mat results, std::string filename, std::string directory);

void writeGeneralMatrixToCSV(arma::mat results, arma::field<std::string> columnLabels, std::string filename, std::string directory);

// Functions for the 1D diffusion equation.

void diffusion1D();

double stabilityConditionExplicit_dt(double dx);

arma::vec ThomasAlgorithm(int n, arma::vec u, double a, double b, bool verbose);

void analytical_solution_1D(int n_x, double tFinal, double tStep, int N_sum);

void explicitScheme(int n, int tFinal, double x_start, double x_end, bool verbose);

void explicitScheme_v2(int Nx, double tFinal, double tStep);

void implicitScheme(int n, int tFinal, double x_start, double x_end, double tStep, bool verbose);

void implicitScheme_v2(int n, double tFinal, double tStep);

void crankNicolsonScheme(int n, int tFinal, double x_start, double x_end, double tStep, bool verbose);

void crankNicolsonScheme_v2(int n, double tFinal, double tStep);


// Solver for the 2D diffusion equation.

void diffusion2D();

int JacobiSolver(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol);

// Functions for after radioactive enrichment.

void diffusion2DAfterEnrichment();

int JacobiSolverAfterEnrichment(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol, double Qt, double k, double beta);

double Qtime(double time);

double Qdepth(int j, double dt);

// Temporary functions for before enrichment.

void diffusion2DBeforeEnrichment();

int JacobiSolverBeforeEnrichment(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol, double eta);

// one dimensional schemes without inputs from terminal.

void run_5c();

#endif