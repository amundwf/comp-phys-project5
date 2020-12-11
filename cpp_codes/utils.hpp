#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <armadillo>

// Saving matrices to file.

void writeGeneralMatrixToCSV_noLabels(arma::mat results, std::string filename, std::string directory);

void writeGeneralMatrixToCSV(arma::mat results, arma::field<std::string> columnLabels, std::string filename, std::string directory);

// Functions for the 1D diffusion equation.

double stabilityConditionExplicit_dt(double dx);

arma::vec ThomasAlgorithm(int n, arma::vec u, double a, double b, bool verbose);

void analytical_solution_1D(int n_x, double x_start, double x_end, double tFinal, double tStep, int N_sum);

void explicitScheme(int n, double x_start, double x_end, double tFinal,  bool verbose);

void explicitScheme_v2(int Nx, double tFinal, double tStep);

void implicitScheme(int n, double x_start, double x_end, double tFinal,  double tStep, bool verbose);

void implicitScheme_v2(int n, double tFinal, double tStep);

void crankNicolsonScheme(int n, double x_start, double x_end, double tFinal,  double tStep, bool verbose);

void crankNicolsonScheme_v2(int n, double tFinal, double tStep);

void diffusion1D();

// Solver for the 2D diffusion equation.

void diffusion2D();

int JacobiSolver(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol);

// Functions for after radioactive enrichment.

double Qtime(double time);

double Qdepth(int j, double dt);

// Temporary functions for before enrichment.

void lithosphere(bool enrichment);

int JacobiSolverLithosphere(int N, double dx, double dt, double t, double alpha, arma::mat &A, arma::mat &A_prev, double abstol, double beta, double gamma, double eta, bool enrichment);

// one dimensional schemes without inputs from terminal.

void run_5c();

#endif