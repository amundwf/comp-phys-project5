#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <armadillo>


//double* Thomas_algo(int n, double x_0, double x_np1, int a, int b)

void writeGeneralMatrixToCSV_noLabels(arma::mat results, std::string filename, std::string directory);

void writeGeneralMatrixToCSV(arma::mat results, arma::field<std::string> columnLabels, std::string filename, std::string directory);

arma::vec ThomasAlgorithm(int n, arma::vec u, double a, double b, bool verbose);

void analytical_solution_1D(int n_x, double tFinal, double tStep, int N_sum);

void explicitScheme(int n, int tFinal, bool verbose);

void explicitScheme_v2(int Nx, double tFinal, double tStep);

void implicitScheme(int n, int tFinal, double tStep, bool verbose);

void implicitScheme_v2(int n, double tFinal, double tStep);

void crankNicolsonScheme(int n, int tFinal, double tStep, bool verbose);

void crankNicolsonScheme_v2(int n, double tFinal, double tStep);

void diffusion1D();

void diffusion2D();

void diffusion2DLithosphere();

int JacobiSolver(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol);

int JacobiSolverLithosphere(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol, double Qt, double k, double beta);

double heatProduction(double time);

double Qdepth(int j, double dt);

void run_5c();

#endif