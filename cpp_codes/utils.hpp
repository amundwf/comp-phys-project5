#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <armadillo>


//double* Thomas_algo(int n, double x_0, double x_np1, int a, int b)

void writeGeneralMatrixToCSV_noLabels(arma::mat results, std::string filename, std::string directory);

void writeGeneralMatrixToCSV(arma::mat results, arma::field<std::string> columnLabels, std::string filename, std::string directory);

void explicitScheme(int n, int tFinal, bool verbose);

void implicitScheme(int n, int tFinal, double tStep, bool verbose);

void crankNicolsonScheme(int n, int tFinal, double tStep, bool verbose);

void analytical_solution_1D(int n_x, double tFinal, double tStep, int N_sum);

void diffusion1D();

void diffusion2D();

void diffusion2DLithosphere();

arma::vec ThomasAlgorithm(int n, arma::vec u, double a, double b, bool verbose);

int JacobiSolver(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol);

int JacobiSolverLithosphere(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol, double Qt, double k, double beta);

double heatProduction(double time);

double Qdepth(int j, double dt);

#endif