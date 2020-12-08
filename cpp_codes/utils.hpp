#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <armadillo>


void writeGeneralMatrixToCSV_noLabels(arma::mat results, std::string filename, std::string directory);

void writeGeneralMatrixToCSV(arma::mat results, arma::field<std::string> columnLabels, std::string filename, std::string directory);

// Functions for the 1D diffusion equation.

double stabilityConditionExplicit_dt(double dx);

void explicitScheme(int n, int tFinal, bool verbose);

void implicitScheme(int n, int tFinal, double tStep, bool verbose);

void crankNicolsonScheme(int n, int tFinal, double tStep, bool verbose);

void analytical_solution_1D(int n_x, double tFinal, double tStep, int N_sum);

void diffusion1D();

arma::vec ThomasAlgorithm(int n, arma::vec u, double a, double b, bool verbose);

// Solver for the 2D diffusion equation.

void diffusion2D();

int JacobiSolver(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol);

// Functions for after radioactive enrichment

void diffusion2DAfterEnrichment();

int JacobiSolverAfterEnrichment(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol, double Qt, double k, double beta);

double Qtime(double time);

double Qdepth(int j, double dt);

// Temporary functions for before enrichment..

void diffusion2DBeforeEnrichment();

int JacobiSolverBeforeEnrichment(int N, double dx, double dt, arma::mat &A, arma::mat &A_prev, double abstol, double eta);

#endif