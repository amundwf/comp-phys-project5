#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <armadillo>


//double* Thomas_algo(int n, double x_0, double x_np1, int a, int b)

void writeGeneralMatrixToCSV_noLabels(arma::mat results, std::string filename, std::string directory);

void writeGeneralMatrixToCSV(arma::mat results, arma::field<std::string> columnLabels, std::string filename, std::string directory);

void explicitScheme(int n, int tFinal);

void implicitScheme(int n, int tFinal, double tStep);

void crankNicolsonScheme(int n, int tFinal, double tStep);

void diffusion1D();

void diffusion2D();

arma::vec ThomasAlgorithm(int n, arma::vec u, double a, double b, bool verbose);

#endif