#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <armadillo>


//double* Thomas_algo(int n, double x_0, double x_np1, int a, int b)

void writeGeneralMatrixToCSV_noLabels(arma::mat results, std::string filename, std::string directory);

void writeGeneralMatrixToCSV(arma::mat results, arma::field<std::string> columnLabels, std::string filename, std::string directory);

void explicitScheme(int n);

void task_5c();

#endif