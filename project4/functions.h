#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void CreateRandomMatrix(int n, mat &L) ;
void CreateMatrix(int n, mat &L) ;
double CalculateEnergy(int n, double J, mat &L) ;
double CalculateMagneticMoment(mat &L) ;
void Metropolis(int n, double &E, double &M, double *w, mat &L, int &accepted_configs) ;

#endif // FUNCTIONS_H
