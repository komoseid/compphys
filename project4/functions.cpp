#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdlib.h>     /* srand, rand */
#include "functions.h"

using namespace std;
using namespace arma;

void CreateRandomMatrix(int n, mat &L){

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){

            int r = rand() % 2;

            if (r == 1){
                L(i, j) = 1.0;
            }
            else {
                L(i, j) = -1.0;
            }
        }
    }

    return;
}

void CreateMatrix(int n, mat &L){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            L(i, j) = 1.0;
        }
    }
}

// function to calculate periodic boundary conditions and energy
double CalculateEnergy(int n, double J, mat &L){

    double E = 0;

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            int ip = (i + 1) % n;
            int jp = (j + 1) % n;
            E += -J*L(i,j)*( L(ip,j) + L(i,jp));
        }
    }

    return E;

    /*
    // vector to keep values calculated with periodic boundary conditions
    vec pbc = zeros<vec>(n*n);
    vec pbc2 = zeros<vec>(n*n);

    // calculating periodic boundary conditions
    for(int i = 0; i < n*n; i++){

        if( i == ( (i/n + 1) * n - 1 ) ){
            pbc(i) = L(i) * L(i-(n-1));
        }
        else {
            pbc(i) = L(i) * L(i+1);
        }

        if ( n + i > (n * n - 1) ){
            pbc2(i) = L(i) * L(i - n*(n-1));
        }
        else {
            pbc2(i) = L(i) * L(n+i);
        }
    }

    // energy of the Ising model
    double energy = - J * accu(pbc + pbc2);

    return energy;
    */

}

// function to calculate magnetic moment
double CalculateMagneticMoment(mat &L){

    // magnetization
    double mag = accu(L);

    return mag;
}


void Metropolis(int n, double &E, double &M, double *w, mat &L, int &accepted_configs){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){

//            int ix = rand() % n;
//            int iy = rand() % n;

            int RandomNumber = rand() % (n*n);
            int ix = RandomNumber / n;
            int iy = RandomNumber % n;

            // periodic boundary conditions
            int ip = (ix + 1) % n;
            int in = (ix - 1 + n) % n;
            int jp = (iy + 1) % n;
            int jn = (iy - 1 + n) % n;

            int delta_E = 2.0 * L(ix, iy) * ( L(ip, iy) + L(in, iy) + L(ix, jp) + L(ix, jn) );

            double r = rand() / ((double) numeric_limits<int>::max());
            if (r <= w[delta_E+8]){

                L(ix, iy) *= -1.0;
                E += (double) delta_E;
                M += (double) 2.0 * L(ix, iy);

                accepted_configs += 1;
            }
        }
    }

    return;
}








