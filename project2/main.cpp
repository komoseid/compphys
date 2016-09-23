// 23.09.2016
//Morten syns vi skal lage vår egen versjon av jacobi.cpp og teste med matrisen vi lager i dette programmet.
#include <iostream>
#include <jacobi.h>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;


void printMatrix(double** A, int n) {
    for (int i=0; i < n; i++) {
        for (int j=0; j < n; j++) {
            cout << A[i][j] << "  ";
        }
        cout << endl;
    }
}

int main() {
    int N = 4;
    double **A = new double*[N];
    double **R = new double*[N];

    //declare variables
    double *rho, h, *V, e, *d;
    rho = new double[N];
    V = new double[N];
    d = new double[N];

    //Boundary conditions
    rho[0] = 0;
    rho[N] = 5;

    //Calculate step length
    h = (rho[0] - rho[N])/N;

    //Fill in rho
    for(int i=1; i<N; i++){
        rho[i] = rho[0] + i*h;
    }

    //fill HO potential V and the diagonal elements d
    for(int i=0;i<N; i++){
        V[i]=pow(rho[i],2);
        d[i]= (2./pow(h,2)) + V[i];
    }

    //non-diagonal matrix element
    e = -(1./pow(h,2));

    //create matrix and fill matrix A
    for(int i=0; i<N; i++) {
        A[i] = new double[N];
        R[i] = new double[N];
        for(int j=0; j<N; j++) {
            if(i == j) {
                A[i][j] = d[j];

            }
            else if (i == (j-1) || j == (i-1)){
                A[i][j] = e;
            }
        }
    }
    //print matrix
    printMatrix(A, N);

    //j = jacobi_method(A,R,10);
    //cout << jacobi_method(A,R,10) << endl;
    return 0;
}
