/*
 argv[1] : number of spins L
 argv[2] : number of Monte Carlo cycles
 argv[3] : initial temperature
*/

#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <fstream>
#include "functions.h"
#include "mpi.h"

using namespace std;
using namespace arma;


int main(int argc, char* argv[]){

    if (argc <= 3){
        cout << "ERROR: Bad usage! Please see the README.mp file for command line arguments.";
        exit(1);
    }

    // files to write to
    ofstream outputfile, plotfile, acceptfile, probfile;

    if ( atoi(argv[1]) >= 1 ){
        string out ("../files_python/mc_cycles");
        out += argv[1];
        out += argv[3];
        out += ".txt";
        plotfile.open(out);
    }
    if ( atoi(argv[3]) >= 1 ){
        string out ("../files_python/accepted");
        string prob ("../files_python/probability");
        out += argv[3];
        //out += "ran";  // for the random matrix
        prob += argv[3];
        out += ".txt";
        prob += ".txt";
        acceptfile.open(out);
        probfile.open(prob);
    }


    // convert input arguments
    int n = atoi( argv[1] );
    //int s = n*n;
    int mc_cycles = atoi( argv[2] );
    double T_init = atof( argv[3] );
    //change for e
    double T_end = atof( argv[4] );
    double T_step = atof( argv[5] );

    double k_B = 1.0;
    double beta = 1.0 / ( k_B * T_init );
    double J = 1.0;
    double w[17];

    // initialize array
    for (int dE = -8; dE <= 8; dE++){ w[dE + 8] = 0; }
    // every eight element is an exponential of energy(-8, -4, 0, 4, 8)
    for (int dE = -8; dE <= 8; dE+=4){ w[dE + 8] = exp(- beta * dE); }

    // counting accepted configurations
    int accepted_configs = 0;

    mat L = zeros<mat>(n, n);

    // make each random number be random every time we call rand()
    srand (time(NULL));

    // create a random matrix
    CreateRandomMatrix(n, L);

    // create an ordered matrix
    //CreateMatrix(n, L);

    // iitialize energy and magnetization
    double E = CalculateEnergy(n, J, L);
    double M = CalculateMagneticMoment(L);

    // values to add to after each monte carlo cycle
    double energy_sum = E;
    double mag_sum = M;
    double energy_sqr = E * E;
    double mag_sqr = M * M;

    cout << E << endl;
    // full monte carlo cycles -> for exercise b and c

    for (int k = 1; k <= mc_cycles; k ++){

        // metropolis algorithm
        Metropolis(n, E, M, w, L, accepted_configs);

        energy_sum += E;
        mag_sum += fabs(M);
        energy_sqr += E * E;
        mag_sqr += M * M;

        // write values to file
        plotfile << E/k << "  " << fabs(M)/k << "  " << k << endl;
        acceptfile << accepted_configs/k << "  " << k << endl;

    }

    cout << E << endl;

    /*
    // we calculate the energy etc. until we know that we are in a steady state -> for exercise d
    for (int k = 1; k <= 1000; k ++){

        // metropolis algorithm
        Metropolis(n, E, M, w, L, accepted_configs);

    }

    // we begin the metropolis sampling
    for (int k = 1001; k <= mc_cycles; k ++){

        // metropolis algorithm
        Metropolis(n, E, M, w, L, accepted_configs);

        probfile << E << " " << k << endl;

        energy_sum += E;
        energy_sqr += E * E;
    }

    */

    // write values to file for the 2x2 lattice -- specialized for exercise b
    if ( atoi(argv[1]) == 2 && atoi(argv[2]) >= 1000000 ){

        // mean energy
        double mean_energy = energy_sum / ((double) (mc_cycles+1));

        // mean absolute magnetic moment
        double mean_mag = mag_sum / ((double) (mc_cycles+1));

        // specific heat
        double mean_energy2 = energy_sqr / ((double) (mc_cycles+1));
        double spec_heat = ( mean_energy2 - mean_energy * mean_energy ) * beta / T_init;

        // magnetic susceptibility
        double mean_mag2 = mag_sqr / ((double) (mc_cycles+1));
        double mag_susc = ( mean_mag2 - mean_mag * mean_mag ) * beta;

        outputfile.open("../files_python/2x2_output.txt");
        outputfile << "Mean energy: " << mean_energy  << endl
                   << "Mean magnetic moment: " << mean_mag<< endl
                   << "Mean energy squared: " << mean_energy2 << endl
                   << "Mean energy x mean energy: " << mean_energy* mean_energy << endl
                   << "Mean magnetic moment squared: " << mean_mag2 << endl
                   << "Mean magnetic moment x mean magnetic moment: " << mean_mag * mean_mag << endl
                   << "Specific heat: " << spec_heat << endl
                   << "Magnetic susceptibility: " << mag_susc << endl;
        outputfile.close();
    }




    cout << "End" << endl;



    return 0;
}
