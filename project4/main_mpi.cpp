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

    int numprocs, my_rank;
    int n, mc_cycles;
    double T_init, T_end, T_step;


    //MPI initialize
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    srand( 237*(my_rank + 13) );

    if ((my_rank == 0) && (argc <= 5)){
        cout << "ERROR: Bad usage! Please see the README.mp file for command line arguments.";
        exit(1);
    }

    // files to write to
    ofstream plotfile, acceptfile, probfile;

    if ((my_rank == 0) && ( atoi(argv[1]) >= 1 )){
        string out ("../files_python/mc_cycles");
        out += argv[1];
        out += argv[3];
        out += "mpi";
        out += ".txt";
        plotfile.open(out);
    }




    // convert input arguments only for master
    if ((my_rank == 0) && (argc > 1)) {
        n = atoi( argv[1] );
        mc_cycles = atoi( argv[2] );
        T_init = atof( argv[3] );
        T_end = atof( argv[4] );
        T_step = atof( argv[5] );
    }

    MPI_Bcast (&mc_cycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&T_init, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&T_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&T_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double k_B = 1.0;
    // counting accepted configurations
    int accepted_configs = 0;


    double J = 1.0;
    double w[17];
    double  TimeStart, TimeEnd, TotalTime;
    TimeStart = MPI_Wtime();
    // TEMP LOOP
    for (double t = T_init; t <=T_end; t+=T_step){
        double beta = 1.0 / (k_B*t);

        if(my_rank == 0) {
            cout << "Working with T = " << t << endl;
        }

        // initialize array
        for (int dE = -8; dE <= 8; dE++){
            w[dE + 8] = 0;
        }
        // every eight element is an exponential of energy(-8, -4, 0, 4, 8)
        for (int dE = -8; dE <= 8; dE+=4){
            w[dE + 8] = exp(- beta * dE);
        }

        mat L = zeros<mat>(n, n);

        // make each random number be random every time we call rand()
        // srand (time(NULL));

        // create a random matrix
        CreateRandomMatrix(n, L);

        // create an ordered matrix
        //CreateMatrix(n, L);

        // initialize energy and magnetization
        double E = CalculateEnergy(n, J, L);
        double M = CalculateMagneticMoment(L);


        for (int k = 0; k <= 10000; k ++){
            // metropolis algorithm
            Metropolis(n, E, M, w, L, accepted_configs);
        }

        // values to add to after each monte carlo cycle
        double energy_sum = E;
        double mag_sum = M;
        double mag_abs_sum = fabs(M);
        double energy_sqr = E * E;
        double mag_sqr = M * M;

        for (int k = 0; k <= mc_cycles; k ++){

            // metropolis algorithm
            Metropolis(n, E, M, w, L, accepted_configs);

            energy_sum += E;
            mag_sum += M;
            mag_abs_sum += fabs(M);
            energy_sqr += E * E;
            mag_sqr += M * M;

        }

        // mean energy
        double energy_sum_total = 0;
        double energy_sqr_sum_total = 0;
        double mag_sum_total = 0;
        double mag_abs_sum_total = 0;
        double mag_sqr_sum_total = 0;
        MPI_Reduce(&energy_sum, &energy_sum_total, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        MPI_Reduce(&energy_sqr, &energy_sqr_sum_total, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        MPI_Reduce(&mag_sum, &mag_sum_total, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        MPI_Reduce(&mag_abs_sum, &mag_abs_sum_total, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        MPI_Reduce(&mag_sqr, &mag_sqr_sum_total, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

        int mc_cycles_total = numprocs * mc_cycles;
        double meanE = energy_sum_total / mc_cycles_total;
        double meanESquared = energy_sqr_sum_total / mc_cycles_total;

        double meanAbsMag = mag_abs_sum_total / mc_cycles_total;
        double meanMagSquared = mag_sqr_sum_total / mc_cycles_total;
        double cv = (meanESquared - meanE*meanE) / (t*t);
        double chi = (meanMagSquared - meanAbsMag*meanAbsMag) / t;

//        double meanE = energy_sum / ((double) (mc_cycles));

//        // mean absolute magnetic moment
//        double meanM = mag_sum / ((double) (mc_cycles));

//        // specific heat
//        double mean_energy2 = energy_sqr / ((double) (mc_cycles));
//        double cv = ( mean_energy2 - meanE * meanE ) / (t*t);

//        // magnetic susceptibility
//        double mean_mag2 = mag_sqr / ((double) (mc_cycles));
//        double chi = ( mean_mag2 - meanM * meanM ) / t;


//        double meanE_tot = 0;
//        MPI_Reduce(&meanE, &meanE_tot, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
//        double meanM_tot = 0;
//        MPI_Reduce(&meanM, &meanM_tot, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
//        double cv_tot = 0;
//        MPI_Reduce(&cv, &cv_tot, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
//        double chi_tot = 0;
//        MPI_Reduce(&chi, &chi_tot, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);


        if(my_rank == 0){
            //plotfile << meanE_tot << "  " << meanM_tot << "  " << cv_tot << "  " << chi << "  " << t << endl;
            plotfile << meanE/mc_cycles << "  " << meanAbsMag/mc_cycles << "  " << cv/mc_cycles << "  " << chi/mc_cycles << "  " << t << endl;
        }
    }


    TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd-TimeStart;
    if ( my_rank == 0) {
        cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
    }

    MPI_Finalize();
    return 0;
}


