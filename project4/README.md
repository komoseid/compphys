PROJECT 4 FYS4150

----------------------
Command Line Arguments
----------------------
For non-MPI case:
arg1 = Lattice dimension
arg2 = nr of MC cycles
arg3 = Initial temperature
Example:  20 1000000 2.0

For MPI case:
Start commandline with "-n 4 ./project4" to specify nr of nodes to use and initialize MPI.
after that
arg1 = lattice dimension
arg2 = nr om MC cycles
arg3 = initial temperature
arg4 = end temperature
arg5 = step length for temperature
Example: -n 4 ./project4 40 100000 2.0 3.0 0.05
