PROJECT 4 FYS4150

----------------------
Command Line Arguments
----------------------
For non-MPI case: <br />
arg1 = Lattice dimension <br />
arg2 = nr of MC cycles <br />
arg3 = Initial temperature <br />
Example:  20 1000000 2.0 <br />

For MPI case: <br />
Start commandline with "-n 4 ./project4" to specify nr of nodes to use and initialize MPI. <br />
After that: <br />
arg1 = lattice dimension <br />
arg2 = nr om MC cycles <br />
arg3 = initial temperature <br />
arg4 = end temperature <br />
arg5 = step length for temperature <br />
Example: -n 4 ./project4 40 100000 2.0 3.0 0.05 <br />
