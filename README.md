# StateToState
Classical mechanics with initial and final positions specified.

The algorithm is described here: https://arxiv.org/abs/2106.15692

The code is written in C and compiles with gcc. It ouputs the coordinates as a function of time and writes them to a json file.  The C code requires the installation of the GNU scientific library (gsl). An example parameters file, example_params is also provided. The name of the executable produced is lmin. 

There are two python scripts that display the results. These require SciPy. 
