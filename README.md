BN_SMC
======

Sequential monte carlo method apply on structure learning of a Bayesian Network.

# Dependencies:

# Compile
Make sure you have the following dependencies installed:
* CMake >= 3.8
* Boost>=1.53
* an OpenMP library

To make the `bn_smc` binary, use cmake under project root directory to generate `Makefile`, then use `Makefile` to compile the binary. 

An example:

```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make 
```

Then you will have `bn_smc` binary under `./build` directory.

If your boost is installed in a custom location, use flag `-DBOOST_PATH=<location>` when calling cmake to help it find boost library.

# Run the program
```
./bn_smc NODES OBS iter chains cutoff cores temperature depth DATA_path PARAM_path output_path
```

Input parameters:
* `NODES`: number of nodes in the data set (a.k.a number of variables)
* `OBS`: number of observations (data size)
* `iter`: maximum iterations will be run for each SMC chain (depends on NODES, typically 1000)
* `chains`: number of SMC chains to make, the more the better, but take more time to finish
* `cutoff`: The cut off point where MIT test determine dependency, typically 0.05, you can set to 0.01 if needed
* `cores`: Number of cores (threads, CPUs) to use, it gets replaced by number of cores available if the later one is smaller.
* `temperature`: used to control the greedy level, 
* `depth`: How many rounds of 3rd stage to perform, dont make it too big, typically 3-5
* `DATA_path`: the path to the data file (please include the full file name)
* `PARAM_path`: the path to the param file (please include the full file name)
* `output_path`: the prefix of output files
* `PRIOR_path`: the path to the prior information file (please include the full file name)

# Data
The data file is a csv file, each row is one observation, while each column is a variable (node).

Sample file (`nodes`=37, `obs`=5):
```
2,2,2,2,2,2,2,2,3,3,2,3,2,2,1,2,2,4,2,3,3,2,2,1,1,4,1,2,3,1,2,4,1,2,3,3,1
2,1,1,2,1,2,2,2,3,2,1,2,2,2,2,2,2,4,2,3,3,3,2,1,1,4,1,2,3,1,2,4,1,2,3,3,3
2,2,2,2,2,2,2,2,3,3,2,3,2,2,3,2,2,1,2,1,1,2,2,1,1,4,2,2,3,2,1,1,3,2,3,3,3
2,2,2,2,2,2,2,2,3,3,2,3,2,2,1,2,2,1,2,1,1,2,2,1,1,2,2,2,3,2,1,1,3,2,3,1,1
2,2,2,2,2,2,2,2,1,1,2,1,2,2,3,2,2,4,2,3,3,2,2,1,1,2,2,3,4,4,2,4,1,1,2,1,2
```

# Param file
This file contains only a column, each row sepcifies how many categories are there in the corresponding variable (node).

Sample file:

```
2
3
3
2
3
2
3
2
3
3
2
3
2
2
3
4
2
4
2
3
3
3
2
2
3
4
2
3
4
4
4
4
3
2
3
3
3
```

# Prior file
this is a csv file with four columns in each row. Each row represents a prior information: from_node, to_node, p_val, probability. The nodes are represented by integers starting from 1. 

Sample file:

```
1, 3, 0.05, 0.95
2, 4, 0.1, 0.99
```