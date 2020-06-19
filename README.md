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
./bn_smc NODES OBS ITER CHAINS CUTOFF CORES TEMPERATURE DEPTH DATA_PATH PARAM_PATH OUTPUT_NAME PRIOR_PATH
```

Input parameters:
* `NODES`: number of nodes in the data set (a.k.a number of variables)
* `OBS`: number of observations (data size)
* `ITER`: maximum iterations will be run for each SMC chain (depends on NODES, typically 1000)
* `CHAINS`: number of SMC chains to make, the more the better, but take more time to finish
* `CUTOFF`: The cut off point where MIT test determine dependency, typically 0.05, you can set to 0.01 if needed
* `CORES`: Number of cores (threads, CPUs) to use, it gets replaced by number of cores available if the later one is smaller.
* `TEMPERATURE`: used to control the greedy level, 
* `DEPTH`: How many rounds of 3rd stage to perform, dont make it too big, typically 3-5
* `DATA_PATH`: the path to the data file (please include the full file name)
* `PARAM_PATH`: the path to the param file (please include the full file name)
* `OUTPUT_NAME`: the prefix of output files
* `PRIOR_PATH`: the path to the prior information file (please include the full file name)

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
This file contains 2 columns, each row sepcifies the name of the node, and how many categories are there in the corresponding variable (node).

Sample file:

```
Node1,2
Node2,3
Node3,3
Node4,2
Node5,3
Node6,2
Node7,3
Node8,2
Node9,3
Node10,3
Node11,2
Node12,3
Node13,2
Node14,2
Node15,3
Node16,4
Node17,2
Node18,4
Node19,2
Node20,3
Node21,3
Node22,3
Node23,2
Node24,2
Node25,3
Node26,4
Node27,2
Node28,3
Node29,4
Node30,4
Node31,4
Node32,4
Node33,3
Node34,2
Node35,3
Node36,3
Node37,3
```

# Prior file
this is a csv file with four columns in each row. Each row represents a prior information: from_node, to_node, p_val, probability. The nodes are represented by integers starting from 1. 

Sample file:

```
Node1, Node3, 0.05, 0.95
Node2, Node4, 0.1, 0.99
```

# Program Outputs

The program will generate a directory `OUTPUT_NAME` under the directory of the executable. It contains the following files or directories:
* `logs` directory: this directory contains the step-by-step debug information for the 2nd step of the algorithm for each monte carlo chain. The file name indicates the random seed used to generate the MC chain.
* `MIT.csv` and `CMIT.csv` are the debug information for the test results of MIT test and CMIT test in the first step.
* `after*.csv` are the debug information for the networks in first step after MIT test/CMIT test/symmetry correction. `afterSymmetryCorrection.csv` is the final network after the first step.
* `summaries` is the directory that contains the final results:
    * `best_BN.txt` contains the model string of the best bayesian network. The model string can be imported to BNlearn package via the [model2network](https://www.bnlearn.com/documentation/man/modelstring.html) function.
    * `bn.tsv` is the tsv file that summarizes all the MC chains during the run. For each chain, it records: 
        1. is this chain the best MC chain of the run? (or being score equivalent to the best one)
        2. the index
        3. random seed used for this chain
        4. prior log-likelihood
        5. log-likelihood from data (BIC)
        6. posterior log-likelihood
        7. final model string
    * other files are as-is since I inherited this program.