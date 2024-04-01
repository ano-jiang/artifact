# Artifact
This library presents a parallel method to accelerate the generation of correctly rounded elementary mathematical functions.
The artifact includes: (1) 6 correctly rounded elementary functions, (2) correctness testing framework for the 6 functions, 
(3) performance testing framework to demonstrate the performance improvements over RLIBM and glibc, and 
(4) the parallel implementation of polynomial generator.
## Requirements
To replicate our experiments, we need a Linux machine with gcc compiler on x86-64 machine. To run the polynomial generator, 
it is recommended to have a machine with at least 16 GB of RAM. We need about 50 GB of storage if all the auxiliary files are unpacked.

## Installation
The installation of PRLIBM is the same as RLIBM-CGO(https://github.com/rutgers-apl/cgo23-artifact). There's one thing to note, 
in step three, the interval file of log2 function that published in RLIBM-CGO(https://github.com/rutgers-apl/cgo23-artifact)
is **NOT correct**( it should be the wrong version that was uploaded ), and you must recreate it using RLIBM-ALL(https://github.com/rutgers-apl/rlibm-all).


## Reproducing Results
###  Polynomial Generations of PRLIBM
We illustrate polynomial generation for the log2 function.
To use the polynomial generator, it needs reduced intervals generated. 
It also requires Soplex installed with the the SOPLEX_INCLUDE and SOPLEX_LIB environment variables set as described in RLIBM-
ALL(https://github.com/rutgers-apl/rlibm-all). After SOPLEX has been installed, `cd` to the root directory of 
the artifact, set the correct path in the file `env.sh`.
Next, to generate the polynomial for Log2, execute the following commands.
```shell
cd artifact/
. ./build.sh
cd polynomial_generator
mpirun.mpich -n 1 ../cmake-build-release-intel/polynomial_generator/polygen-estrin-fma ./cfg_files/log2-estrin-fma.txt <INTERVALS>/Log2Intervals
```
The configurations for Log2 polynomials that we generate are in log2-estrin-fma.txt file. At the end, the 
polynomial generator prints out the polynomial and detail execution time.

###  Polynomial Generations of PRLIBM-naive
To compare the performance of PRLIBM, we have also implemented the naive parallelization of RLIBM, called PRLIBM-naive.
To generate the polynomial for Log2, execute the following commands.
```shell
cd artifact/
. ./build.sh
cd prlibm-naive/polynomial_generator
mpirun.mpich -n 1 ../../cmake-build-release-intel/prlibm-naive/polynomial_generator/polygen-estrin-fma ./cfg_files/log2-estrin-fma.txt <INTERVALS>/Log2Intervals
```
Again, the polynomial generator prints out the polynomial and detail execution time.
### Correctness test
You can test the correctness of all the functions generated from our libraries using the correctness test infrastructure as follows. Firstly,
`cd` to the root directory of the artifact, set the correct path in the file `env.sh`. Next, run 
```shell
. ./build.sh
```
To test out the PRLIBM function, you can execute the following command, which checks if the implementation 
produces correctly rounded results for all inputs.
```shell
cd artifact/
./cmake-build-release-intel/correctness_test/correct_test_Log2 <ORACLE>/Log2Oracle
```
You should see an output like the following:
```
Wrong results: 000 (0)
prlibm-latest   wrong result: 0
```

### Performance
We provide an automated script to test the performance of the PRLIBM. To run the performance testing framework, 
execute the following command
```shell
cd artifact/performance_test
. ./runPRLIBM.sh
```
It automatically executes the implementations of the 6 functions and creates text files with the timing data in `logs` directory. 


You can also test the performance speedup of our functions over RLIBM and glibc with the following command. It is
advised to not execute other programs simultaneously with the script.
```shell
. ./runRLIBM.sh
. ./runGLIBC.sh
```
