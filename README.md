# Direct ACE
----

Direct ACE is a C++ based algorithm to obtain numerically approximate
solutions of the direct Ising problem, that is, to compute the free energy and
the equilibrium observables of spin systems with arbitrary two-spin interactions.
It is described in the publication (available at https://link.springer.com/article/10.1140/epjb/e2019-100313-9)

Cocco S., Croce G., Zamponi F., Adaptive cluster expansion for Ising spin models, EPJB (2019)


The directACE code is based on the code used for the inverse Ising problem in https://github.com/johnbarton/ACE.

References for the inverse Ising problem:
- John P Barton, Eleonora De Leonardis, Alice Coucke, and Simona Cocco. Ace: adaptive cluster
expansion for maximum entropy graphical model inference. Bioinformatics, 32(20):3089–
3097, 2016.
- Simona Cocco and Rémi Monasson. Adaptive cluster expansion for the inverse ising problem:
convergence, algorithm and tests. Journal of Statistical Physics, 147(2):252–314, 2012.


*Please note that a real documentation is not available yet, but we report here some usage examples to get started.*


## Download and install

1) Download your local copy of the repository
2) Run the following commands in the terminal 
```bash
$ ./run_make.sh
```
2) this should compile the codes and produce two output files "directACE_spin_pm1.out" and "directACE_spin01.out" respectively for Ising and boolean spins variables.


# Input file
The input file should contain the hamiltonian of the system you want to study.

Here we focus on Ising spin variables, if you wish to use Boolean representation just replace everywhere "spin_pm1" with "spin01"

As an example, let's consider a system of N spins interacting via the
following hamiltonian:

H=- \sum_{i=1}^{N} J_{ij} s_i s_j - \sum h_i s_i 		 s_i={-1,1}

the input file (which we name "input.hj") must have the following structure:
```bash
 N (the number of spins)
 h_1
 h_2
 ...
 h_N
 J_{12}
 J_{13}
 ...
 J_{(N-1)N}
```
Pleas note the the input file should have extension *.hj*  

# Running the algorithm 
To run the directACE code just type:
```bash
./directACE_spin_pm1.out -d $PWD -i $INPUT_FILE  
```

where $PWD is the curent working directory directory, if you want to specify a different one just replace $PWD with your working directory.
$INPUT_FILE is the name of the input file *WITHOUT the extensions .hj*.

For example if my Hamiltonian is stored in *input_random_ER.hj*, so $INPUT_FILE should be *input_random_ER*

This will produce three files:

- out.sce: for the free energy
- out.mc: for the magnetizations
- out.ccs: for the connected correlations

optional flags are
```
-v: produce verbose output
-b: to specify the inverse temperature 
-kmax: to specify a max cluster length 
-kmin: to specify a minimal cluster length (the algorithm will construct all clusters up to kmin)
-tmin: to specify a minimal theta that when reached the algorithm will stop.
```
eg.
```bash
./directACE_spin_pm1.out -d . -i input_random_ER  -v -b 0.7 -kmin 4
```
