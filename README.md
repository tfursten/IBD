IBD
===
Isolation by distance simulation

A spatially explicit individual-based simulation to model dispersal on a lattice using different dispersal distribution functions.

Author
------
Tara Furstenau  
Biodesign Institute  
Center for Evolutionary Medicine and Informatics  
Arizona State University  
[Website](http://tfursten.github.io)  

Description
-----------
In the simulation, a population exists on a NxN rectangular lattice with either periodic boundaries (a torus) or absorbing boundaries. Individuals are uniformly distributed on the lattice with a single individual per node. Individuals are haploid and contain a single neutral genetic locus. (TODO: diploid option)

In the initial generation, the population contains the maximum number of individuals allowed by the landscape and each individual is assigned a unique allele (represented by an integer). During every discrete generation cycle, all individuals reproduce by producing a set number of clonal offspring.  These offspring experinece mutations according to the infinite alleles model at rate mu.  A burn-in period may be set to allow the population to reach a drift/mutation equilibrium. 

The offspring will next disperse from their orgininal cell according to a set dispersal distribution.  As offspring land on their destination cell they are immediately accepted or rejected using a reservoir sampling method, which allows the offspring to be uniformly sampled at each location as they arrive instead of storing them all in memory (Vitter, 1985).  When dispersal is complete, there is a maximum of one offspring per cell and that offspring becomes a parent in the next generation. 

There are currently 9 different dispersal distributions: exponential, gamma, half-normal, Pareto, Rayleigh, Rice, ring, and uniform. Some distributions have two implementations where one version is faster than the other.  The faster version is used by setting the --fast flag to true which is default. All uniform pseudo-random numbers are generated using an efficient xorshift algorithm (Marsaglia 2003).

###Exponential
The exponential dispersal function takes a single argument (sigma) and returns polar coordinates with exponetially distributed distances with rate 1/sigma and a uniform angle.  The distance values are drawn using an implementation of the ziggurat rejection sampling algorithm for the exponential distribution (Marsaglia and Tsang,2000).
###Gamma
The gamma dispersal function takes two arguments, sigma and alpha, and the beta parameter is calculated so that the second moment of the distribution is equal to 2*sigma^2. This function returns polar coordinates with gamma distributed distances and a uniform angle. The distance values are generated using Marsaglia's (2000) rejection sampling algorithm which takes advantage of the fast ziggurat procedure for generating normally distributed PRNs (Marsaglia and Tsang, 2000).    
###Half-Normal
The half-normal dispersal function takes a single argument, sigma, and returns polar coordinates with half-normal distributed distances with variance parameter sigma*sqrt(2) and a uniform angle. The distance values are the absolute value of draws from a normal distribution using an implementation of the ziggurat rejection sampling algorithm (Marsaglia and Tsang, 2000).
###Pareto
The Pareto dispersal function takes two arguments, sigma and alpha and the Xmin parameter is calculated so that the second moment of the distribution is equal to 2*sigma^2. The function returns polar coordinates with Pareto distributed distances and a uniform angle. The distance values are generated using inverse transform sampling, which, due to the relationship between the Pareto and exponential distributions, takes advantage of the fast ziggurat procedure for generating exponentially distributed PRNs.  
###Rayleigh
The Rayleigh dispersal function may take one or two arguments, sigma_x and sigma_y.  If only a single argument is given, sigma will be the same for each dimension and result in an isometric 2-dimensional normal distribution.  For the sake of uniformity, this was called the Rayleigh dispersal function because it originally returned polar coordinates with Rayleigh distributed distances and a uniform angle; however, this required an inefficient converion from polar to Cartesian coordinates. We instead draw axial offset distances from a normal distribution with variance sigma_x^2 (or sigma_y^2) using an implementation of the ziggurat rejection sampling algorithm.  There is a "fast" version of this disperal function available which precalculates the probability of dispersing on a discrete lattice.  Once the dispersal probabilities are calculated we use the efficient Alias method to sample from this discrete probability distribution.  
###Rice

###Ring
###Triangular
The Triangular dispersal function takes a si
###Uniform






Compiling from Source Code
--------------------------
IBD requires [CMake 2.8](http://www.cmake.org/) to build from source. 

1. Download the source code.  
2. Decompress the tar-bzip archive  
  ```
  tar xvzf IBD-*.tar.bz2
  ```
3. Change to the build directory.  
  ```
  cd IBD-*/build
  ```
4. Run the CMake build system.  
  ```
  cmake ..
  ```  
5. Compile  
  ```
  make
  ```

Dependencies
-------------
The Boost c++ Library is required for compilation and usage.
* Foreach  
* Program Options  

Run
----
Usage:
```
$ ./ibd config.txt
$ ./ibd --help
Allowed Options:

General Options:
  --help                Produce help message

Configuration:
  -x [ --maxX ] arg (=100)              Set X dimension
  -y [ --maxY ] arg (=100)              Set Y dimension
  -g [ --generations ] arg (=10)        Set number of Generations to run after 
                                        burn-in
  -o [ --offspring ] arg (=10)          Set number of offspring per individual
  -m [ --mut ] arg (=0)                 Set mutation rate
  -d [ --distribution ] arg (=triangular)
                                        Set Dispersal Distribution
  -s [ --sigma ] arg (=2)               Set dispersal parameter
  -b [ --burn ] arg (=0)                Set Burn-in Period
  -t [ --sample ] arg (=1)              Sample every n generations after 
                                        burn-in
  -f [ --output_file ] arg (=data)      Output File Name
  --seed arg (=0)                       Set PRNG seed, 0 to create random seed
  --landscape arg (=torus)              Set boundary conditions: torus or 
                                        rectangular
  --transect arg (=0)                   Set position of transect in X axis.
  --verbose arg (=0)                    Print data to screen
  --sparam arg (=0)                     Extra Parameter for dispersal
  --fast arg (=1)                       Use fast dispersal when available


