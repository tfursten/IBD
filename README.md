IBD
===
Isolation by distance simulation

A spatially explicit individual-based simulation which will model dispersal on a toroidal lattice using different dispersal distribution functions.


Author
------
Tara Furstenau  
Biodesign Institute  
Center for Evolutionary Medicine and Informatics  
Arizona State University  
[Website](http://tfursten.github.io)  

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

Run
----
Usage:
```
./ibd config.txt
./ibd --help

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
  -d [ --distribution ] arg (=exponential)
                                        Set Dispersal Distribution
  -s [ --sigma ] arg (=2)               Set dispersal parameter
  -b [ --burn ] arg (=0)                Set Burn-in Period
  -t [ --sample ] arg (=1)              Sample every n generations after 
                                        burn-in
  -f [ --output_file ] arg (=data)      Output File Name
  --seed arg (=0)                       Set PRNG seed, 0 to create random seed
  --transect arg (=0)                   Set position of transect in X axis.
```
