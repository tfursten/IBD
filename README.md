IBD
===
Isolation by distance simulation

A spatially explicit individual-based simulation which will model dispersal on a toroidal lattice using different dispersal distribution functions.


Author
------
Tara Furstenau  
Biodesign Institute  
Center for Evolutionary Medicine and Informatics  
[Website](http://tfursten@github.io)  

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
5.Compile
```
make
```
6. Install
```
make install
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
```
