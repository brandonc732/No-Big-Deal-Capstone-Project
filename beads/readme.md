Author: Brandon

# Overview

This folder contains all work done with the beads.exe program described in [![DOI](https://img.shields.io/badge/DOI-10.1214%2F09--AOAS299-blue)](https://doi.org/10.1214/09-AOAS299)

The program is a 2010 32bit c++ application retrieved from John Hughes. The original zip folder can be downloaded from this [link](https://drive.google.com/file/d/14dYG__a6HE4OaJqqPYbris_M_763Y-hw/view?usp=sharing)

I've made multiple changes to the folder and my c++ packages to allow the program to compile and run on Windows 10. I also added print outs throughout the program's main function to give an idea of where it is within it the process when running.

Overall, the program takes too long to run and doesn't produce useable results for our application. For NIST test data images, I had to remove extreme single-pixel peaks from the image by averaging them with the neighbors. This brought the runtime down from multiple hours to around one hour. And even with this, the program still chooses an optimal of 1 particle, so I created a version (small change) that continues the estimation until a given mininum number of particles. This of course defeats the purpose of the original algorithm and was only used to see if it would return reasonable position estimations at higher particle counts. Out of curiosity, it did return what seemed to be valid position estimations, but it's still not useable for the project. 

**note:** This program is single-threaded, which is partially why it's very slow

# Usage

### Original version

```bash
beads.exe width file [dumpfile]
```
- width : pixel width in nanometers
- file : path to tiff image
- dumpfile : optional file location to dump grid search information


### Minimum beads version

```bash
beads_min.exe width file min_beads [dumpfile]
```
- width : pixel width in nanometers
- file : path to tiff image
- min_beads : minimum number of beads for estimation
- dumpfile : optional file location to dump grid search information 




# Compiling

Since the software was so old, we had to add multiple c++ libraries to the compilation command to correctly compile. Since I don't program in c++, I opted to use vcpkg as an easy package manager for this. This gives me the compile command: 

```bash
g++ beads.cpp -o beads.exe -I"C:\Users\user\vcpkg\installed\x64-windows\include" -L"C:\Users\user\vcpkg\installed\x64-windows\lib" -ltiff -lgsl -lgslcblas
```

my vcpkg folder includes:
- gsl.lib
- gslcblas.lib
- jpeg.lib
- lzma.lib
- tiff.lib
- turbojpeg.lib
- zlib.lib










