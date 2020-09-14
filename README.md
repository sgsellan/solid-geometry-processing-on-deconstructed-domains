# Solid Geometry Processing on Deconstructed Domains
Public code release for the paper "[Solid Geometry Processing on Deconstructed Domains](http://dgp.toronto.edu/~sgsellan/pdf/overlapping.pdf)", presented at SGP 2019 and authored by Silvia Sellán, Herng Yi Cheng, Yuming Ma, Mitchell Dembowski and Alec Jacobson. Released under MIT licence.

## Installation
To install this library, please start by cloning the repository recursively
```
git clone --recursive https://github.com/sgsellan/solid-geometry-processing-on-deconstructed-domains.git
```
After this, we will build the mex functions in the `gptoolbox` directory:
```
cd solid-geometry-processing-on-deconstructed-domains/gptoolbox/mex
mkdir build
cd build
cmake ..
make
```

## Use
To replicate the results in the paper, start by adding `include-matlab` and `gptoolbox` to your Matlab path, for example using `addpath(genpath('include-matlab'));addpath(genpath('gptoolbox''));`. Once this is done, you can replicate each result in our paper by running each of the scripts in the `scripts` directory. The first commented line in each of them specifies by number which figure of the paper it corresponds to. The convergence tests in Figures 8 and 10 can be replicated by running the scripts in `scrips/convergence-tests`.

If you want to use our method on your own example meshes and equations for comparisons or applications, you can take a look at our  `example.m` file and substitute with your own code.


## Known Issues
Please do not hesitate to contact
[sgsellan@cs.toronto.edu](mailto:sgsellan@cs.toronto.edu) if you find any issues
or bugs in this code, or you struggle to run it in any way.
