# The diffusion equation: Analytical and numerical solutions in one and two spatial dimensions applied to studying the temperature distribution of the lithosphere

## Description

The github repository was made to solve the one dimensional and two dimensional diffusion equation. The methods for solving the one dimensional case are the: Explicit Scheme, Implicit Scheme and Crank Nicolson. The two dimensional case cannot be solved with these methods and so the Jacobi iterative solver is used. The code script is used to model the Earth's lithosphere (upper crust, lower crust and mantle) after the mantle is enriched with Uranium, Thorium and Potassium approximately one Giga years ago. The diffusion equation therefore models the temperature distribution in this problem. 

## Usage

The c++ scripts can be run from the make file within ./cpp_codes. There you will be asked the kind of problem you would like to run.


Below is an example...
```bash
Q1: Please enter if you want to run in 1D or 2D diffusion equation (int)
2

Q2: Do you want to run the Lithosphere problem [Y] or Unit Test [N]? [Y/N]
Y

Q3: Do you want to run before enrichment? [Y/N]?
Y

Q4: Starting from 1 Giga year ago, how long do you want to run for in Gy (double)
1.0

Q5: Please enter a time step (dt) in Gy (double)
0.01

Q6: Please enter max depth in km. Note that there are changes to heat production between 0 and 120 km
120

Q7: Please enter the distance step (dx) in km
0.01

Code running ...
```

The plotting scripts are written in python. They can be found in the folder ./py_codes. There are several options for one dimensional plotting or two dimensional plotting.
## Contributors

James Claxton and Amund Fredriksen
