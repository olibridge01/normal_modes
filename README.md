# Part II Chemistry Programming Course
Oli Bridge, <ob344@cam.ac.uk>

St Catharine's College, Cambridge

## Exercise 2: Data processing and plotting

In this exercise, a program was written that calculates the equilibrium geometry and normal modes for a triatomic molecule (H2O or H2S). It also plots the potential energy surface as a function of bond length and bond angle. The program takes as an input the name of the directory containing Gaussian output files for the triatomic at various geometries.

The program will display a PE surface, and a quadratic fit plot (for determining the normal modes) for small oscillations in r and theta. These figures are also saved as .pdf files to the directory containing the python file.

The program requires the output file directories (H2Ooutfiles and H2Soutfiles) from Moodle as an input (these can be found in the Programming Practical section on Moodle).

## Libraries Required
The program runs on Python3 with the following libraries:
- numpy
- matplotlib
- sys
- os

## User Input
To run the program, execute the following command in the directory containing the python file:
```
python3 ex2.py /[Gaussian output file directory]/
```
NB. Make sure the Gaussian file directory is in the same folder as the python file.

## Test Example
```
python3 ex2.py /H2Ooutfiles/

Equilibrium geometry for H2O:
r / Angstroms        Theta / Degrees      Energy / Hartrees   
-------------------------------------------------------------
    0.95                  105.0            -76.0242882171 
-------------------------------------------------------------

Symmetric stretch: v1 = 4043.9 cm-1
Symmetric bend: v2 = 1711.6 cm-1

```

## Notes 
The program uses the matplotlib library to plot PE surfaces for the triatomics as well as plot the quadratic fit functions that are used to determine the normal modes.
