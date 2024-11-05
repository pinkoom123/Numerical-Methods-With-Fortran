# Numerical-Methods-With-Fortran
Using Fortran to solve Numerical problems 

#Monte Carlo Simulations in Fortran

This repository contains a collection of Fortran routines that use Monte Carlo methods to perform numerical      estimations. It includes subroutines for estimating the value of π and evaluating definite integrals, leveraging both the Hit-or-Miss and Crude Monte Carlo methods. The code also demonstrates the Rejection method for generating points from a specific distribution function.

   #Key Features:

    #Pi Estimation:
        Hit-or-Miss Method: Estimates π by calculating the proportion of random points falling within a unit circle.
        Crude Monte Carlo Method: Provides an alternative approach to estimate π by leveraging integration.

    #Integral Estimation:
        Hit-or-Miss Method: Evaluates the integral of 3x23x2 over a defined range.
        Crude Monte Carlo: Calculates the integral of e−x2e−x2 over a defined range by averaging function values at random points.

    Rejection Method for Distribution Sampling:
        Generates points following the distribution f(x)=3x2f(x)=3x2 in the interval [0,1][0,1] using the Rejection #Sampling technique.
        Outputs sampled points to a file, distribution_points.txt, for further analysis or plotting.

#Files:

    monte_carlo.f90: Main module with subroutines for Monte Carlo estimations and sampling methods.
    distribution_points.txt: Output file storing sampled points generated using the Rejection method for the distribution f(x)=3x2f(x)=3x2.

#Instructions for Use:

    Compile and run the program in a Fortran-compatible environment.
    Adjust the value of N in the main program to control the number of iterations or samples for each estimation.

#Sample Output:

The program prints estimated values of π and integrals to the console, while distribution_points.txt saves points for the specified distribution.
