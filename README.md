To compile just run "make" after cloning. Several executables will be generated:

`./estimation <input filename> <output filename>`

It runs a single estimation based on given model parameters in the input file.
This executable could be used by an external program which tries to find optimial values for the parameters of the model.
The program takes a file cobtaining the model's parameter and output a single value which represent how well these parameters fit the model inthe file provided as the output file.
The framework that was used for the actual estimation is appspack: https://software.sandia.gov/appspack/version3/index.html

However, any optimization framework that has similar file interface would work.

`./estimation_test <input filename> <output filename>`

It runs a single estimation based on given model parameters in the input file, similarly to "estimation" however, it also output to screen many statistics on the model, that could eb used for manual estimation of the fit of the model.

`./estimation_sim <input filename> <output filename> <simulation type 0-3,5,6> [percent 0-100]`

Run policy simulation based on a set of inputs, as well as the type and parameters for the simulation. These are the possible policies that could be simulated:
* 0 - None (similar to estimation_test)
* 1 - Rent
* 2 - Wage
* 3 - Travel Cost
* 5 - Future Interest
* 6 - Married Only

Note that simulation of type 4 (lump sum) has its own executable: estimation_sim4

Other executables are used as utility programs mainly to log all internal state transitions of the model, for debugging purposes.
