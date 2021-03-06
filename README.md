#PSR_MC
Required sw: astropy, numpy, tempo2

Update (14/06/2017): Ran code overnight, all the chi_squareds its outputting are way too large (500-1000).  However, after some troubleshooting it seems like the problem is either (1) in the values we are gridding over (I am especially unsure of appropriate values for px and h3) OR in the fact that tempo2 is fitting only for the parameters we calculate in the program (PBDOT, XDOT, OMDOT, M2, SINI), since it seems to be reading the updated .par files perfectly fine.

NOTE: The manualcon_tempo2mc.py file does a manual conversion for the units (there may be some issues)... while tempo2_mc.py uses astropy.units to convert everything possible.

For PSR J1640+2224 there are complementary python and bash programs that work together to perform a Monte Carlo simulation over 4 variables (px, h3, omega, cosi) to try to find a better timing solution that gives a smaller error on the value of the companion's mass.

NOTE: Changing the number of parameters (of h3, parallax, omega, cosi) that the program will grid over requires changing both programs in superficial ways right now.

The bulk of the program occurs in the python file tempo2_mc.py.  Input the names of the .par and .tim files you will use with tempo2.  The .par file should not contain the parameters that are impacted by your choice of h3, px, cosi, or omega, since these values will need to be updated in a temporary .par file before calling tempo2.  It then defines various constants and several functions to relate between h3, px, cosi, omega and the parameters tempo2 needs for its T2 binary model, and creates arrays of values to grid over.  All of these can and should be updated as better measurements become available.

After this, the for loops begin.  The program iterates over each of the variables we want to grid over, calculates T2 parameters from these values, writes them to a temporary .par file, temp.par, then calls a bash script to run tempo2.  The bash script calls tempo2 with temp.par and the .tim file you specified at the start of the python program.  It then records the chi_squared/dof for the fit, and writes the variables in the MC as well as this chi_squared to an output text file, mc_output.txt, and removes the temporary par file.
