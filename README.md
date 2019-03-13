# STEAM
Table-based device modeling with polynomial interpolants. This code implements the 2017 ASPDAC paper titled "STEAM: Spline-based tables for efficient and accurate device modelling". Several other features were added to the code after the paper was published and are described in greater detail in a technical report "Table-based Device Modeling: Methods and Applications".

Links:
1. Paper: https://ieeexplore.ieee.org/document/7858366 
2. Tech Report: https://www2.eecs.berkeley.edu/Pubs/TechRpts/2018/EECS-2018-66.html 

====
Using the package.
====
I have only tried setting up and running this package on Linux. It relies on
softlinks and some other features which might not be readily available on
Windows.

1. Get the most recent version of this package from github.com

            $git clone https://github.com/architgupta93/STEAM 

2. Get the submodules needed for polynomial interpolation. I developed them as
   stand-alone packages that could be used on their own and maintain separate
   version control for them.

            $git submodule update --init

   This should fetch two pacakges (as of 2019/03/13), polynomial-interpolation,
   and device-models. The first, as it name states, is a polynomial
   interpolation package, implementing splines, Lagrange and
   Barycentric-Lagrange interpolants.

3. Setup MAPP. Run the following to setup MAPP

            $autoconf
            $./configure
            $make

4. Open MATLAB (Octave might work too, but it doesn't work very well with
   MATLAB's class based programs.)

            >> start_MAPP

5. Goto the directory containing STEAM examples. It should be under
   STEAM/examples. Open/Run any of the scripts to see how STEAM and MAPP work.

            >> STEAM_demo.m

