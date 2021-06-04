**Bjorken VH (c) Mike McNelis**

Created on 11/12/2017 by Mike McNelis\
Last edited on 6/4/2021 by Mike McNelis

## Summary

Anisotropic hydrodynamics for a quark-gluon plasma subject to (0+1)-dimensional Bjorken expansion.

This C++ code evolves nonconformal anisotropic hydrodynamics for a Bjorken expanding system with a QCD equation of state. We used this code to compare anisotropic and viscous hydrodynamics in the paper 

    M. McNelis, D. Bazow and U. Heinz, Phys. Rev. C97, 054912 (2018)


## Running the code

To compile and run the code, do

    sh bjorken_vah.sh
    
The output is written to `results`. You can transfer the data for plotting by doing

    sh copy_results.sh equilibrium    # or glasma
   
You can compare the results to second-order viscous hydrodynamics (computed by [Bjorken VH](https://github.com/mjmcnelis/bjorken_vhydro-)) with the Mathematica notebooks in `plot`.

Note: the `GLASMA` macro in `bjorken.cpp` controls the Bjorken initial conditions.

