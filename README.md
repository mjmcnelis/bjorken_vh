**Bjorken VH (c) Mike McNelis**

Created on 11/12/2017 by Mike McNelis\
Last edited on 6/4/2021 by Mike McNelis

## Summary

Second-order viscous hydrodynamics for a quark-gluon plasma subject to (0+1)-dimensional Bjorken expansion.

This repository is complimentary to the [bjorken_vah](https://github.com/mjmcnelis/bjorken_vah) code used in the paper

    M. McNelis, D. Bazow and U. Heinz, Phys. Rev. C97, 054912 (2018)


## Running the code

To compile and run the code, do

    sh bjorken_vh.sh
    
The output is written to `results`. 
   
The macros `GLASMA` and `KINETIC` in `bjorken.cpp` control the Bjorken initial conditions and kinetic theory model, respectively.

