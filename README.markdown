# General many-body library

This general many-body library enables you to calculate the energy (and energy derived properties) to any order when you interface it properly to a quantum chemistry code.

## Requirements

Obviously, any quantum chemistry code that can provide SCF capability. Furthermore, either this quantum chemistry code or an external library to provide the integrals are needed. An example of such a library would be the [GEN1INT library](https://repo.ctcc.no/projects/gen1int "library") (a [paper describing GEN1INT](http://onlinelibrary.wiley.com/doi/10.1002/qua.22886/abstract "paper") is also available).

## FAQ

### What codes are this library interfaced to?

Currently, only the [DALTON](http://www.daltonprogram.org/ "DALTON") program in a development version. The expected public release is sometime in 2014. We do have a plan to provide a patch to the 2013 release before the end of this year which will enable you to try it out yourself.

### How does it work?

Because mblib is a generalized library, it does *not* care too much about the details quantum chemistry. What it cares about is the distribution of quantum chemistry jobs which in turn is run by a code such as DALTON. Properties from this quantum chemistry code are then stored appropriately in mblib for later use.

Please **NOTE** that the mblib API is still young and should not be considered stable because changes can be made as more quantum chemistry codes are added.
