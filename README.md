# This is NHDS - The New Hampshire Dispersion relation Solver

Copyright (C) 2025 - Daniel Verscharen (d.verscharen@ucl.ac.uk)\
Mullard Space Science Laboratory, University College London\
Space Science Center, University of New Hampshire

NHDS solves the plasma dispersion relation. In its standard mode, it solves the hot-plasma dispersion relation based on the Vlasov-Maxwell set of equations assuming the plasma consists of a combination of species with drifting bi-Maxwellian background distributions. NHDS can also solve the cold-plasma dispersion relation (or treat the susceptibilities of individual species with cold-plasma theory).


## Requirements

NHDS has the following requirements:

- make
- Fortran 90 compiler
- HDF5 (although it can now also be compiled without HDF5)

## Compiling NHDS

Compile the code by executing
```
    ./configure                                  (if you do not want to use HDF5 capabilities)
    ./configure --with-hdf5=/path/to/hdf5        (this could be: "/usr/local/hdf5" on a Mac)

    make
```
in the main directory. Then run NHDS with
```
    ./src/NHDS <input_file.in>
```

It is possible to install NHDS in the main binary directory through
```
    sudo make install
```

The plasma and numerical parameters are set in the file <input_file.in>.
The code requires only the folder src in the main directory.
