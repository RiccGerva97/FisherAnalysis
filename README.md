# Fisher calc

This code is created to evaluate the Fisher Information Matrix in its approximated form.

## Input data requirements

The files containing the coeffients are axpected to be in .hdf5 extention, created with the h5py Python library.
The structure of the dataset names must be: ```{cosmo}_{i}_{type}```, where:

- ```{cosmo}```: is the cosmology model abbreviation
- ```{i}```: is the number of the realization
- ```{type}```: is the name of the data;
    + `wst`: if wst coefficient
    + `mono`: if monopole of power spectrum
    + `quad`: if quadrupole of power spectrum
    + `hexa`: if hexadecapole of power spectrum
    + `k_coord`: if k coordinates of power spectrum
    + `sn`: if shot noise of power spectrum

NOTE: when using different estimator, be aware to modify all the `if`/`elif` conditions while reading the file and also the corrisponding arrays.

## Descriprion of folders content

### dicts

Here's a collection of useful dictionaries and arrays to describe the parameters that are used in the Quijote simulations.
The dictionaries are meant to make it easier to put and get information about a cosmological model or parameter always with the same index.

### functions

#### math

The algorithms are:

- inverse matrix

- Hartlap formalism

- gaussian

#### parsers

Get useful info from file name (only if it follows the same formalism).

#### pcg

Implements the Melissa O'Neil prng. You can use the class (and its method) PCG to get random numbers or use the `randind` function to get a list of integers (that starts from 0)

#### phys

Implements a dumb periodic condition (`PacMan`). Implements the evaluation of Hubble parameter as a function of redshift.

#### plotting

There are two functions to get the width, height, and angle of a generic ellipse from its covariance matrix and the desired confidence level. The results are quite identical; the difference is that one uses algorithms to get eigenvalues, the other uses analytical formulae.

## main files

### main.py

Main code to evaluate the inverse of the FIM if two estimators are used in real space; there's anoption to evaluate togheter the redshift space.

### main_slim.py

The simplier version of main: no space info and only one estimator. It is assumed there is not a structure in 'multipoles' such as with the power spectrum.

### main_pk.py

As the main.py file, but only about power spectrum, no info about space.
The commented lines are about the multipoles stacked into a single array (useful for redshift space).
