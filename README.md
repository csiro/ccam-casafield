# Casafield (CASA carbon cycle for CCAM)

Casafield is used to create carbon cycle datasets for the conformal cubic grid
employed by the Conformal Cubic Atmospheric Model (CCAM).

## Website

For documentation, see our website at

[https://confluence.csiro.au/display/CCAM/CCAM]

## Dependencies

Casafield requires the NetCDF C library.

## Building casafield

To build casafield with intel, gnu and cray fortran compiler use

```
make
make GFORTRAN=yes
make CRAY=yes
```

Debugging is also enabled with

```
make TEST=yes
```
