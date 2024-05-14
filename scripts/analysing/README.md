This folder contains two scripts used to analyse the data produced by AURORA.jl
simulations:

- `Downsampling_fluxes.jl`: Downsamples the electron fluxes saved in the
`IeFlickering-*.mat` files. Useful when data needs to be transfered between 
computers and time resolution can be coarser.

- `Calculate_densities.jl`: Reads the `IeFlickering-*.mat` files and do some
operations on the electron fluxes to calculate the electron densities.