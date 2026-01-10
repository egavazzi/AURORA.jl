# Troubleshooting

## IRI and MSIS data
You might run into some issues when trying to load/generate some msis or iri 
data for the first time. This is because AURORA uses the `pymsis` and `iri2016` 
python packages as dependencies. These packages use in turn some fortran
dependencies. If your system is complaining about missing fortran compiler or
missing cmake install, try running the following commands in a terminal 
(on Linux):

```
$> sudo apt update
$> sudo apt install cmake
$> sudo apt install gfortran
```