# Troubleshooting

## Segmentation fault

If you get a segmentation fault when the code is calling Matlab (when the 
cross-section are being loaded), you want to check what version of Matlab is
being called by Julia. By default, the package `MATLAB.jl` which is taking care 
of calling Matlab is using the Matlab installation with highest version number 
installed on your machine.

The problem is that some of the newer version of Matlab don't work well with the 
package `MATLAB.jl` and cause segmentation fault when being called. If that is 
the case, you need to switch back to a Matlab version with a version number 
older than R2022a. To do this, you can follow the steps at <https://github.com/JuliaInterop/MATLAB.jl?tab=readme-ov-file#changing-matlab-version>.

Don't forget to login and activate your Matlab license on that older version 
before trying to call it from Julia.