# Helper to regenerate reference results after voluntary numerical breaking changes.
# Run this script from the AURORA.jl package root after any change that intentionally
# alters numerical output.  The generated NC files are then committed to the repository
# and used by test/regression/test_regression.jl for future comparisons.

# ======================================================================================== #
#                                  STEADY STATE                                            #
# ======================================================================================== #

using AURORA
altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
E_max = 500;                   # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

msis_file = "test/regression/reference_results/msis_20051008-2200_70N-19E.txt"
iri_file = "test/regression/reference_results/iri_20051008-2200_70N-19E.txt"

## Build the model
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

## Define where to save the results
savedir = mktempdir()

## Define input flux
flux = InputFlux(FlatSpectrum(1e-2; E_min=E_max - 100); beams=1:2)

## Run the simulation
sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())
initialize!(sim; force_recompute=true)
run!(sim)

## Analyze the results
make_volume_excitation_file(sim)

## Overwrite the reference results
mkpath("test/regression/reference_results/SS")
source_file = joinpath(savedir, "analysis", "volume_excitation.nc")
dest_file = "test/regression/reference_results/SS/volume_excitation.nc"
cp(source_file, dest_file; force = true)
println("Saved SS reference to $dest_file")




# ======================================================================================== #
##                                 TIME DEPENDENT                                          #
# ======================================================================================== #

using AURORA
altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
E_max = 500;                   # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

msis_file = "test/regression/reference_results/msis_20051008-2200_70N-19E.txt"
iri_file = "test/regression/reference_results/iri_20051008-2200_70N-19E.txt"

## Build the model
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

## Define where to save the results
savedir = mktempdir()

## Define input flux
flux = InputFlux(FlatSpectrum(1e-2; E_min=100.0), SinusoidalFlickering(5.0);
                 beams=1, z_source=1000.0)

## Run the simulation
sim = AuroraSimulation(model, flux, savedir;
                       mode=TimeDependentMode(duration = 0.2, dt = 0.01,
                                              CFL_number = 128, n_loop = 2))
initialize!(sim; force_recompute=true)
run!(sim)

## Analyze the results
make_volume_excitation_file(sim)

## Overwrite the reference results
mkpath("test/regression/reference_results/TD")
source_file = joinpath(savedir, "analysis", "volume_excitation.nc")
dest_file = "test/regression/reference_results/TD/volume_excitation.nc"
cp(source_file, dest_file; force = true)
println("Saved TD reference to $dest_file")
