# Helper to regenerate reference results after voluntary numerical breaking changes.
# Run this script from the AURORA.jl package root after any change that intentionally
# alters numerical output:
#
#     julia --project test/regression/update_reference_results.jl
#
# It reruns the reference simulations, packages the results as a tarball, uploads it as a
# new asset on the "test-data" GitHub release (requires the `gh` CLI, authenticated with
# push access, i.e. must be done by a maintainer), and rewrites/updates the
# test/regression/Artifacts.toml to point at it. Commit the Artifacts.toml change to make
# the new baselines effective. Assets are append-only: never overwrite an existing one,
# or older commits will fail their checksum.

using AURORA
using Dates
using SHA
import Pkg: GitTools, PlatformEngines

repo = "egavazzi/AURORA.jl"
release_tag = "regression-data"
stage_dir = mktempdir()

# The simulation parameters below must stay in sync with test/regression/test_regression.jl
altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
E_max = 500;                    # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

msis_file = joinpath(@__DIR__, "input_data", "msis_20051008-2200_70N-19E.txt")
iri_file = joinpath(@__DIR__, "input_data", "iri_20051008-2200_70N-19E.txt")

# Steady state
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
savedir = mktempdir()
flux = InputFlux(FlatSpectrum(1e-2; E_min=E_max - 100); beams=1:2)
sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())
initialize!(sim; force_recompute=true)
run!(sim)
make_volume_excitation_file(sim)
make_column_excitation_file(sim)
mkpath(joinpath(stage_dir, "SS"))
cp(joinpath(savedir, "analysis", "volume_excitation.nc"), joinpath(stage_dir, "SS", "volume_excitation.nc"))
cp(joinpath(savedir, "analysis", "column_excitation.nc"), joinpath(stage_dir, "SS", "column_excitation.nc"))
println("Generated SS references")

# Time dependent
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
savedir = mktempdir()
flux = InputFlux(FlatSpectrum(1e-2; E_min=100.0), SinusoidalFlickering(5.0);
                 beams=1, z_source=1000.0)
sim = AuroraSimulation(model, flux, savedir;
                       mode=TimeDependentMode(duration = 0.2, dt = 0.01,
                                              CFL_number = 128, n_loop = 2))
initialize!(sim; force_recompute=true)
run!(sim)
make_volume_excitation_file(sim)
make_column_excitation_file(sim)
mkpath(joinpath(stage_dir, "TD"))
cp(joinpath(savedir, "analysis", "volume_excitation.nc"), joinpath(stage_dir, "TD", "volume_excitation.nc"))
cp(joinpath(savedir, "analysis", "column_excitation.nc"), joinpath(stage_dir, "TD", "column_excitation.nc"))
println("Generated TD references")

# Package, upload as release asset, and update Artifacts.toml
tree_hash = bytes2hex(GitTools.tree_hash(stage_dir))
asset = "reference_results_$(Dates.format(today(), "yyyy-mm-dd"))_$(first(tree_hash, 8)).tar.gz"
tarball = joinpath(mktempdir(), asset)
PlatformEngines.package(stage_dir, tarball)
sha256_hex = bytes2hex(open(io -> SHA.sha256(io), tarball))

run(`gh release upload $release_tag $tarball --repo $repo`)

write(joinpath(@__DIR__, "Artifacts.toml"), """
    [reference_results]
    git-tree-sha1 = "$tree_hash"
    lazy = true

        [[reference_results.download]]
        url = "https://github.com/$repo/releases/download/$release_tag/$asset"
        sha256 = "$sha256_hex"
    """)
println("Uploaded $asset and updated test/regression/Artifacts.toml — commit the Artifacts.toml change")
