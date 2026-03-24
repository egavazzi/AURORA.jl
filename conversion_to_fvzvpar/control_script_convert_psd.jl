# =============================================================================
# control_script_convert_psd.jl
#
# Batch conversion of AURORA IeFlickering-NN.mat files into phase-space density.
#
# For each file it produces:
#   psd-NN.mat containing the common grids and, depending on `compute_mode`,
#   either or both of:
#     f            [Nz, n_μ, Nt, nE] full phase-space density [s^3 m^-6]
#     F            [Nvpar, Nz, Nt]   reduced distribution F(v‖) [s m^-4]
#     vpar_edges   [Nvpar+1]         v‖ bin edges [m/s]
#     vpar_centres [Nvpar]           v‖ bin centres [m/s]
#     dvpar        [Nvpar]           v‖ bin widths [m/s]
#
#   It always writes:
#     v         [nE]                 speed at energy bin centres [m/s]
#     v_par     [n_μ, nE]            v‖ at each (pitch, energy) bin centre [m/s]
#     v_perp    [n_μ, nE]            v⊥ at each (pitch, energy) bin centre [m/s]
#     dE_J      [nE]                 energy bin widths [J]
#     dOmega    [n_μ]                solid-angle bin widths [sr]
#     mu_center [n_μ]                average cosine per pitch bin
#     E         [nE]                 energy left-edge grid [eV]
#     t_run     [Nt]                 time [s]
#     h_atm     [Nz]                 altitude [m]
#     mu_lims   [n_μ+1]              cosine pitch-angle limits
# =============================================================================

using AURORA
using Printf
using MAT

include(joinpath(@__DIR__, "convert_Ie_to_psd.jl"))

# ---------------------------------------------------------------------------
# Configuration — edit these paths before running
# ---------------------------------------------------------------------------
# main_dir = "/nfs/revontuli/data/etienne/paper-1/simulating/results/"
main_dir = "data/Visions2/"

# Select subdirectories to process (adjust the filter string as needed)
dir_filter = "335s_pchip"

# Pattern matching simulation output files
file_pattern = r"IeFlickering\-[0-9]+\.mat"

# Subdirectory name inside each run directory where results are saved
output_subdir = "psd_new"

# What to compute and save: :f_only, :F_only, or :both
compute_mode = :F_only

# Optional custom v_parallel grid for the F reduction.
data = matread("data/Visions2/Alfven_train_335s_pchip/IeFlickering-01.mat")
vpar_edges = range(0, v_of_E.(data["E"][end]), length=100)
vpar_edges = vcat(-reverse(vpar_edges), vpar_edges[2:end])  # symmetric edges around zero
# vpar_edges = nothing


function write_psd_result(outfile::AbstractString, res)
    matopen(outfile, "w") do io
        if hasproperty(res, :f)
            write(io, "f", res.f)
        end

        if hasproperty(res, :F)
            write(io, "F",            res.F)
            write(io, "vpar_edges",   res.vpar_edges)
            write(io, "vpar_centres", res.vpar_centres)
            write(io, "dvpar",        res.Δvpar)
        end

        write(io, "v",         res.v)
        write(io, "v_par",     res.v_par)
        write(io, "v_perp",    res.v_perp)
        write(io, "dE_J",      res.ΔE_J)
        write(io, "dOmega",    res.ΔΩ)
        write(io, "mu_center", res.μ_center)
        write(io, "E",         res.E)
        write(io, "t_run",     res.t_run)
        write(io, "h_atm",     res.h_atm)
        write(io, "mu_lims",   res.μ_lims)
    end
end

# ---------------------------------------------------------------------------
# Processing loop
# ---------------------------------------------------------------------------
directories = readdir(main_dir, join=true)
directories_to_process = filter(d -> isdir(d) && contains(basename(d), dir_filter), directories)

for current_dir in directories_to_process
    println("Processing $(basename(current_dir))...")

    files = readdir(current_dir, join=true)
    input_files = filter(f -> !isnothing(match(file_pattern, basename(f))), files)
    input_files = sort(input_files, by=f -> parse(Int, match(r"(\d+)", basename(f)).captures[1]))

    if isempty(input_files)
        @warn "  No IeFlickering files found in $current_dir, skipping."
        continue
    end

    savedir = joinpath(current_dir, output_subdir)
    isdir(savedir) || mkdir(savedir)

    for (j, filepath) in enumerate(input_files)
        outfile = joinpath(savedir, @sprintf("psd-%02d.mat", j))
        println("  $(basename(filepath)) → $(basename(outfile))")

        res = make_psd_from_AURORA(filepath; compute=compute_mode, vpar_edges=vpar_edges)
        write_psd_result(outfile, res)
    end
end

println("Done.")




# =============================================================================
# Interactive / scratch section
# =============================================================================
##
using AURORA

E_grid, _ = make_energy_grid(5e3)
v_grid = v_of_E.(E_grid)
vpar_edges = vcat(-reverse(v_grid), v_grid)
# Quick test on a single file:
include("convert_Ie_to_psd.jl")
@time res = make_psd_from_AURORA("data/backup/20251104-1634/IeFlickering-01.mat"; compute=:both);
@time res_F = make_psd_from_AURORA("data/backup/20251104-1634/IeFlickering-01.mat"; compute=:F_only);
@profview res_F = make_psd_from_AURORA("data/backup/20251104-1634/IeFlickering-01.mat"; compute=:F_only);

## Example plot of f(v‖, v⊥) at one altitude and time
using CairoMakie
iz = 50  # altitude index
it = 40  # time index

# f(v‖, v⊥) at one altitude and time
fig, ax, hm = voronoiplot(res.v_perp[:, :], res.v_par[:, :],
                      res.f[iz, :, it, :];
                      strokewidth = 0,
                      colorscale = log10,
                      show_generators = false,
                      axis=(xlabel="v‖ [m/s]", ylabel="v⊥ [m/s]",
                            title="f(v‖, v⊥), z=$(round(res.h_atm[iz]/1e3, digits=1)) km"))
Colorbar(fig[1,2], hm)
display(fig)


## Example plot of F(v‖) at one altitude and time
iz = 350  # altitude index
it = 100   # time index

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], xlabel="v‖ [m/s]", ylabel="F(v‖) [s m⁻⁴]",
            title="F(v‖), z=$(round(res_F.h_atm[iz]/1e3, digits=1)) km, t=$(round(res_F.t_run[it], digits=2)) s",
            yscale=log10)
l = scatter!(res_F.vpar_centres, res_F.F[:, iz, it], markersize=4)
ylims!(maximum(res_F.F[:, iz, it]) / 1e8, maximum(res_F.F[:, iz, it]) * 10)
display(fig)
