# =============================================================================
# control_script_convert_psd.jl
#
# Batch conversion of AURORA IeFlickering-NN.mat files into phase-space density.
#
# For each file it produces:
#   psd-NN.mat  containing:
#     F        [Nvpar, Nz, Nt]      reduced distribution F(v‖) [s m⁻⁴]
#     vpar_centres [Nvpar]          v‖ bin centres [m/s]
#     v        [nE]                 speed at energy bin centres [m/s]
#     v_par    [n_μ, nE]            v‖ at each (pitch, energy) bin centre [m/s]
#     v_perp   [n_μ, nE]            v⊥ at each (pitch, energy) bin centre [m/s]
#     dE_J     [nE]                 energy bin widths [J]
#     dOmega   [n_μ]                solid-angle bin widths [sr]
#     mu_center [n_μ]               average cosine per pitch bin
#     E        [nE]                 energy left-edge grid [eV]
#     t_run    [Nt]                 time [s]
#     h_atm    [Nz]                 altitude [m]
#     mu_lims  [n_μ+1]             cosine pitch-angle limits
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
output_subdir = "psd"

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

        res = make_psd_from_AURORA(filepath; compute=:F_only)

        matopen(outfile, "w") do io
            write(io, "F",            res.F)            # [Nvpar, Nz, Nt]
            write(io, "vpar_centres", res.vpar_centres) # [Nvpar]
            write(io, "v",            res.v)
            write(io, "v_par",        res.v_par)
            write(io, "v_perp",       res.v_perp)
            write(io, "dE_J",         res.ΔE_J)
            write(io, "dOmega",       res.ΔΩ)
            write(io, "mu_center",    res.μ_center)
            write(io, "E",            res.E)
            write(io, "t_run",        res.t_run)
            write(io, "h_atm",        res.h_atm)
            write(io, "mu_lims",      res.mu_lims)
        end
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
@time res = make_psd_from_AURORA("data/backup/20251104-1634/IeFlickering-01.mat"; vpar_edges);
@time res = make_psd_from_AURORA("data/backup/20251104-1634/IeFlickering-01.mat"; compute = :F_only);
@profview res = make_psd_from_AURORA("data/backup/20251104-1634/IeFlickering-01.mat"; compute = :F_only);
##
# Access fields:
#   res.f         → [n_μ, Nz, Nt, nE]   phase-space density
#   res.v_par     → [n_μ, nE]            v‖ coordinates
#   res.v_perp    → [n_μ, nE]            v⊥ coordinates
#   res.F               → [Nvpar, Nz, Nt]  [s m⁻⁴]
#   res.vpar_centres    → [Nvpar]          [m/s]

## Example plot of f(v‖, v⊥) at one altitude and time
using CairoMakie
iz = 50  # altitude index
it = 40  # time index

# f(v‖, v⊥) at one altitude and time
fig, ax, hm = voronoiplot(res.v_perp[:, :], res.v_par[:, :],
                      res.f[:, iz, it, :];
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
            title="F(v‖), z=$(round(res.h_atm[iz]/1e3, digits=1)) km, t=$(round(res.t_run[it], digits=2)) s",
            yscale=log10)
l = scatter!(res.vpar_centres, res.F[:, iz, it], markersize=4)
ylims!(maximum(res.F[:, iz, it]) / 1e8, maximum(res.F[:, iz, it]) * 10)
display(fig)


##
using MAT

folder = "/nfs/revontuli/data/etienne/paper-1/simulating/results/aw1/psd/"
psd_files = sort(filter(f -> endswith(f, ".mat"), readdir(folder, join=true)))

# Load all files and concatenate F along the time dimension (3rd dim)
chunks = [matread(f) for f in psd_files]
F      = cat([c["F"] for c in chunks]...; dims=3)        # [Nvpar, Nz, Nt_total]
vpar   = chunks[1]["vpar_centres"]                        # [Nvpar]
z_grid = chunks[1]["h_atm"]                              # [Nz]
t_grid = vcat([c["t_run"] for c in chunks]...)           # [Nt_total]

##
iz = 140  # altitude index
it = 280   # time index

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], xlabel="v‖ [m/s]", ylabel="F(v‖) [s m⁻⁴]",
            title="F(v‖), z=$(round(z_grid[iz]/1e3, digits=1)) km, t=$(round(t_grid[it], digits=2)) s",
            yscale=log10)
l = scatter!(vpar, F[:, iz, it], markersize=4)
ylims!(maximum(F[:, iz, it]) / 1e8, maximum(F[:, iz, it]) * 10)
display(fig)



##
using AURORA
using Printf
using MAT

include(joinpath(@__DIR__, "convert_Ie_to_psd.jl"))


file = "data/Visions2/Alfven_train_335s_pchip/IeFlickering-20.mat"
data = load_Ie(file)

vpar_edges = range(0, v_of_E.(data.E[end]), length=100)
vpar_edges = vcat(-reverse(vpar_edges), vpar_edges[2:end])  # symmetric edges around zero
res = make_psd_from_AURORA(file, compute=:F_only, vpar_edges=vpar_edges)


conserved_density = dropdims(sum(data.Ie ./ reshape(res.v, 1, 1, 1, :); dims=(2, 4)), dims=(2, 4))
reduced_density = dropdims(sum(res.F .* reshape(res.Δvpar, :, 1, 1); dims=1), dims=1)

##
i_t = 40
i_z = 300

using CairoMakie
f = Figure(size = (1600, 1200))
ax = Axis(f[1, 1], xlabel="v‖ [m/s]", ylabel="F(v‖) [s m⁻⁴]",
            title="F(v‖), z=$(round(res.h_atm[i_z]/1e3, digits=1)) km, t=$(round(res.t_run[i_t], digits=2)) s",
            yscale=log10)
scatter!(res.vpar_centres, res.F[:, i_z, i_t])
f



comparison = matread("data/Visions2/Alfven_train_335s_pchip/psd/psd-20.mat")
scatter!(comparison["vpar_centres"], comparison["F"][:, i_z, i_t] .* m_e, color=:red)
f
