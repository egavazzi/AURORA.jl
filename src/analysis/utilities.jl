using MAT: matopen

# ======================================================================================== #
#                              RESULT FILE HELPERS                                       #
# ======================================================================================== #

"""
    list_result_files(sim_dir) → Vector{String}

Return the sorted list of `IeFlickering-*.mat` result files in `sim_dir`, ordered
numerically by loop index (not lexicographically).

# Example
```julia
files = list_result_files("data/my_run/")
```
"""
function list_result_files(sim_dir::AbstractString)
    files = readdir(sim_dir; join=true)
    result_files = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    sort!(result_files,
          by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]))
    return result_files
end

"""
    read_result(file) → NamedTuple

Read one `IeFlickering-*.mat` result file produced by AURORA and return a named tuple
with all stored fields.

# Returns
Named tuple with fields:
- `Ie_ztE`         : electron flux `[n_z·n_μ, n_t, n_E]` (#e⁻/m²/s)
- `t_run`          : time axis for this file `[n_t]` (seconds)
- `h_atm`          : altitude grid `[n_z]` (metres)
- `E_centers`      : energy bin centres `[n_E]` (eV)
- `E_edges`        : energy bin edges `[n_E+1]` (eV)
- `dE`             : energy bin widths `[n_E]` (eV)
- `mu_lims`        : pitch-angle cosine limits `[n_μ+1]`
- `I0`             : state at the end of this loop (initial condition for next loop)
- `solver_type`    : `"time_dependent"` or `"steady_state"`
- `mode_type`      : `"time_dependent"`, `"steady_state_multi_step"`, or `"steady_state_single_step"`
- `mu_scatterings` : Dict with scattering geometry (`P_scatter`, `BeamWeight`, etc.)

# Example
```julia
result = read_result("data/my_run/IeFlickering-01.mat")
result.t_run      # time axis
result.Ie_ztE     # flux array
```
"""
# MAT.jl reads a 1×1 MATLAB matrix as a scalar Float64; wrap it back into a vector.
MATvec(x::AbstractArray) = vec(x)
MATvec(x::Number)        = [Float64(x)]

function read_result(file::AbstractString)
    f = matopen(file)
    Ie_ztE         = read(f, "Ie_ztE")
    t_run          = MATvec(read(f, "t_run"))
    h_atm          = vec(read(f, "h_atm"))
    E_centers      = vec(read(f, "E_centers"))
    E_edges        = vec(read(f, "E_edges"))
    dE             = vec(read(f, "dE"))
    mu_lims        = vec(read(f, "mu_lims"))
    I0             = read(f, "I0")
    solver_type    = read(f, "solver_type")
    mode_type      = read(f, "mode_type")
    mu_scatterings = read(f, "mu_scatterings")
    close(f)
    return (
        Ie_ztE         = Ie_ztE,
        t_run          = t_run,
        h_atm          = h_atm,
        E_centers      = E_centers,
        E_edges        = E_edges,
        dE             = dE,
        mu_lims        = mu_lims,
        I0             = I0,
        solver_type    = solver_type,
        mode_type      = mode_type,
        mu_scatterings = mu_scatterings,
    )
end
