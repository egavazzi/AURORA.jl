using HCubature: hcubature, hcubature_buffer

# ======================================================================================== #
#           CASCADING SPEC — Species specific cascading physics definition                 #
# ======================================================================================== #


struct CascadingSpec{F}
    name::String
    ionization_thresholds::Vector{Float64}
    n_secondaries::Vector{Int}   # secondary e⁻ ejected per channel: 1 = single, 2 = double ionization
    secondary_law::F   # callable (E_secondary, E_primary) -> Float64
    function CascadingSpec(name::AbstractString, thresholds::AbstractVector, law;
                           n_secondaries::AbstractVector{<:Integer} = ones(Int, length(thresholds)))
        require_reproducible(law, "cascading secondary_law")
        length(n_secondaries) == length(thresholds) || throw(ArgumentError(
            "n_secondaries (length $(length(n_secondaries))) must match ionization_thresholds \
             (length $(length(thresholds)))"))
        all(n -> 1 <= n <= 2, n_secondaries) || throw(ArgumentError(
            "n_secondaries entries must be 1 (single) or 2 (double ionization); got $(collect(n_secondaries))"))
        return new{typeof(law)}(String(name), collect(Float64, thresholds),
                                collect(Int, n_secondaries), law)
    end
end

# Secondary-electron distribution for atomic O. It needs external parameters (A and B) so we
# cannot just wrap it in @law. Instead we build a custom callable struct that contains the
# parameters in its fields.
struct OSecondaryLaw
    energy_params::Vector{Float64}
    A_params::Vector{Float64}
    B_params::Vector{Float64}
end

function (law::OSecondaryLaw)(E_s, E_p)
    A_factor = interp_flat(law.energy_params, law.A_params, E_p) * 1.25  # empirical correction
    B_factor = interp_flat(law.energy_params, law.B_params, E_p)
    return B_factor / (1 + (E_s / A_factor)^(5 / 3))
end

# Small custom linear interpolator, with flat extrapolation beyond the x range.
function interp_flat(x::AbstractVector, y::AbstractVector, xq)
    xq <= x[1]   && return float(y[1])
    xq >= x[end] && return float(y[end])
    i = searchsortedlast(x, xq)
    t = (xq - x[i]) / (x[i + 1] - x[i])
    return y[i] + t * (y[i + 1] - y[i])
end

function DefaultCascadingSpecN2()
    ionization_thresholds = [15.581, 16.73, 18.75, 24.0, 42.0]
    n_secondaries         = [1, 1, 1, 1, 2]
    law = @law (E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)
    return CascadingSpec("N2", ionization_thresholds, law; n_secondaries)
end

function DefaultCascadingSpecO2()
    ionization_thresholds = [12.072, 16.1, 16.9, 18.2, 18.9, 32.51]
    n_secondaries         = [1, 1, 1, 1, 1, 2]
    law = @law (E_s, E_p) -> 1.0 / (15.2^2 + E_s^2)
    return CascadingSpec("O2", ionization_thresholds, law; n_secondaries)
end

function DefaultCascadingSpecO()
    ionization_thresholds = [13.618, 16.9, 18.6, 28.5]
    n_secondaries         = [1, 1, 1, 2]
    energy_params = [100.0, 200, 500, 1000, 2000]  # eV
    A_params = [12.6, 13.7, 14.1, 14.0, 13.7]
    B_params = [7.18, 4.97, 2.75, 1.69, 1.02] .* 1e-22
    law = OSecondaryLaw(energy_params, A_params, B_params)
    return CascadingSpec("O", ionization_thresholds, law; n_secondaries)
end


# ======================================================================================== #
#                      CASCADING CACHE — Per-species in-memory cache                       #
# ======================================================================================== #

# Species specific cascading cache container
mutable struct SpeciesCascadingCache{S<:CascadingSpec}
    spec::S
    primary_transfer_matrix::Array{Float64, 3}
    secondary_transfer_matrix::Array{Float64, 3}
    E_edges::Vector{Float64}
    ionization_thresholds::Vector{Float64}
end

# Initialization constructor
function SpeciesCascadingCache(spec::S) where {S<:CascadingSpec}
    return SpeciesCascadingCache{S}(spec, zeros(0, 0, 0), zeros(0, 0, 0), Float64[], Float64[])
end








# ======================================================================================== #
#                                CALCULATION FUNCTIONS                                     #
# ======================================================================================== #

#=
The implementation here can be a bit hard to follow, so here are some explanations.

For each primary energy bin (outer loop), we want to calculate the distribution of
probabilities of degraded primary. To do this, we loop over allowed degraded primary bins
(inner loop). For each degraded primary bin, we fix a degraded primary energy inside of it,
and integrate over the allowed primary energies that can produce this degraded primary energy.
For each of these allowed primary energies and fixed degraded primary energy, the secondary
energy is determined by energy conservation, and we can evaluate the secondary distribution
law for that specific secondary energy to get the probability of this specific combination
of primary and degraded primary energies. Then we fix another degraded primary energy inside
the degraded bin, and repeat until we have covered the whole degraded primary bin.
Mathematically, all of this translates into a double integral.

We do the same thing for the secondary transfer matrix, except that we now inner loop over
allowed secondary bins, and for each fixed secondary energy we integrate over the allowed
primary energies that can produce this secondary energy.
=#

"""
    DegradedCascadingIntegrand(E_primary_bin_min, E_primary_bin_max, threshold, secondary_law)

Callable integrand used by `hcubature` to integrate the degraded-primary transfer matrix.
"""
struct DegradedCascadingIntegrand{FT, F}
    E_primary_bin_min::FT
    E_primary_bin_max::FT
    threshold::FT
    secondary_law::F
end

"""
    DegradedCascadingIntegrand(E_primary_bin_min, E_primary_bin_max, threshold, secondary_law)

Callable integrand used by `hcubature` to integrate the secondary transfer matrix.
"""
struct SecondaryCascadingIntegrand{FT, F}
    E_primary_bin_min::FT
    E_primary_bin_max::FT
    threshold::FT
    secondary_law::F
end



function (integrand::DegradedCascadingIntegrand)(vars)
    # integration variables passed by hcubature; u_primary ∈ [0, 1] is the mapped primary-energy coordinate
    E_degraded, u_primary = vars

    # The lower bound of E_primary depends on E_degraded (it must exceed E_degraded + threshold
    # so that the secondary electron energy is positive).
    # The upper bound of E_primary also depends on E_degraded (it cannot exceed
    # 2 * E_degraded + threshold, otherwise the secondary would be more energetic than
    # the degraded primary.
    # This comes from
    #       E_secondary = E_primary - threshold - E_degraded
    # and
    #       E_secondary ≤ E_degraded
    #   =>  E_primary - threshold - E_degraded ≤ E_degraded
    #   =>  E_primary ≤ 2 * E_degraded + threshold
    #
    # This makes the 2D integration domain non-rectangular (trapezoidal).
    # To handle this with hcubature, we map E_primary onto u_primary ∈ [0, 1] where
    # u_primary = 0 corresponds to E_primary_lower and u_primary = 1 corresponds to
    # E_primary_upper. The Jacobian dE_primary/du_primary of this transformation is
    # included below.
    E_primary_lower = max(integrand.E_primary_bin_min, E_degraded + integrand.threshold)
    E_primary_upper = min(integrand.E_primary_bin_max, integrand.threshold + 2 * E_degraded)
    jacobian = E_primary_upper - E_primary_lower
    jacobian <= 0 && return 0.0 # not really necessary as
                                # (E_primary_lower ≤ E_primary_upper)
                                # but one is never too safe

    E_primary = E_primary_lower + u_primary * jacobian
    E_secondary = E_primary - integrand.threshold - E_degraded

    return jacobian * integrand.secondary_law(E_secondary, E_primary)
end

function (integrand::SecondaryCascadingIntegrand)(vars)
    E_secondary, u_primary = vars

    # For a fixed secondary energy, we enforce E_secondary <= E_degraded with
    # E_degraded = E_primary - threshold - E_secondary, which implies
    # E_primary >= threshold + 2 * E_secondary.
    E_primary_lower = max(integrand.E_primary_bin_min,
                          integrand.threshold + 2 * E_secondary)
    E_primary_upper = integrand.E_primary_bin_max
    jacobian = E_primary_upper - E_primary_lower
    jacobian <= 0 && return 0.0

    E_primary = E_primary_lower + u_primary * jacobian

    return jacobian * integrand.secondary_law(E_secondary, E_primary)
end


#=
DOUBLE IONIZATION
-----------------
A double-ionization event ejects TWO secondary electrons, so the excess energy
    W = E_primary - threshold
is shared by three outgoing electrons: the degraded primary E_d and two secondaries E_s1, E_s2,
with E_d + E_s1 + E_s2 = W. The two secondaries are drawn independently from the same `secondary_law`,
so their joint (unnormalized) density is `secondary_law(E_s1, E_p) * secondary_law(E_s2, E_p)`.
By convention the degraded primary is the most energetic electron (E_d ≥ E_s1 and E_d ≥ E_s2),
which restricts the integration to
    R = { E_s1, E_s2 ≥ 0 ,  2·E_s1 + E_s2 ≤ W ,  E_s1 + 2·E_s2 ≤ W }.

These two integrands add a third integration variable (the partner electron) on top of the
E_primary mapping used by the single-ionization integrands. The resulting matrices are built so
that, summed over their respective output bins, both equal the same event count
    Z₂ = ∬_R secondary_law(E_s1)·secondary_law(E_s2) dE_s1 dE_s2 ,
exactly like the single-ionization pair. `compute_ionization_spectra!` then normalizes by
`sum_primary` (≈ Z₂, the degraded primary being the highest-energy electron and thus on-grid)
and multiplies the secondary spectrum by `n_secondary = 2`. With the secondary matrix carrying
the PER-secondary marginal (∫ over the partner), this reproduces ⟨E_d⟩ + 2·⟨E_s⟩ = W exactly,
so energy is conserved with no change to `compute_ionization_spectra!`.

Small caveat: `secondary_law` is a SINGLE-ionization fit. The forms used here (the Cauchy/Lorentzian
(of course I had to place Cauchy) `1 / (a² + E_s²)` for N₂/O₂ and the Opal–Peterson–Beaty form for O)
describe the ejected-electron energy spectrum of *single* ionization. There is no equally well-established
law (to our knowledge) for the two electrons freed in a double-ionization event, so we do the approximation
of treating the two secondaries as independent, identically-distributed draws from the single-ionization law.
This neglects both the (likely different) true double-ionization secondary spectrum and the energy + scattering
correlation between the two ejected electrons. It is the best we can do for now without a dedicated
double-ionization differential cross section.
=#

"""
    DoubleSecondaryCascadingIntegrand(E_primary_bin_min, E_primary_bin_max, threshold, secondary_law)

Integrand for the per-secondary marginal of a double-ionization event. Outer variable is the
binned secondary energy `E_s1`; `u_partner` maps the partner secondary `E_s2`, and `u_primary`
maps `E_primary` (each over [0, 1]).
"""
struct DoubleSecondaryCascadingIntegrand{FT, F}
    E_primary_bin_min::FT
    E_primary_bin_max::FT
    threshold::FT
    secondary_law::F
end

"""
    DoublePrimaryCascadingIntegrand(E_primary_bin_min, E_primary_bin_max, threshold, secondary_law)

Integrand for the degraded-primary distribution of a double-ionization event. Outer variable is
the binned degraded energy `E_d`; `u_secondary` maps the free secondary `E_s1` (the partner
`E_s2 = W - E_d - E_s1` follows from energy conservation), and `u_primary` maps `E_primary`.
"""
struct DoublePrimaryCascadingIntegrand{FT, F}
    E_primary_bin_min::FT
    E_primary_bin_max::FT
    threshold::FT
    secondary_law::F
end

function (integrand::DoubleSecondaryCascadingIntegrand)(vars)
    E_s1, u_partner, u_primary = vars

    # The degraded primary is the most energetic electron, so E_d ≥ E_s1 requires
    #   W - E_s1 - E_s2 ≥ E_s1  (using E_d = W - E_s1 - E_s2)
    # In the worst case E_s2 = 0, leaving E_d = W - E_s1 ≥ E_s1, i.e. W ≥ 2·E_s1, i.e.
    #   E_primary ≥ threshold + 2·E_s1.
    E_primary_lower = max(integrand.E_primary_bin_min, integrand.threshold + 2 * E_s1)
    E_primary_upper = integrand.E_primary_bin_max
    jac_primary = E_primary_upper - E_primary_lower
    jac_primary <= 0 && return 0.0
    E_primary = E_primary_lower + u_primary * jac_primary
    W = E_primary - integrand.threshold

    # With E_s1 fixed, find the valid range for the partner E_s2. Energy conservation gives
    #   E_d = W - E_s1 - E_s2,
    # and the two ordering constraints E_d ≥ E_s1 and E_d ≥ E_s2 translate to:
    #
    #   (1) E_d ≥ E_s1:  W - E_s1 - E_s2 ≥ E_s1  =>  E_s2 ≤ W - 2·E_s1
    #   (2) E_d ≥ E_s2:  W - E_s1 - E_s2 ≥ E_s2  =>  E_s2 ≤ (W - E_s1)/2
    #
    # Both bounds are ≥ 0 because E_primary ≥ threshold + 2·E_s1 guarantees W ≥ 2·E_s1.
    # Which one is tighter depends on W relative to E_s1:
    #   W - 2·E_s1 ≤ (W - E_s1)/2  ⟺  2W - 4·E_s1 ≤ W - E_s1  ⟺  W ≤ 3·E_s1.
    # So constraint (1) wins when the primary is barely above threshold; (2) wins at high energy.
    E_s2_upper = min(W - 2 * E_s1, (W - E_s1) / 2)
    E_s2_upper <= 0 && return 0.0
    E_s2 = u_partner * E_s2_upper   # E_s2 ∈ [0, E_s2_upper], mapped from u_partner ∈ [0, 1]
    jac_partner = E_s2_upper

    return jac_primary * jac_partner *
           integrand.secondary_law(E_s1, E_primary) * integrand.secondary_law(E_s2, E_primary)
end

function (integrand::DoublePrimaryCascadingIntegrand)(vars)
    E_degraded, u_secondary, u_primary = vars

    # E_d is the most energetic electron, so E_d ≥ E_s1 and E_d ≥ E_s2. Since E_s1, E_s2 ≥ 0
    # and E_d + E_s1 + E_s2 = W, the minimum W is achieved when both secondaries are zero:
    #   W_min = E_d  =>  E_primary_lower = threshold + E_d.
    # The maximum W is achieved when both secondaries are as large as possible while still
    # satisfying E_d ≥ E_s1 and E_d ≥ E_s2, i.e. E_s1 = E_s2 = E_d:
    #   W_max = 3·E_d  =>  E_primary_upper = threshold + 3·E_d.
    E_primary_lower = max(integrand.E_primary_bin_min, integrand.threshold + E_degraded)
    E_primary_upper = min(integrand.E_primary_bin_max, integrand.threshold + 3 * E_degraded)
    jac_primary = E_primary_upper - E_primary_lower
    jac_primary <= 0 && return 0.0
    E_primary = E_primary_lower + u_primary * jac_primary
    W = E_primary - integrand.threshold

    # With E_d and E_primary fixed, E_s2 = W - E_d - E_s1 is determined by E_s1 (Jacobian 1).
    # The four constraints on E_s1 are:
    #   E_s1 ≥ 0                (obvious)
    #   E_s2 ≥ 0  =>  E_s1 ≤ W - E_d
    #   E_s1 ≤ E_d              (E_d is the maximum)
    #   E_s2 ≤ E_d  =>  W - E_d - E_s1 ≤ E_d  =>  E_s1 ≥ W - 2·E_d
    # Combined: E_s1 ∈ [max(0, W - 2·E_d),  min(E_d, W - E_d)].
    E_s1_lower = max(0.0, W - 2 * E_degraded)
    E_s1_upper = min(E_degraded, W - E_degraded)
    jac_secondary = E_s1_upper - E_s1_lower
    jac_secondary <= 0 && return 0.0
    E_s1 = E_s1_lower + u_secondary * jac_secondary
    E_s2 = W - E_degraded - E_s1

    return jac_primary * jac_secondary *
           integrand.secondary_law(E_s1, E_primary) * integrand.secondary_law(E_s2, E_primary)
end


"""
    calculate_cascading_matrices(spec, E_edges; verbose=true)

Calculate the energy-degradation transfer matrix for a species defined by its `CascadingSpec`.
The matrices can later be used to directly get the degraded primary and secondary electron
distributions for any given primary energy and ionization threshold.

The outer loop structure is identical for all species. The only species-specific ingredients are
- `spec.secondary_law(E_secondary, E_primary) -> Float64`, which describes how secondary
    electrons distribute in energy given a primary electron at the current integration energy.
- `spec.ionization_thresholds`, which define the ionization thresholds for the species and thus
    the number of transfer matrices to calculate.

Single-ionization channels are integrated with adaptive cubature (`rtol = 1e-4`).
Double-ionization channels use the numerical-CDF method with fixed Gauss-Legendre rules
(see `fill_double_ionization_bin_cdf!`): on rows carrying real weight it agrees with a
tight adaptive reference to a few 1e-4, and its cost stays bounded on large grids where
adaptive integration of the 3-D integrands becomes intractable (E_max ≳ 50 keV). Rows just
above the threshold (~1e-5 of the weight) are less accurate (up to ~1e-1); the test suite
validates the method against the adaptive reference implementation.

# Arguments
- `spec::CascadingSpec` contains species name, ionization thresholds and secondary distribution law
- `E_edges`: Energy grid edges to match (eV)

# Returns
- `(primary_transfer_matrix, secondary_transfer_matrix, E_edges, ionization_thresholds)`:
    degraded-primary transfer matrix [n_E, n_E, n_thresholds], secondary transfer matrix
    [n_E, n_E, n_thresholds], energy grid edges, and ionization thresholds
"""
function calculate_cascading_matrices(spec::CascadingSpec, E_edges; verbose = true,
                                       progress_interval::Real = 10)
    # An ExprLaw's callable is evaluated code (possibly in a newer world age than this
    # frame, e.g. rebuilt during a JLD2 cache reload) behind an untyped field, so calling
    # it through the ExprLaw wrapper is a dynamic call with boxed arguments — very bad with
    # the billions of evaluations the integrations below perform. We unwrap it once and enter
    # the loop body through a single invokelatest: the body runs in the latest world age
    # (where the law's method exists) and every inner law call dispatches statically.
    law = runtime_law(spec.secondary_law)
    return Base.invokelatest(calculate_cascading_matrices, spec, law, E_edges;
                             verbose, progress_interval)
end

# Barrier method: `law` is `spec.secondary_law` with any ExprLaw wrapper already removed,
# so it has a concrete type and the integration loops below specialize on it.
function calculate_cascading_matrices(spec::CascadingSpec, law, E_edges; verbose = true,
                                      progress_interval::Real = 10)
    E_left = @view(E_edges[1:end-1])
    n_E = length(E_left) # number of energy bins is one less than number of edges

    ionization_thresholds = spec.ionization_thresholds
    n_secondaries = spec.n_secondaries
    n_thresholds = length(ionization_thresholds)
    primary_transfer_matrix = zeros(n_E, n_E, n_thresholds)
    secondary_transfer_matrix = zeros(n_E, n_E, n_thresholds)

    verbose && print("Calculating energy-degradation transfer matrices for e⁻ - $(spec.name) ionizing collisions...")

    # Throttle progress updates. Rewrite a TTY line in place, but append log-friendly lines
    # when stdout is redirected.
    t_start = time()
    last_report = Threads.Atomic{Float64}(t_start)
    printed_progress = Threads.Atomic{Bool}(false)
    use_tty = stdout isa Base.TTY

    # Pre-allocate hcubature work buffers, one per thread (heap located). Single-ionization
    # integrands are 2-D. The double-ionization path uses fixed Gauss-Legendre rules and
    # needs no buffers.
    primary_bufs = [hcubature_buffer(DegradedCascadingIntegrand(0.0, 1.0, 0.0, law),
                                     (0.0, 0.0), (1.0, 1.0)) for _ in 1:Threads.maxthreadid()]
    secondary_bufs = [hcubature_buffer(SecondaryCascadingIntegrand(0.0, 1.0, 0.0, law),
                                       (0.0, 0.0), (1.0, 1.0)) for _ in 1:Threads.maxthreadid()]

    # Loop over ionization thresholds
    for (i_pass, i_threshold) in enumerate(n_thresholds:-1:1)
        threshold = ionization_thresholds[i_threshold]
        single_secondary = n_secondaries[i_threshold] == 1

        # Find the first primary bin whose left edge is above the threshold and can
        # therefore contribute to ionization collisions.
        i_min_primary = searchsortedfirst(E_left, threshold)
        n_bins = n_E - i_min_primary + 1
        bins_done = Threads.Atomic{Int}(0)
        pass_printed = Threads.Atomic{Bool}(false)
        status(done, t_now) = "  $(spec.name) cascading: threshold pass $(i_pass)/$(n_thresholds), " *
                              "$(done)/$(n_bins) primary bins ($(round(Int, 100 * done / n_bins)) %), " *
                              "$(round(Int, t_now - t_start)) s elapsed"

        # Loop over primary electron energy bins
        Threads.@threads :static for i_primary in i_min_primary:n_E
            tid = Threads.threadid()
            if single_secondary
                fill_single_ionization_bin!(primary_transfer_matrix, secondary_transfer_matrix,
                                            E_edges, E_left, threshold, i_primary, i_threshold,
                                            law, primary_bufs[tid], secondary_bufs[tid])
            else
                fill_double_ionization_bin_cdf!(primary_transfer_matrix, secondary_transfer_matrix,
                                                E_edges, E_left, threshold, i_primary, i_threshold,
                                                law)
            end
            done = Threads.atomic_add!(bins_done, 1) + 1
            if verbose
                t_now = time()
                t_last = last_report[]
                # Report at most once per interval.
                if t_now - t_last >= progress_interval &&
                   Threads.atomic_cas!(last_report, t_last, t_now) === t_last
                    # The opening "Calculating..." print has no newline; break it once.
                    prefix = Threads.atomic_xchg!(printed_progress, true) ? "" : "\n"
                    Threads.atomic_xchg!(pass_printed, true)
                    if use_tty
                        # \r + clear-line: rewrite the status in place
                        print(prefix * "\r\e[2K" * status(done, t_now))
                    else
                        println(prefix * status(done, t_now))
                    end
                    flush(stdout)
                end
            end
        end

        # Finalize this pass's status as a single completed line (on a TTY this overwrites
        # the in-place updates, leaving exactly one line per threshold pass).
        if verbose && pass_printed[]
            line = status(n_bins, time())
            println(use_tty ? "\r\e[2K" * line : line)
            flush(stdout)
        end
    end

    if verbose
        # If progress lines were printed, the opening print was already line-broken.
        closing = printed_progress[] ? "  $(spec.name) cascading done" : " done"
        println("$closing ($(round(time() - t_start; digits = 1)) s).")
    end
    return primary_transfer_matrix, secondary_transfer_matrix, E_edges, ionization_thresholds
end


# Single-ionization (one secondary) contribution for one primary bin `i_primary` and threshold
# index `i_threshold`.
function fill_single_ionization_bin!(primary_transfer_matrix, secondary_transfer_matrix,
                                     E_edges, E_left, threshold, i_primary, i_threshold,
                                     secondary_law, primary_buf, secondary_buf)
    E_primary_bin_min = E_edges[i_primary]      # left edge
    E_primary_bin_max = E_edges[i_primary + 1]  # right edge

    # Define the integrands for this primary bin and ionization threshold.
    degraded_integrand = DegradedCascadingIntegrand(E_primary_bin_min, E_primary_bin_max,
                                                    threshold, secondary_law)
    secondary_integrand = SecondaryCascadingIntegrand(E_primary_bin_min, E_primary_bin_max,
                                                      threshold, secondary_law)

    # For a fixed primary energy, the secondary/degraded boundary sits at half of the excess
    # energy. Using the lower edge of the primary bin gives the lowest such boundary.
    E_secondary_boundary_lower = (E_primary_bin_min - threshold) / 2
    # First bin that can receive a degraded primary electron (right edge above that boundary).
    i_min_degraded = searchsortedlast(E_left, E_secondary_boundary_lower)
    # If the lowest secondary/degraded boundary lies below the grid minimum, skip this bin.
    i_min_degraded == 0 && return
    # Loop over degraded primary electron energy bins
    for i_degraded in i_min_degraded:(i_primary - 1)
        E_degraded_bin_min = E_edges[i_degraded]
        E_degraded_bin_max = E_edges[i_degraded + 1]
        E_degraded_lower = max(E_degraded_bin_min, E_secondary_boundary_lower)
        E_degraded_upper = min(E_degraded_bin_max, E_primary_bin_max - threshold)
        # Integrate only if limits are physical
        if E_degraded_upper > E_degraded_lower
            result, _ = hcubature(degraded_integrand, (E_degraded_lower, 0.0),
                                 (E_degraded_upper, 1.0);
                                 rtol = 1e-4, buffer = primary_buf)
            primary_transfer_matrix[i_primary, i_degraded, i_threshold] = result
        end
    end

    # Using the upper edge of the primary bin gives the highest secondary/degraded boundary.
    E_secondary_boundary_upper = (E_primary_bin_max - threshold) / 2
    # Last bin that can receive a secondary electron (left edge below that boundary).
    i_max_secondary = searchsortedlast(E_left, E_secondary_boundary_upper)
    # On a coarse grid the boundary could lie above the left edge of the primary bin; clamp.
    i_max_secondary = min(i_max_secondary, i_primary - 1)
    # Loop over the secondary electron energy bins
    for i_secondary in 1:i_max_secondary
        E_secondary_bin_min = E_edges[i_secondary]
        E_secondary_bin_max = E_edges[i_secondary + 1]
        E_secondary_upper = min(E_secondary_bin_max, E_secondary_boundary_upper)
        # Integrate only if limits are physical
        if E_secondary_upper > E_secondary_bin_min
            result, _ = hcubature(secondary_integrand, (E_secondary_bin_min, 0.0),
                                 (E_secondary_upper, 1.0);
                                 rtol = 1e-4, buffer = secondary_buf)
            secondary_transfer_matrix[i_primary, i_secondary, i_threshold] = result
        end
    end
    return
end


# REFERENCE implementation of the double-ionization (two secondaries) contribution for one
# primary bin, via adaptive 3-D cubature of the joint integrands defined above. NOT used by
# the production path (`fill_double_ionization_bin_cdf!` below): its cost per row grows
# without bound on large grids (minutes per row at E_max ≳ 100 keV). It is kept as the
# ground truth that the test suite validates the production path against, since its rtol
# gives it a controlled error for any secondary law.
# Both matrices sum (over their output bins) to the same event count Z₂, so
# `compute_ionization_spectra!` conserves energy unchanged. Bin ranges differ from single
# ionization: the degraded primary is the highest-energy electron, so it reaches down only
# to W/3 (not W/2), while a single secondary still maxes out at W/2.
function fill_double_ionization_bin!(primary_transfer_matrix, secondary_transfer_matrix,
                                     E_edges, E_left, threshold, i_primary, i_threshold,
                                     secondary_law, primary_buf, secondary_buf)
    E_primary_bin_min = E_edges[i_primary]      # left edge
    E_primary_bin_max = E_edges[i_primary + 1]  # right edge

    degraded_integrand = DoublePrimaryCascadingIntegrand(E_primary_bin_min, E_primary_bin_max,
                                                        threshold, secondary_law)
    secondary_integrand = DoubleSecondaryCascadingIntegrand(E_primary_bin_min, E_primary_bin_max,
                                                          threshold, secondary_law)

    # Degraded primary: E_d ≥ W/3, so the lowest reachable degraded energy for this primary bin
    # is (E_primary_bin_min − threshold)/3 and the highest is W = E_primary_bin_max − threshold.
    E_degraded_boundary_lower = (E_primary_bin_min - threshold) / 3
    E_degraded_max = E_primary_bin_max - threshold
    i_min_degraded = max(1, searchsortedlast(E_left, E_degraded_boundary_lower))
    for i_degraded in i_min_degraded:(i_primary - 1)
        E_degraded_bin_min = E_edges[i_degraded]
        E_degraded_bin_max = E_edges[i_degraded + 1]
        E_degraded_lower = max(E_degraded_bin_min, E_degraded_boundary_lower)
        E_degraded_upper = min(E_degraded_bin_max, E_degraded_max)
        if E_degraded_upper > E_degraded_lower
            result, _ = hcubature(degraded_integrand, (E_degraded_lower, 0.0, 0.0),
                                 (E_degraded_upper, 1.0, 1.0);
                                 rtol = 1e-3, buffer = primary_buf)
            primary_transfer_matrix[i_primary, i_degraded, i_threshold] = result
        end
    end

    # Secondary: a single secondary still reaches at most W/2 = (E_primary_bin_max − threshold)/2.
    E_secondary_boundary_upper = (E_primary_bin_max - threshold) / 2
    i_max_secondary = searchsortedlast(E_left, E_secondary_boundary_upper)
    i_max_secondary = min(i_max_secondary, i_primary - 1)
    for i_secondary in 1:i_max_secondary
        E_secondary_bin_min = E_edges[i_secondary]
        E_secondary_bin_max = E_edges[i_secondary + 1]
        E_secondary_upper = min(E_secondary_bin_max, E_secondary_boundary_upper)
        if E_secondary_upper > E_secondary_bin_min
            result, _ = hcubature(secondary_integrand, (E_secondary_bin_min, 0.0, 0.0),
                                 (E_secondary_upper, 1.0, 1.0);
                                 rtol = 1e-3, buffer = secondary_buf)
            secondary_transfer_matrix[i_primary, i_secondary, i_threshold] = result
        end
    end
    return
end


# ======================================================================================== #
#          NUMERICAL-CDF DOUBLE IONIZATION (generic secondary-law fallback)                 #
# ======================================================================================== #

# Four-, eight-, and sixteen-point Gauss-Legendre rules on [-1, 1].  The outer
# primary/output-bin integrations use the four-point rule.  The cumulative-law table uses
# the eight-point rule, and the constrained self-convolution uses sixteen points because it
# may span the sharp low-energy part of a secondary law.
const _GL4_X = (-0.8611363115940526, -0.3399810435848563,
                 0.3399810435848563,  0.8611363115940526)
const _GL4_W = ( 0.3478548451374539,  0.6521451548625461,
                 0.6521451548625461,  0.3478548451374539)
const _GL8_X = (-0.9602898564975363, -0.7966664774136267,
                -0.5255324099163290, -0.1834346424956498,
                 0.1834346424956498,  0.5255324099163290,
                 0.7966664774136267,  0.9602898564975363)
const _GL8_W = ( 0.1012285362903763,  0.2223810344533745,
                 0.3137066458778873,  0.3626837833783620,
                 0.3626837833783620,  0.3137066458778873,
                 0.2223810344533745,  0.1012285362903763)
const _GL16_X = (-0.9894009349916499, -0.9445750230732326,
                 -0.8656312023878318, -0.7554044083550030,
                 -0.6178762444026438, -0.4580167776572274,
                 -0.2816035507792589, -0.0950125098376374,
                  0.0950125098376374,  0.2816035507792589,
                  0.4580167776572274,  0.6178762444026438,
                  0.7554044083550030,  0.8656312023878318,
                  0.9445750230732326,  0.9894009349916499)
const _GL16_W = ( 0.0271524594117541,  0.0622535239386479,
                  0.0951585116824928,  0.1246289712555339,
                  0.1495959888165767,  0.1691565193950025,
                  0.1826034150449236,  0.1894506104550685,
                  0.1894506104550685,  0.1826034150449236,
                  0.1691565193950025,  0.1495959888165767,
                  0.1246289712555339,  0.0951585116824928,
                  0.0622535239386479,  0.0271524594117541)

@inline function _checked_secondary_law(law, E_secondary, E_primary)
    value = law(E_secondary, E_primary)
    isfinite(value) || throw(DomainError(value, "secondary law must be finite"))
    value >= 0 || throw(DomainError(value, "secondary law must be nonnegative"))
    return value
end

@inline function _gauss(f, a, b, X, W)
    b > a || return 0.0
    midpoint = (a + b) / 2
    halfwidth = (b - a) / 2
    result = 0.0
    @inbounds for k in eachindex(X)
        result += W[k] * f(midpoint + halfwidth * X[k])
    end
    return halfwidth * result
end
_gauss4(f, a, b) = _gauss(f, a, b, _GL4_X, _GL4_W)
_gauss8(f, a, b) = _gauss(f, a, b, _GL8_X, _GL8_W)

"""
Numerical cumulative integral of a custom secondary law at one fixed primary energy.

`cumulative[k]` is the integral from zero through `edges[k]`.  The energy-grid edges are
used as knots, so the importance map never needs to resolve structure finer than an output
energy bin.  The bin masses themselves are evaluated with eight-point Gauss-Legendre
quadrature rather than assuming a particular analytical law.
"""
struct NumericalSecondaryCDF{F}
    edges::Vector{Float64}
    cumulative::Vector{Float64}
    E_primary::Float64
    law::F
end

"""
Build a `NumericalSecondaryCDF` for one fixed primary energy, using the energy-grid edges (up
to `E_max`) as CDF knots and eight-point Gauss-Legendre quadrature for each knot-interval mass.

Named separately from the struct's implicit field constructor (which has the same arity) so
that call sites unambiguously request the built-from-scratch cumulative table.
"""
function build_secondary_cdf(E_edges, E_max, E_primary, law)
    edges = Float64[0.0]
    for edge in E_edges
        edge <= 0 && continue
        edge >= E_max && break
        push!(edges, edge)
    end
    E_max > 0 && push!(edges, Float64(E_max))

    cumulative = zeros(length(edges))
    for i in 1:(length(edges) - 1)
        mass = _gauss8(edges[i], edges[i + 1]) do E_secondary
            _checked_secondary_law(law, E_secondary, E_primary)
        end
        cumulative[i + 1] = cumulative[i] + mass
    end
    return NumericalSecondaryCDF(edges, cumulative, Float64(E_primary), law)
end

# Accurate cumulative-law value.  Completed knot intervals come from the prefix table; only
# the final partial interval is integrated here.
@inline function cumulative_law(cdf::NumericalSecondaryCDF, energy)
    energy <= 0 && return 0.0
    energy >= cdf.edges[end] && return cdf.cumulative[end]
    i = searchsortedlast(cdf.edges, energy)
    prefix = cdf.cumulative[i]
    partial = _gauss8(cdf.edges[i], energy) do E_secondary
        _checked_secondary_law(cdf.law, E_secondary, cdf.E_primary)
    end
    return prefix + partial
end

# Piecewise-linear cumulative coordinate used solely as an importance map.  Returning its
# local slope makes the subsequent change of variables exact for this map: the physical law
# is still evaluated explicitly and divided by dC/dE.
@inline function cumulative_map(cdf::NumericalSecondaryCDF, energy)
    energy <= 0 && return 0.0
    energy >= cdf.edges[end] && return cdf.cumulative[end]
    i = searchsortedlast(cdf.edges, energy)
    fraction = (energy - cdf.edges[i]) / (cdf.edges[i + 1] - cdf.edges[i])
    return cdf.cumulative[i] +
           fraction * (cdf.cumulative[i + 1] - cdf.cumulative[i])
end

@inline function inverse_cumulative_map(cdf::NumericalSecondaryCDF, cumulative_value)
    cumulative_value <= 0 && return (0.0, Inf)
    cumulative_value >= cdf.cumulative[end] && return (cdf.edges[end], Inf)

    i = searchsortedlast(cdf.cumulative, cumulative_value)
    i = min(i, length(cdf.cumulative) - 1)
    # Flat cumulative intervals represent zero-mass bins.  Interior Gauss nodes should not
    # land there, but skip them defensively for custom laws with compact/disjoint support.
    while i < length(cdf.cumulative) && cdf.cumulative[i + 1] <= cumulative_value
        i += 1
    end
    i >= length(cdf.cumulative) && return (cdf.edges[end], Inf)

    mass = cdf.cumulative[i + 1] - cdf.cumulative[i]
    mass > 0 || return (cdf.edges[i], Inf)
    width = cdf.edges[i + 1] - cdf.edges[i]
    fraction = (cumulative_value - cdf.cumulative[i]) / mass
    return cdf.edges[i] + fraction * width, mass / width
end

@inline function double_secondary_density(E_secondary, E_primary, threshold,
                                           law, cdf::NumericalSecondaryCDF)
    W = E_primary - threshold
    (0 <= E_secondary <= W / 2) || return 0.0
    partner_upper = min(W - 2 * E_secondary, (W - E_secondary) / 2)
    partner_upper > 0 || return 0.0
    return _checked_secondary_law(law, E_secondary, E_primary) *
           cumulative_law(cdf, partner_upper)
end

@inline function double_primary_density_cdf(E_degraded, E_primary, threshold,
                                             law, cdf::NumericalSecondaryCDF)
    W = E_primary - threshold
    (W / 3 <= E_degraded <= W) || return 0.0
    secondary_sum = W - E_degraded
    secondary_sum > 0 || return 0.0

    E_lower = max(0.0, secondary_sum - E_degraded)
    E_mid = secondary_sum / 2
    C_lower = cumulative_map(cdf, E_lower)
    C_mid = cumulative_map(cdf, E_mid)
    C_width = C_mid - C_lower
    C_width > 0 || return 0.0

    # The convolution interval is symmetric about secondary_sum/2.  Integrate its lower
    # half in the numerical cumulative coordinate and double it.  Mapping y from [-1,1]
    # to [C_lower,C_mid] and applying the symmetry factor leaves the prefactor C_width.
    result = 0.0
    @inbounds for k in eachindex(_GL16_X)
        C = (C_lower + C_mid) / 2 + (C_width / 2) * _GL16_X[k]
        E_secondary, map_density = inverse_cumulative_map(cdf, C)
        isfinite(map_density) || continue
        physical_density = _checked_secondary_law(law, E_secondary, E_primary)
        partner_density = _checked_secondary_law(law,
                                                  secondary_sum - E_secondary,
                                                  E_primary)
        result += _GL16_W[k] * physical_density * partner_density / map_density
    end
    return C_width * result
end

"""
Double-ionization row builder based on numerical marginalization.

The primary-bin integral is a four-point Gauss-Legendre sum.  At each primary-energy node,
a numerical CDF of the arbitrary secondary law removes the partner dimension from the
secondary marginal and supplies an importance coordinate for the degraded-primary
self-convolution.  Output-bin integrations are also four-point Gauss-Legendre sums.
"""
function fill_double_ionization_bin_cdf!(primary_transfer_matrix,
                                         secondary_transfer_matrix,
                                         E_edges, E_left, threshold,
                                         i_primary, i_threshold, law)
    E_primary_min = E_edges[i_primary]
    E_primary_max = E_edges[i_primary + 1]
    E_primary_mid = (E_primary_min + E_primary_max) / 2
    E_primary_halfwidth = (E_primary_max - E_primary_min) / 2

    # Accuracy limit of this fixed-order primary-bin rule: the output-bin clamps
    # (max(edge, W/3), min(edge, W) for the degraded primary) introduce kinks in the
    # integrand as E_primary sweeps the bin, which the four-point Gauss rule does not resolve.
    # The least accurate rows are therefore the near-threshold ones, where W varies strongly
    # across the bin. Measured worst per-row degraded-primary sum error is ~1.3e-3 vs the
    # rtol=1e-3 hcubature reference on the standard 3 keV grid — physically negligible, since
    # those rows carry little weight.
    @inbounds for k_primary in eachindex(_GL4_X)
        E_primary = E_primary_mid + E_primary_halfwidth * _GL4_X[k_primary]
        primary_weight = E_primary_halfwidth * _GL4_W[k_primary]
        W = E_primary - threshold
        W > 0 || continue
        cdf = build_secondary_cdf(E_edges, W / 2, E_primary, law)

        # Degraded primary marginal: d in [W/3, W].
        i_min_degraded = max(1, searchsortedlast(E_left, W / 3))
        i_max_degraded = min(i_primary, searchsortedlast(E_left, W))
        for i_degraded in i_min_degraded:i_max_degraded
            lower = max(E_edges[i_degraded], W / 3)
            upper = min(E_edges[i_degraded + 1], W)
            upper > lower || continue
            value = _gauss4(lower, upper) do E_degraded
                double_primary_density_cdf(E_degraded, E_primary, threshold, law, cdf)
            end
            primary_transfer_matrix[i_primary, i_degraded, i_threshold] +=
                primary_weight * value
        end

        # Per-secondary marginal: s in [0, W/2].  Energies below the grid floor remain
        # unbinned exactly as in the HCubature implementation, but they are retained in the
        # partner CDF and therefore still influence the on-grid marginal.
        i_max_secondary = min(i_primary - 1, searchsortedlast(E_left, W / 2))
        for i_secondary in 1:i_max_secondary
            lower = E_edges[i_secondary]
            upper = min(E_edges[i_secondary + 1], W / 2)
            upper > lower || continue
            value = _gauss4(lower, upper) do E_secondary
                double_secondary_density(E_secondary, E_primary, threshold, law, cdf)
            end
            secondary_transfer_matrix[i_primary, i_secondary, i_threshold] +=
                primary_weight * value
        end
    end
    return
end


# Load the secondary electron distribution, for a given initial primary energy index
# and ionization threshold.
function secondary_spectrum(cache::SpeciesCascadingCache, i_primary::Integer,
                            E_ionization_threshold)

    i_threshold = findmin(x -> abs(x - E_ionization_threshold), cache.ionization_thresholds)[2]
    return @view(cache.secondary_transfer_matrix[i_primary, :, i_threshold])
end

function secondary_spectrum(cache::SpeciesCascadingCache, E_primary_energy,
                            E_ionization_threshold)

    i_primary = searchsortedlast(cache.E_edges, E_primary_energy) - 1
    i_primary = clamp(i_primary, 1, size(cache.secondary_transfer_matrix, 1))
    return secondary_spectrum(cache, i_primary, E_ionization_threshold)
end


# Load the degraded primary electron distribution, for a given initial primary energy index
# and ionization threshold.
function primary_spectrum(cache::SpeciesCascadingCache, i_primary::Integer,
                          E_ionization_threshold)

    i_threshold = findmin(x -> abs(x - E_ionization_threshold), cache.ionization_thresholds)[2]
    return @view(cache.primary_transfer_matrix[i_primary, :, i_threshold])
end

function primary_spectrum(cache::SpeciesCascadingCache, E_primary_energy,
                          E_ionization_threshold)

    i_primary = searchsortedlast(cache.E_edges, E_primary_energy) - 1
    i_primary = clamp(i_primary, 1, size(cache.primary_transfer_matrix, 1))
    return primary_spectrum(cache, i_primary, E_ionization_threshold)
end
