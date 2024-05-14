#################################################################################
#               Multi-threaded part of add_ionization_collisions!               #
#################################################################################

using BenchmarkTools
n_E = 337;
n_z = 5742;
n_t = 51;
# n_t = 401;
iE = 337;
Q = zeros(n_z, n_t, n_E);
Ionization = rand(size(Q, 1), size(Q, 2));
Ionizing = rand(size(Q, 1), size(Q, 2));
secondary_e_spectra = rand(n_E);
primary_e_spectra = rand(n_E);
##
function splitter(n, nchunks, ichunk)
    n_per_chunk = ceil(Int, n / nchunks) # only works for multiples
    first = (ichunk - 1) * n_per_chunk + 1
    if ichunk * n_per_chunk <= n
        last = ichunk * n_per_chunk
    else
        last = n
    end
    return first:last
end
function test_threads(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
    for n in 1:1 # number of E_levels in total
        Threads.@threads for iI in 1:(iE - 1)
            @view(Q[:, :, iI]) .+= Ionization .* secondary_e_spectra[iI] .+
                                    Ionizing .* primary_e_spectra[iI]
        end
    end
end
function test_threads_chunks(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE, Nthreads)
    for n in 1:1 # number of E_levels in total
        Threads.@threads for i in 1:Nthreads
            for iI in splitter(iE - 1, Nthreads, i)
                @views(Q[:, :, iI]) .+=  Ionization  .* secondary_e_spectra[iI] .+
                                                Ionizing .* primary_e_spectra[iI]
            end
        end
    end
end
function test_serial(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
    for n in 1:1 # number of E_levels in total
        for iI in 1:(iE - 1)
            @view(Q[:, :, iI]) .+= sin.(Ionization) .* secondary_e_spectra[iI] .+
                                    Ionizing .* primary_e_spectra[iI]
        end
    end
end
using LoopVectorization
function test_turbo(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
    for n in 1:1 # number of E_levels in total
        @turbo thread=4  for iI in 1:(iE - 1)
            for j in axes(Q, 2)
                for k in axes(Q, 1)
                    Q[k, j, iI] += Ionization[k, j] * secondary_e_spectra[iI] +
                                            Ionizing[k, j] * primary_e_spectra[iI]
                end
            end
        end
    end
end
function test_spawn(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
    for n in 1:1 # number of E_levels in total
        @sync for iI in 1:(iE - 1)
            Threads.@spawn @view(Q[:, :, iI]) .+= sin.(Ionization) .* secondary_e_spectra[iI] .+
                                    Ionizing .* primary_e_spectra[iI]
        end
    end
end
using Polyester
function test_polyester(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE, Nthreads)
    for n in 1:1 # number of E_levels in total
        nbatch = Int(floor(iE / Nthreads)) - 1 # will run on Nthreads threads, 6 seems to be optimal on my machine
        nbatch < 1 ? nbatch = 1 : nothing
        @batch minbatch=nbatch for iI in 1:(iE - 1)
            @view(Q[:, :, iI]) .+= Ionization .* secondary_e_spectra[iI] .+
                                    Ionizing .* primary_e_spectra[iI]
        end
    end
end
##
@btime test_serial(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
@btime test_turbo(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
@btime test_threads(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
@btime test_threads_chunks(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, 100, 20)
@btime test_spawn(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
@btime test_polyester(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, 100, 20)



Q = zeros(n_z, n_t, n_E);

test_turbo(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)

Q_save = copy(Q);

Q = zeros(n_z, n_t, n_E);

test_polyester(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE, 20)

Q .== Q_save
Q ≈ Q_save

Q[:, :, 1]
Q_save[:, :, 1]









#################################################################################
#                           several steps at once                               #
#################################################################################

## investigating if adding all ionization to Q in one go is faster or not
using BenchmarkTools
n_E = 680;
n_z = 5742;
n_t = 51;
iE = 680;
Q = zeros(n_z, n_t, n_E);

Ionization = rand(size(Q, 1), size(Q, 2));
Ionizing = rand(size(Q, 1), size(Q, 2));
secondary_e_spectra = rand(n_E);
primary_e_spectra = rand(n_E);

Ionization_matrix = [rand(size(Q, 1), size(Q, 2)) for _ in 1:15];
Ionizing_matrix = [rand(size(Q, 1), size(Q, 2)) for _ in 1:15];
secondary_vector = [rand(n_E) for _ in 1:15];
primary_vector = [rand(n_E) for _ in 1:15];
##
using Polyester
function splitter(n, nchunks, ichunk)
    n_per_chunk = ceil(Int, n / nchunks) # only works for multiples
    first = (ichunk - 1) * n_per_chunk + 1
    if ichunk * n_per_chunk <= n
        last = ichunk * n_per_chunk
    else
        last = n
    end
    return first:last
end

function test_bandwidth1(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE, Nthreads)
    for n in 1:5 # number of E_levels in total
        nbatch = Int(floor(iE / Nthreads)) - 1 # will run on Nthreads threads, 6 seems to be optimal on my machine
        nbatch < 1 ? nbatch = 1 : nothing
        @inbounds @batch minbatch=nbatch for iI in 1:(iE - 1)
        # Threads.@threads for i in 1:Nthreads
        #     for iI in splitter(iE - 1, Nthreads, i)
                @view(Q[:, :, iI]) .+= Ionization .* secondary_e_spectra[iI] .+
                                        Ionizing .* primary_e_spectra[iI]
            # end
        end
    end
end
using LoopVectorization
function test_bandwidth_turbo(Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector,
    Q, iE)
    for n in 1:1 # number of E_levels in total
        # @turbo inline=false thread=10 for iI in 1:(iE - 1)
        @tturbo inline=false for iI in 1:(iE - 1)
        # @inbounds Threads.@threads for iI in 1:(iE - 1)
            for j in axes(Q, 2)
                for k in axes(Q, 1)
                    Q[k, j, iI] += Ionization_matrix[1][k, j] * secondary_vector[1][iI] +
                                        Ionizing_matrix[1][k, j] * primary_vector[1][iI] +
                                        Ionization_matrix[2][k, j] * secondary_vector[2][iI] +
                                        Ionizing_matrix[2][k, j] * primary_vector[2][iI] +
                                        Ionization_matrix[3][k, j] * secondary_vector[3][iI] +
                                        Ionizing_matrix[3][k, j] * primary_vector[3][iI]  +
                                        Ionization_matrix[4][k, j] * secondary_vector[4][iI] +
                                        Ionizing_matrix[4][k, j] * primary_vector[4][iI] +
                                        Ionization_matrix[5][k, j] * secondary_vector[5][iI] +
                                        Ionizing_matrix[5][k, j] * primary_vector[5][iI]
                                        # Ionization_matrix[6][k, j] * secondary_vector[6][iI] +
                                        # Ionizing_matrix[6][k, j] * primary_vector[6][iI]
                                        # Ionization_matrix[7][k, j] * secondary_vector[7][iI] +
                                        # Ionizing_matrix[7][k, j] * primary_vector[7][iI]
                                        # Ionization_matrix[8][k, j] * secondary_vector[8][iI] +
                                        # Ionizing_matrix[8][k, j] * primary_vector[8][iI]
                                        # Ionization_matrix[9][k, j] * secondary_vector[9][iI] +
                                        # # Ionizing_matrix[9][k, j] * primary_vector[9][iI] +
                end
            end
        end
    end
end
##
@benchmark test_bandwidth1(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE, 20) seconds=5
b = @benchmark test_bandwidth_turbo($Ionization_matrix, $Ionizing_matrix, $secondary_vector, $primary_vector,
                                    $Q, $iE) seconds=5

# max size 5742*101*680 ~ 400 000 000 values (3 GiB)
# max size 5742*101 ~ 575 000 values (4.4 MiB slices)











#################################################################################
#                    automatic investigation of ideal size of Q                 #
#################################################################################
##
function test_size_Q(nt_to_benchmark)
    benchmark_results = zeros(length(nt_to_benchmark), 2)

    for (i, n_t) in enumerate(nt_to_benchmark)
        n_E = 680;
        iE = 680;
        n_z = 5742;
        # n_t = i * 50 + 1;
        Q_local = zeros(n_z, n_t, n_E);
        Ionization_matrix_local = [rand(size(Q_local, 1), size(Q_local, 2)) for _ in 1:15];
        Ionizing_matrix_local = [rand(size(Q_local, 1), size(Q_local, 2)) for _ in 1:15];
        secondary_vector_local = [rand(n_E) for _ in 1:15];
        primary_vector_local = [rand(n_E) for _ in 1:15];

        # println(size(Q_local))
        b = @benchmark test_bandwidth_turbo($Ionization_matrix_local, $Ionizing_matrix_local, $secondary_vector_local, $primary_vector_local,
                $Q_local, $iE) samples=200 seconds=20

        display(b)

        benchmark_results[i, 1] = mean(b).time * 1e-9   # convert from ns to s
        benchmark_results[i, 2] = n_z * n_t
    end

    return benchmark_results
end

X = test_size_Q([11, 21, 41, 51, 101, 201, 401, 501, 601, 801, 1001]);

##
using CairoMakie
f = Figure()
ax = Axis(f[1, 1], xscale = identity)
scatter!(ax, X[:, 2] ./ 5742, X[:, 1])
f
lines!(ax, [51, 1001], [X[4, 1], X[4, 1] * (1001 / 51)])
lines!(ax, [41, 1001], [X[3, 1], X[3, 1] * (1001 / 41)])
f
##

using AURORA
h_atm = AURORA.make_altitude_grid(500)
a, b = CFL_criteria(range(0, 0.35, 106)[1:11], h_atm, v_of_E(3000), 128)
size(a)

















##
#################################################################################
#               Multi-threaded part of add_inelastic_collisions!               #
#################################################################################
using AURORA
using BenchmarkTools

altitude_max = 600;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 7000;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith
t_sampling = 0:0.001:0.01;   # (s) time-array over which data will be saved
n_loop = 1;                 # number of loops to run
msis_file = find_nrlmsis_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );
iri_file = find_iri_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );

h_atm, ne, Te, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup_new(altitude_max, θ_lims, E_max, msis_file, iri_file);

t, CFL_factor = CFL_criteria(t_sampling, h_atm, v_of_E(E_max), 64)

phaseN2e, phaseN2i = phase_fcn_N2(μ_scatterings.theta1, E);
phaseO2e, phaseO2i = phase_fcn_O2(μ_scatterings.theta1, E);
phaseOe, phaseOi = phase_fcn_O(μ_scatterings.theta1, E);
phase_fcn_neutrals = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));
cascading_neutrals = (cascading_N2, cascading_O2, cascading_O) # tuple of functions

Q  = zeros(length(h_atm) * length(μ_center), length(t), length(E));
Ie = rand(length(h_atm) * length(μ_center), length(t), length(E));

B, B2B_inelastic_neutrals = make_B(n_neutrals, σ_neutrals, E_levels_neutrals,
                                    phase_fcn_neutrals, dE, iE, μ_scatterings.Pmu2mup,
                                    μ_scatterings.BeamWeight_relative, μ_scatterings.theta1);

##
n = n_neutrals[3];                          # Neutral density
σ = σ_neutrals[3];                          # Array with collision cross sections
E_levels = E_levels_neutrals[3];            # Array with collision enery levels and number of secondary e-
B2B_inelastic = B2B_inelastic_neutrals[3];  # Array with the probablities of scattering from beam to beam
cascading = cascading_neutrals[3];          # Cascading function for the current i-th species

##
function fff!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE)
    Ie_degraded = Matrix{Float64}(undef, size(Ie, 1), size(Ie, 2))

    # Multiply each element of B2B with n (density vector) and resize to get a matrix that can
    # be multiplied with Ie
    AB2B =  AURORA.make_big_B2B_matrix(B2B_inelastic, n, h_atm)
    Ie_scatter = AB2B * @view(Ie[:, :, iE])

    # Loop over the energy levels of the collisions with the i-th neutral species
    for i_level in 2:size(E_levels, 1)
        if E_levels[i_level, 2] <= 0  # these collisions should not produce secondary e-
            # The flux of e- degraded from energy bin [E[iE], E[iE] + dE[iE]] to any lower energy
            # bin by excitation of the E_levels[i_level] state of the current species.
            # The second factor corrects for the case where the energy loss is maller than the width
            # in energy of the energy bin. That is, when dE[iE] > E_levels[i_level,1], only the
            # fraction E_levels[i_level,1] / dE[iE] is lost from the energy bin [E[iE], E[iE] + dE[iE]].

            Ie_degraded .= (σ[i_level, iE] * min(1, E_levels[i_level, 1] / dE[iE])) .* Ie_scatter

            # Find the energy bins where the e- in the current energy bin will degrade when losing
            # E_levels[i_level, 1] eV
            i_degrade = intersect(  findall(x -> x > E[iE] - E_levels[i_level, 1], E + dE),     # find lowest bin
                                    findall(x -> x < E[iE] + dE[iE] - E_levels[i_level,1], E))  # find highest bin

            partition_fraction = zeros(size(i_degrade)) # initialise
            if !isempty(i_degrade) && i_degrade[1] < iE
                # Distribute the degrading e- between those bins
                partition_fraction[1] = min(1, (E[i_degrade[1]] .+ dE[i_degrade[1]] .-
                                                E[iE] .+ E_levels[i_level, 1]) / dE[iE])
                if length(i_degrade) > 2
                    partition_fraction[2:end-1] = min.(1, dE[i_degrade[2:end-1]] / dE[iE])
                end
                partition_fraction[end] = min(1, (E[iE] .+ dE[iE] .- E[i_degrade[end]] .-
                                                    E_levels[i_level, 1]) / dE[iE])
                if i_degrade[end] == iE
                    partition_fraction[end] = 0
                end

                # normalise
                partition_fraction = partition_fraction / sum(partition_fraction)

                # println(partition_fraction)
                # println(findall(x -> x != 0, partition_fraction))
                # println(eachindex(partition_fraction))
                # println()

                # and finally calculate the flux of degrading e-
                # idx = findall(x -> x != 0, partition_fraction)
                Threads.@threads for i_u in findall(x -> x != 0, partition_fraction)
                # Threads.@threads for i_u in eachindex(findall(x -> x != 0, partition_fraction))
                # for i_u in findall(x -> x != 0, partition_fraction)
                # for i_u in eachindex(idx)
                    @view(Q[:, :, i_degrade[i_u]]) .+=  max.(0, Ie_degraded) .* partition_fraction[i_u]
                end
            end
        end
    end
end

fill!(Q, 0.0);

@time fff!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE)
@btime fff!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE)
@profview fff!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE)

Q_save = copy(Q);

using LoopVectorization
function ggg!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE)
    Ie_degraded = Matrix{Float64}(undef, size(Ie, 1), size(Ie, 2))

    # Multiply each element of B2B with n (density vector) and resize to get a matrix that can
    # be multiplied with Ie
    AB2B =  AURORA.make_big_B2B_matrix(B2B_inelastic, n, h_atm)
    Ie_scatter = AB2B * @view(Ie[:, :, iE])

    # Loop over the energy levels of the collisions with the i-th neutral species
    for i_level in 2:size(E_levels, 1)
        if E_levels[i_level, 2] <= 0  # these collisions should not produce secondary e-
            # The flux of e- degraded from energy bin [E[iE], E[iE] + dE[iE]] to any lower energy
            # bin by excitation of the E_levels[i_level] state of the current species.
            # The second factor corrects for the case where the energy loss is maller than the width
            # in energy of the energy bin. That is, when dE[iE] > E_levels[i_level,1], only the
            # fraction E_levels[i_level,1] / dE[iE] is lost from the energy bin [E[iE], E[iE] + dE[iE]].

            Ie_degraded .= (σ[i_level, iE] * min(1, E_levels[i_level, 1] / dE[iE])) .* Ie_scatter

            # Find the energy bins where the e- in the current energy bin will degrade when losing
            # E_levels[i_level, 1] eV
            i_degrade = intersect(  findall(x -> x > E[iE] - E_levels[i_level, 1], E + dE),     # find lowest bin
                                    findall(x -> x < E[iE] + dE[iE] - E_levels[i_level,1], E))  # find highest bin

            partition_fraction = zeros(size(i_degrade)) # initialise
            if !isempty(i_degrade) && i_degrade[1] < iE
                # Distribute the degrading e- between those bins
                partition_fraction[1] = min(1, (E[i_degrade[1]] .+ dE[i_degrade[1]] .-
                                                E[iE] .+ E_levels[i_level, 1]) / dE[iE])
                if length(i_degrade) > 2
                    partition_fraction[2:end-1] = min.(1, dE[i_degrade[2:end-1]] / dE[iE])
                end
                partition_fraction[end] = min(1, (E[iE] .+ dE[iE] .- E[i_degrade[end]] .-
                                                    E_levels[i_level, 1]) / dE[iE])
                if i_degrade[end] == iE
                    partition_fraction[end] = 0
                end

                # normalise
                partition_fraction = partition_fraction / sum(partition_fraction)


                # and finally calculate the flux of degrading e-
                idx = findall(x -> x != 0, partition_fraction)
                # println(eachindex(idx))
                # println(eachindex(findall(x -> x != 0, partition_fraction)))
                # for i in idx
                #     println(i)
                # end
                # println(findall(x -> x != 0, partition_fraction))
                # Threads.@threads for i_u in findall(x -> x != 0, partition_fraction)
                # @turbo inline=false thread=1 for i_u in eachindex(idx)
                @turbo inline=false thread=true for i_u in eachindex(findall(x -> x != 0, partition_fraction))
                    for j in axes(Q, 2)
                        for k in axes(Q, 1)
                            Q[k, j, i_degrade[i_u]] +=  max(0, Ie_degraded[k, j]) * partition_fraction[i_u]
                            # Q[k, j, i_degrade[i_u]] +=  Ie_degraded[k, j] * partition_fraction[i_u]
                        end
                    end
                end
            end
        end
    end
end

@time ggg!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE)
@benchmark ggg!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE) seconds=10
@profview ggg!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE)

Q ≈ Q_save
