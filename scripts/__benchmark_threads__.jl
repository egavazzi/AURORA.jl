## Multi-threaded part of add_ionization_collisions!
using BenchmarkTools
n_E = 680;
n_z = 5742;
n_t = 51;
# n_t = 401;
iE = 680;
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


















#################################################################################
#                           several steps at once                               #
#################################################################################

## investigating if adding all ionization to Q in one go is faster or not
using BenchmarkTools
# n_E = 680;
n_E = 340;
n_z = 5742;
# n_t = 51;
n_t = 401;
iE = 340;
Q = zeros(n_z, n_t, n_E);

Ionization = rand(size(Q, 1), size(Q, 2));
Ionizing = rand(size(Q, 1), size(Q, 2));
secondary_e_spectra = rand(n_E);
primary_e_spectra = rand(n_E);

# Ionization2 = rand(size(Q, 1), size(Q, 2));
# Ionizing2 = rand(size(Q, 1), size(Q, 2));
# secondary_e_spectra2 = rand(n_E);
# primary_e_spectra2 = rand(n_E);

# Ionization3 = rand(size(Q, 1), size(Q, 2));
# Ionizing3 = rand(size(Q, 1), size(Q, 2));
# secondary_e_spectra3 = rand(n_E);
# primary_e_spectra3 = rand(n_E);

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
        @batch minbatch=nbatch for iI in 1:(iE - 1)
        # Threads.@threads for i in 1:Nthreads
        #     for iI in splitter(iE - 1, Nthreads, i)
                @view(Q[:, :, iI]) .+= Ionization .* secondary_e_spectra[iI] .+
                                        Ionizing .* primary_e_spectra[iI]
            # end
        end
    end
end
function test_bandwidth3(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra,
    Ionization2, Ionizing2, secondary_e_spectra2, primary_e_spectra2,
    Ionization3, Ionizing3, secondary_e_spectra3, primary_e_spectra3,
    Q, iE, Nthreads)
    for n in 1:1 # number of E_levels in total
        # nbatch = Int(floor(iE / Nthreads)) - 1 # will run on Nthreads threads, 6 seems to be optimal on my machine
        # nbatch < 1 ? nbatch = 1 : nothing
        # @batch minbatch=nbatch for iI in 1:(iE - 1)
        Threads.@threads for i in 1:Nthreads
            for iI in splitter(iE - 1, Nthreads, i)
            @view(Q[:, :, iI]) .+= Ionization .* secondary_e_spectra[iI] .+
                                    Ionizing .* primary_e_spectra[iI] .+
                                    Ionization2 .* secondary_e_spectra2[iI] .+
                                    Ionizing2 .* primary_e_spectra2[iI] .+
                                    Ionization3 .* secondary_e_spectra3[iI] .+
                                    Ionizing3 .* primary_e_spectra3[iI]
            end
        end
    end
end
function test_bandwidth15(Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector,
    Q, iE, Nthreads)
    for n in 1:1 # number of E_levels in total
        nbatch = Int(floor(iE / Nthreads)) - 1 # will run on Nthreads threads, 6 seems to be optimal on my machine
        nbatch < 1 ? nbatch = 1 : nothing
        @batch minbatch=nbatch for iI in 1:(iE - 1)
        # Threads.@threads for i in 1:Nthreads
        #     for iI in splitter(iE - 1, Nthreads, i)
                @view(Q[:, :, iI]) .+= Ionization_matrix[1] .* secondary_vector[1][iI] .+
                                    Ionizing_matrix[1] .* primary_vector[1][iI] .+
                                    Ionization_matrix[2] .* secondary_vector[2][iI] .+
                                    Ionizing_matrix[2] .* primary_vector[2][iI] .+
                                    Ionization_matrix[3] .* secondary_vector[3][iI] .+
                                    Ionizing_matrix[3] .* primary_vector[3][iI] .+
                                    Ionization_matrix[4] .* secondary_vector[4][iI] .+
                                    Ionizing_matrix[4] .* primary_vector[4][iI]
                                    # Ionization_matrix[5] .* secondary_vector[5][iI] .+
                                    # Ionizing_matrix[5] .* primary_vector[5][iI]
                                    # Ionization_matrix[6] .* secondary_vector[6][iI] .+
                                    # Ionizing_matrix[6] .* primary_vector[6][iI] .+
                                    # Ionization_matrix[7] .* secondary_vector[7][iI] .+
                                    # Ionizing_matrix[7] .* primary_vector[7][iI] .+
                                    # Ionization_matrix[8] .* secondary_vector[8][iI] .+
                                    # Ionizing_matrix[8] .* primary_vector[8][iI] .+
                                    # Ionization_matrix[9] .* secondary_vector[9][iI] .+
                                    # Ionizing_matrix[9] .* primary_vector[9][iI] .+
                                    # Ionization_matrix[10] .* secondary_vector[10][iI] .+
                                    # Ionizing_matrix[10] .* primary_vector[10][iI] .+
                                    # Ionization_matrix[11] .* secondary_vector[11][iI] .+
                                    # Ionizing_matrix[11] .* primary_vector[11][iI] .+
                                    # Ionization_matrix[12] .* secondary_vector[12][iI] .+
                                    # Ionizing_matrix[12] .* primary_vector[12][iI] .+
                                    # Ionization_matrix[13] .* secondary_vector[13][iI] .+
                                    # Ionizing_matrix[13] .* primary_vector[13][iI] .+
                                    # Ionization_matrix[14] .* secondary_vector[14][iI] .+
                                    # Ionizing_matrix[14] .* primary_vector[14][iI] .+
                                    # Ionization_matrix[15] .* secondary_vector[15][iI] .+
                                    # Ionizing_matrix[15] .* primary_vector[15][iI]
            # end
        end
    end
end
using LoopVectorization
function test_bandwidth_turbo(Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector,
    Q, iE)
    for n in 1:1 # number of E_levels in total
        @turbo thread=20 for iI in 1:(iE - 1)
        # @tturbo for iI in 1:(iE - 1)
            for j in axes(Q, 2)
                for k in axes(Q, 1)
                    Q[k, j, iI] += Ionization_matrix[1][k, j] * secondary_vector[1][iI] +
                                        Ionizing_matrix[1][k, j] * primary_vector[1][iI] +
                                        Ionization_matrix[2][k, j] * secondary_vector[2][iI] +
                                        Ionizing_matrix[2][k, j] * primary_vector[2][iI] +
                                        Ionization_matrix[3][k, j] * secondary_vector[3][iI] +
                                        Ionizing_matrix[3][k, j] * primary_vector[3][iI] +
                                        Ionization_matrix[4][k, j] * secondary_vector[4][iI] +
                                        Ionizing_matrix[4][k, j] * primary_vector[4][iI] +
                                        Ionization_matrix[5][k, j] * secondary_vector[5][iI]
                                        # Ionizing_matrix[5][k, j] * primary_vector[5][iI] +
                                        # Ionization_matrix[6][k, j] * secondary_vector[6][iI] +
                                        # Ionizing_matrix[6][k, j] * primary_vector[6][iI] +
                                        # Ionization_matrix[7][k, j] * secondary_vector[7][iI] +
                                        # Ionizing_matrix[7][k, j] * primary_vector[7][iI] +
                                        # Ionization_matrix[8][k, j] * secondary_vector[8][iI] +
                                        # Ionizing_matrix[8][k, j] * primary_vector[8][iI] +
                                        # Ionization_matrix[9][k, j] * secondary_vector[9][iI] +
                                        # Ionizing_matrix[9][k, j] * primary_vector[9][iI] +
                                        # Ionization_matrix[10][k, j] * secondary_vector[10][iI] +
                                        # Ionizing_matrix[10][k, j] * primary_vector[10][iI]
                end
            end
        end
    end
end
##
@benchmark test_bandwidth1(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE, 20) seconds=15
# @btime test_bandwidth3(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra,
                                    # Ionization2, Ionizing2, secondary_e_spectra2, primary_e_spectra2,
                                    # Ionization3, Ionizing3, secondary_e_spectra3, primary_e_spectra3,
                                    # Q, iE, 40)
@btime test_bandwidth15(Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector,
                                    Q, iE, 10) seconds=15
@benchmark test_bandwidth_turbo(Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector,
                                    Q, iE) seconds=15
