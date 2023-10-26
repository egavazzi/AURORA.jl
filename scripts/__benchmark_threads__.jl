## Multi-threaded part of add_ionization_collisions!
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
# using LoopVectorization
# function test_turbo(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
#     for n in 1:1 # number of E_levels in total
#         @turbo thread=true  for iI in 1:(iE - 1)
#             for j in axes(Q, 2)
#                 for k in axes(Q, 1)
#                     Q[k, j, iI] += Ionization[k, j] * secondary_e_spectra[iI] +
#                                             Ionizing[k, j] * primary_e_spectra[iI]
#                 end
#             end
#         end
#     end
# end
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
@btime test_serial(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
# @btime test_turbo(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
@btime test_threads(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
@btime test_threads_chunks(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE, 6)
@btime test_spawn(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE)
@btime test_polyester(Ionization, Ionizing, secondary_e_spectra, primary_e_spectra, Q, iE, 6);
