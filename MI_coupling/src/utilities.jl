using Distributed
using SharedArrays
using MAT
function load_fzvzmu_parallel(path_to_vlasov_initial_file, path_to_vlasov_simulation, index_specie, NPROCS)
    # find the fzvzmuXXXX.mat files
    dDir = readdir(joinpath(path_to_vlasov_simulation, "outp"), join=true)
    fzvzmu_files = dDir[contains.(dDir, "fzvzmu")]
    # load the initial fzvzmuXXXX.mat file
    initial_file = matopen(path_to_vlasov_initial_file)
        fzvzmustruct = read(initial_file, "fzvzmustruct")
    close(initial_file)

    # call the workers
    println("Call the workers...")
    addprocs(NPROCS - 1)
    # instantiate environment in all processes
    @everywhere @eval begin
        using Pkg; Pkg.activate("/mnt/data/etienne/Julia/AURORA");
    end
    # load the dependencies in all processes
    @everywhere @eval using MAT, SharedArrays

    # initialise fzvzmu as a shared array over all the processes, [n_files x Nvz x Nmu x Nz]
    Nvz = size(fzvzmustruct["f"][1], 1)
    Nmu = size(fzvzmustruct["f"][1], 2)
    Nz  = size(fzvzmustruct["f"][1], 3)
    fzvzmu = SharedArray{Float64}(length(fzvzmu_files) + 1, Nvz, Nmu, Nz)
    # save the initial fzvzmu into the shared array
    fzvzmu[1, :, :, :] = fzvzmustruct["f"][index_specie];
    # now let's goooo
    @sync @distributed for i in eachindex(fzvzmu_files)
        @async begin
            file = matopen(fzvzmu_files[i])
                fzvzmustruct = read(file, "fzvzmustruct")
            close(file)
            fzvzmu[i + 1, :, :, :] = fzvzmustruct["f"][index_specie]
        end
    end

    # and send the workers back to sleep...
    rmprocs(workers())

    return fzvzmu
end

function load_fzvzmu_serial(path_to_vlasov_initial_file, path_to_vlasov_simulation, index_specie, first_run=0)
    if first_run != 0
        # in that case we will load only the initial file (not sure this is relevant here...)
        fzvzmu_files = []
    else
        # find the fzvzmuXXXX.mat files
        dDir = readdir(joinpath(path_to_vlasov_simulation, "outp"), join=true)
        fzvzmu_files = dDir[contains.(dDir, "fzvzmu")]
    end
    # load the initial fzvzmuXXXX.mat file
    initial_file = matopen(path_to_vlasov_initial_file)
        fzvzmustruct = read(initial_file, "fzvzmustruct")
    close(initial_file)

    # initialise fzvzmu [n_files x Nvz x Nmu x Nz]
    Nvz = size(fzvzmustruct["f"][1], 1)
    Nmu = size(fzvzmustruct["f"][1], 2)
    Nz  = size(fzvzmustruct["f"][1], 3)
    fzvzmu = zeros(length(fzvzmu_files) + 1, Nvz, Nmu, Nz)

    # save the initial fzvzmu into the array
    fzvzmu[1, :, :, :] = fzvzmustruct["f"][index_specie]
    println("Loading fzvzmu 1/", size(fzvzmu, 1)," done")

    # and save the other fzvzmu into the array
    for i in eachindex(fzvzmu_files)
        file = matopen(fzvzmu_files[i])
            fzvzmustruct = read(file, "fzvzmustruct")
        close(file)
        fzvzmu[i + 1, :, :, :] = fzvzmustruct["f"][index_specie]
        println("Loading fzvzmu ", i + 1,"/", size(fzvzmu, 1), " done")
    end

    return fzvzmu
end

function load_fzvzmuIB_serial(path_to_vlasov_initial_file, path_to_vlasov_simulation, index_specie, first_run=0)
    if first_run != 0
        # in that case we will load only the initial file
        fzvzmu_files = []
    else
        # find the fzvzmuIBXXXX.mat files
        dDir = readdir(joinpath(path_to_vlasov_simulation, "outp"), join=true)
        fzvzmu_files = dDir[contains.(dDir, "fzvzmuIB")]
    end
    # load the initial fzvzmuIBXXXX.mat file
    initial_file = matopen(path_to_vlasov_initial_file)
        fzvzmustruct = read(initial_file, "fzvzmustruct")
    close(initial_file)

    # initialise fzvzmu [n_files x Nvz x Nmu x Nz]
    Nvz = size(fzvzmustruct["f"][1], 1)
    Nmu = size(fzvzmustruct["f"][1], 2)
    Nz  = size(fzvzmustruct["f"][1], 3)
    fzvzmu = zeros(length(fzvzmu_files) + 1, Nvz, Nmu, Nz)

    # save the initial fzvzmu into the array
    fzvzmu[1, :, :, :] = fzvzmustruct["f"][index_specie]
    println("Loading fzvzmu 1/", size(fzvzmu, 1)," done")

    # and save the other fzvzmu into the array
    for i in eachindex(fzvzmu_files)
        file = matopen(fzvzmu_files[i])
            fzvzmustruct = read(file, "fzvzmustruct")
        close(file)
        fzvzmu[i + 1, :, :, :] = fzvzmustruct["f"][index_specie]
        println("Loading fzvzmu ", i + 1,"/", size(fzvzmu, 1), " done")
    end

    return fzvzmu
end


function load_Bfield(path_to_vlasov_simulation)
    file = matopen(joinpath(path_to_vlasov_simulation, "outp", "Bfield.mat"))
        B = read(file, "B")
        dB = read(file, "dB")
        Nz = read(file, "Nz")
        z = read(file, "z")
        zcorn = read(file, "zcorn")
        dz = read(file, "dz")
        dt = read(file, "dt")
        particle = read(file, "particle") # will come out as a Dict
    close(file)

    # convert particle into a NamedTuple
    names = Tuple(Symbol.(keys(particle)))
    val = collect(values(particle))

    Nspecies = length(val[1])
    particle = Array{NamedTuple}(undef, Nspecies)
    for index_species in 1:Nspecies
        val_for_species = [] # initialise
        for i_line in eachindex(val)
            if val[i_line][index_species] isa Matrix
                # this is because Matlab save vectors as matrices, and we should not have matrices here
                push!(val_for_species, vec(val[i_line][index_species]))
            else
                push!(val_for_species, val[i_line][index_species])
            end
        end
        particle[index_species] = NamedTuple{names}(val_for_species)
    end

    return B, dB, Nz, z, zcorn, dz, dt, particle
end


## ====================================================================================== ##

"""
    findnearestindex(X, Y)

This function returns the index of the nearest element of a value Y in a **sorted** array X.

# Calling
index = findnearestindex(X, Y)
"""
function findnearestindex(X, Y)
    # X is an Array
    # Y is the value for which we want to find the index in the array X
    idx1 = searchsortedfirst(X, Y)
    idx2 = searchsortedlast(X, Y)
    # case when Y is bigger than the maximum value of X (to avoid out of bound error)
    if idx1 > length(X) || idx2 > length(X)
        if idx1 < idx2
            return idx1
        else
            return idx2
        end

    # case when Y is smaller than the minimum value of X (to avoid out of bound error)
    elseif idx1 < 1 || idx2 < 1
        if idx1 > idx2
            return idx1
        else
            return idx2
        end
        print(b)
    # normal case when Y is contained inside the values of X
    else
        if abs(X[idx1] - Y) < abs(X[idx2] - Y)
            return idx1
        else
            return idx2
        end
    end
end
