function read_input(path_to_my_data)
    mydata = readlines(path_to_my_data)

    ## General parameters
    # defaults
    dt = 1.0
    Niter = 1
    dump_period_fields = 1
    fields_per_file = 1
    dump_period_distr = 1
    dump_period_distr_IonoBoundary = 1
    dump_period_distr_1v = 10000000
    dump_start = 1
    shift_test_period = 1
    resistance = 0.0
    Nz = 10
    zmin = 0.0
    zmax = 1.0
    Nspecies = 0
    const_a = 10.0
    BC_Poisson=1
    voltage = 0.0
    initialiser = 1
    voltage_init = 0.0
    E0 = 0.0
    startfromdumpfile = false
    dump_period_dump = 1000
    exit_after_dump = false
    transffilename = "transfb2.dat"
    voltagefilename = "No file name here, sorry."

    # then read the file
    finish = findfirst(contains.(mydata, "%END"))
    for i_lines in eachindex(mydata)[1:finish]
        i = findfirst("=", mydata[i_lines])
        j = findfirst(";", mydata[i_lines])
        # check if there are "=" and ";" signs on the current line, and that the line is not a comment
        if !isnothing(i) && !isnothing(j) && (mydata[i_lines][1] != '%')
            i = i[1] # convert to integer
            j = j[1] # convert to integer
            if contains(mydata[i_lines][1:i], "dt")
                dt = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "Niter")
                Niter = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "dump_period_fields")
                dump_period_fields = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "fields_per_file")
                fields_per_file = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "dump_period_distr")
                dump_period_distr = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "dump_period_distr_IonoBoundary")
                dump_period_distr_IonoBoundary = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "dump_period_distr_1v")
                dump_period_distr_1v = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "dump_start")
                dump_start = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "shift_test_period")
                shift_test_period = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "resistance")
                resistance = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "Nz")
                Nz = tryparse(Int64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "zmin")
                zmin = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "zmax")
                zmax = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "Nspecies")
                Nspecies = tryparse(Int64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "const_a")
                const_a = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "BC_Poisson")
                BC_Poisson = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "voltage ")
                voltage = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "initialiser")
                initialiser = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "voltage_init")
                voltage_init = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "E0")
                E0 = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "startfromdumpfile")
                if contains(mydata[i_lines][1:i], "'yes'")
                    startfromdumpfile = true
                end
            elseif contains(mydata[i_lines][1:i], "dump_period_dump")
                dump_period_dump = tryparse(Float64, mydata[i_lines][i+1:j-1])
            elseif contains(mydata[i_lines][1:i], "exit_after_dump")
                if contains(mydata[i_lines][1:i], "'yes'")
                    exit_after_dump = true
                end        
            elseif contains(mydata[i_lines][1:i], "transffilename")
                transffilename = mydata[i_lines][i+1:j-1]
            elseif contains(mydata[i_lines][1:i], "voltagefilename")
                voltagefilename = mydata[i_lines][i+1:j-1]
            else
                println("Input quantity ", mydata[i_lines][1:i], " is unknown.")
            end
        end
    end

    ## Species specific parameters
    # defaults
    Nvz = 2
    vzmin = 1
    vzmax = 1
    Nmu = 2
    mumin = 0
    mumax = 1
    muexp = 1
    mass = 1.674927211e-27
    charge = 0
    n0 = 1e-99
    vz0 = 0
    kTz = 1
    kTp = 1
    n0L = 0
    vz0L = 0
    kTzL = 1
    kTpL = 1
    n0R = 0
    vz0R = 0
    kTzR = 1
    kTpR = 1

    # then read the file
    particle = Array{NamedTuple}(undef, Int64(Nspecies))
    for ii in 1:Nspecies
        start = findnext(contains.(mydata, "%SPEC"), finish)
        finish = findnext(contains.(mydata, "%END"), start)
        for i_lines in eachindex(mydata)[start:finish]
            i = findfirst("=", mydata[i_lines])
            j = findfirst(";", mydata[i_lines])
            # check if there are "=" and ";" signs on the current line, and that the line is not a comment
            if !isnothing(i) && !isnothing(j) && (mydata[i_lines][1] != '%')
                i = i[1] # convert to integer
                j = j[1] # convert to integer
                if contains(mydata[i_lines][1:i], "Nvz")
                    Nvz = tryparse(Int64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "vzmin")
                    vzmin = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "vzmax")
                    vzmax = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "Nmu")
                    Nmu = tryparse(Int64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "mumin")
                    mumin = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "mumax")
                    mumax = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "muexp")
                    muexp = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "mass")
                    mass = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "charge")
                    charge = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "n0")
                    n0 = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "vz0")
                    vz0 = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "kTz")
                    kTz = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "kTp")
                    kTp = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "n0L")
                    n0L = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "vz0L")
                    vz0L = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "kTzL")
                    kTzL = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "kTpL ")
                    kTpL = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "n0R")
                    n0R = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "vz0R")
                    vz0R = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "kTzR")
                    kTzR = tryparse(Float64, mydata[i_lines][i+1:j-1])
                elseif contains(mydata[i_lines][1:i], "kTpR")
                    kTpR = tryparse(Float64, mydata[i_lines][i+1:j-1])
                end
            end
        end
        particle[ii] = (; Nvz, vzmin, vzmax, Nmu, mumin, mumax, muexp, mass, charge,
                            n0, vz0, kTz, kTp, n0L, vz0L, kTzL, kTpL, n0R, vz0R, kTzR, kTpR)
    end

    return dt, Niter, dump_period_fields, fields_per_file, dump_period_distr, dump_period_distr_IonoBoundary,
            dump_period_distr_1v, dump_start, shift_test_period, resistance, Nz, zmin, zmax,
            Nspecies, const_a, BC_Poisson, voltage, initialiser, voltage_init, E0, startfromdumpfile,
            transffilename, voltagefilename, particle
end


        













# function try_find_variable(file_line, i, j, variable)
#     if contains(file_line[1:i], variable)
#         var = tryparse(Float64, file_line[i+1:j-1])
#         return var
#     end
# end

# tryparse(Float64, mydata[2][4+1:13-1])
# try_find_variable(mydata[i_lines], i, j, "dt")







# mydata = readdlm(path_to_my_data, '=', comments=true, comment_char = '%')
# strip.(mydata[1,2], [';'])
# findfirst("%END", mydata)

# name, val = read_input(path_to_my_data)
# using DataFrames
# Input = DataFrame(; name, val);

## -------------------------------------------------------------------------------------- ##



# dir = "/mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7/Nz_3000_150to300s_results/"

# ## ----- Input data ----- ##
# path_to_my_data = joinpath(dir, "inputb6.m"); 
# mydata2 = read(path_to_my_data, String); # read inputb6.m as a continuous string

# # Get the general parameters
# m = collect(eachmatch(r"\n.*= *(\d+.*);", mydata2[1:findfirst("%END", mydata2)[end]]))
# dt =                    tryparse(Float64, m[1].captures[])
# Niter =                 tryparse(Int64, m[2].captures[])
# dump_period_fields =    tryparse(Int64, m[3].captures[])
# fields_per_file =       tryparse(Int64, m[4].captures[])
# dump_period_distr =     tryparse(Int64, m[5].captures[])
# dump_period_distr_1v =  tryparse(Int64, m[6].captures[])
# dump_start =            tryparse(Int64, m[7].captures[])
# shift_test_period =     tryparse(Int64, m[8].captures[])
# resistance =            tryparse(Float64, m[9].captures[])
# Nz =                    tryparse(Int64, m[10].captures[])
# zmin =                  tryparse(Float64, m[11].captures[])
# zmax =                  tryparse(Float64, m[12].captures[])
# Nspecies =              tryparse(Int64, m[13].captures[])
# const_a =               tryparse(Float64, m[14].captures[])
# BC_Poisson =            tryparse(Int64, m[15].captures[])
# voltage =               tryparse(Int64, m[16].captures[])
# initialiser =           tryparse(Int64, m[17].captures[])
# E0 =                    tryparse(Float64, m[18].captures[])
# dump_period_dump =      tryparse(Int64, m[19].captures[])

# # Get the species specific parameters
# particle = Array{NamedTuple}(undef, 4)
# for i in 1:Nspecies
#     start = findfirst("%SPEC", mydata2)[end]
#     finish = findnext("%END", mydata2, start)[end]
#     m = collect(eachmatch(r"\n[^%].*= *(-?\d+.*);", mydata2[start:finish]))
#     Nvz =       tryparse(Int64, m[1].captures[])
#     vzmin =     tryparse(Float64, m[2].captures[])
#     vzmax =     tryparse(Float64, m[3].captures[])
#     Nmu =       tryparse(Int64, m[4].captures[]) 
#     mumin =     tryparse(Int64, m[5].captures[])
#     mumax =     tryparse(Float64, m[6].captures[])
#     muexp =     tryparse(Int64, m[7].captures[])
#     mass =      tryparse(Float64, m[8].captures[])
#     charge =    tryparse(Float64, m[9].captures[])
#     n0 =        tryparse(Float64, m[10].captures[])
#     vz0 =       tryparse(Float64, m[11].captures[])
#     kTz =       tryparse(Float64, m[12].captures[])
#     kTp =       tryparse(Float64, m[13].captures[])
#     n0L =       tryparse(Float64, m[14].captures[])
#     vz0L =      tryparse(Float64, m[15].captures[])
#     kTzL =      tryparse(Float64, m[16].captures[])
#     kTpL =      tryparse(Float64, m[17].captures[])
#     n0R =       tryparse(Float64, m[18].captures[])
#     vz0R =      tryparse(Float64, m[19].captures[])
#     kTzR =      tryparse(Float64, m[20].captures[])
#     kTpR =      tryparse(Float64, m[21].captures[])

#     particle[i] = (; Nvz, vzmin, vzmax, Nmu, mumin, mumax, muexp, mass, charge,
#                             n0, vz0, kTz, kTp, n0L, vz0L, kTzL, kTpL, n0R, vz0R, kTzR, kTpR)
# end

# # Find the number of precessors used from the job.bat file, if there is one
# # I don't understand that part ?

# # Construct XI-vectors
# dxi = 1 / Nz
# xicorn = dxi * 0:Nz
# xi = 0.5 * (xicorn[1:end-1] + xicorn[2:end])
# # Compute transformations, i.e. z-vectors and g'


# #=
# There is much more to do. Problem is that we would need to translate the poleslevel scripts and more...
# Too much time to lose here. Probably 2-3 days...
# =#