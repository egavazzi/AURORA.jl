# This script extract Ie from simulation results, downsample it in time, and save the
# downsampled version in a new file

using AURORA
using MAT

directory_to_process = "backup/20231027-0958"
downsampling = 10

full_path_to_directory = pkgdir(AURORA, "data", directory_to_process)
files = readdir(full_path_to_directory, join=true)
files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]

for j in files_to_process
    f = matopen(j)
        Ie = read(f, "Ie_ztE") # [n_μ * nz, nt, nE]
        E = read(f, "E")
        t_run = read(f, "t_run")
        μ_lims = read(f, "mu_lims")
        h_atm = read(f, "h_atm")
        I0 = read(f, "I0")
        μ_scatterings = read(f, "mu_scatterings")
    close(f)

    # downsample the data
    dt = diff(t_run)[1]
    println("The time-step from simulation is ", dt, "s.")
    new_dt = dt * downsampling
    println("The time-step of the new file will be ", new_dt, "s.")
    Ie = Ie[:, 1:downsampling:end, :]
    t_run = t_run[:, 1:downsampling:end]

    # create new subdir if it doesn't exist
    new_subdir = "downsampled_" * string(downsampling) * "x"
    full_path_to_new_subdir = joinpath(full_path_to_directory, new_sub_dir)
    mkpath(full_path_to_new_subdir)

    # create new file
    new_filename = splitdir(j)[2][1:end-4] * "d.mat" # add a 'd' after the number to indicate that it is downsampled
    full_path_to_new_filename = joinpath(full_path_to_new_subdir, new_filename)
    file = matopen(full_path_to_new_filename, "w")
        write(file, "Ie_ztE", Ie)
        write(file, "E", E)
        write(file, "t_run", t_run)
        write(file, "mu_lims", μ_lims)
        write(file, "h_atm", h_atm)
        write(file, "I0", I0)
        write(file, "mu_scatterings", μ_scatterings)
    close(file)
end




# X = 0:0.001:0.1
# X = (0:0.001:0.1) .+ 0.1
# dt = 0.001
# new_dt = dt * 10
# X[X .∈ [0:new_dt:X[end]]]

# X[mod.(X, 0.01) .== 0]
