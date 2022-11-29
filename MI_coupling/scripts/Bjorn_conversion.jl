# This script converts Ie into fzvzv⟂ over altitude and time.
using Aurora
using Printf
using MAT

HMR_E = 40  # How Much to Refine E
HMR_MU = 80 # How Much to Refine μ (= cos(θ))

main_dir = "/mnt/data/bjorn/Campaigns/Etrp/"
directories = readdir(main_dir, join=true)
directories_to_process = contains.(directories, "AWA")
# loop through the directories
for i in eachindex(directories[directories_to_process])
    current_directory = directories[directories_to_process][i]
    println("Processing $current_directory...")
    files = readdir(current_directory, join=true)
    files_to_process = contains.(files, r"IeFlickering\-[0-9]+\.mat")
    savedir = joinpath(current_directory, "f_distr")
    if !isdir(savedir)
        mkdir(savedir) # create subfolder where results will be saved
    end
    # loop through the IeFlickering files
    for j in eachindex(files[files_to_process])
        fzvzvperp, E_grid, t_run, h_atm, μ_pitch_grid,
            vz, Δvz, vperp, Δvperp = make_f_from_Aurora(files[files_to_process][j], HMR_E, HMR_MU);

        save_file = joinpath(savedir, (@sprintf( "fzvzvperp-%02d.mat", j)))
        f = matopen(save_file, "w")
            write(f, "fzvzvperp", fzvzvperp)
            write(f, "E", E_grid)
            write(f, "t_run", t_run)
            write(f, "h_atm", h_atm)
            write(f, "mu_lims", μ_pitch_grid)

            write(f, "vz", vz)
            write(f, "vperp", vperp)
            write(f, "dvz", Δvz)
            write(f, "dvperp", Δvperp)
        close(f)
    end
 end
