using AURORA
using MAT
using GLMakie
GLMakie.activate!()
# using CairoMakie
# CairoMakie.activate!()

# directory_to_process1 = "CFL_simulations/7keV_CFL-128"
# directory_to_process1 = "test_matlab_comparisons/PnLh-05-3000-3kev-FA"
directory_to_process1 = "Visions2/Alfven_536s_correct_msis_and_scattering"
# directory_to_process2 = "CFL_simulations/7keV_CFL-32"
# directory_to_process2 = "../../../test_matlab_AURORA/PnLh-05-3000-3keV-FA"
directory_to_process2 = "Visions2/Alfven_536s_correct_msis"
directory_to_process3 = "CFL_simulations/7keV_CFL-256"


##
# Load the J data (parallel current)
# full_path_to_directory = pkgdir(AURORA, "data", directory_to_process1)
full_path_to_directory = REVONTULI_MOUNT * "/mnt/data/etienne/Julia/AURORA.jl/data/" * directory_to_process1
J_file = joinpath(full_path_to_directory, "J.mat")
data = matread(J_file)
J_up1 = vec(data["J_up"])
J_down1 = vec(data["J_down"])
IeEup1 = vec(data["IeEup"])
IeEdown1 = vec(data["IeEdown"])
t1 = vec(data["t"])

# Load the J data (parallel current)
# full_path_to_directory = REVONTULI_MOUNT * "/mnt/data/etienne/Julia/AURORA.jl/data/" * directory_to_process2
full_path_to_directory = REVONTULI_MOUNT * "/mnt/data/etienne/Julia/AURORA.jl/data/" * directory_to_process2
J_file = joinpath(full_path_to_directory, "J_top-1.mat")
data = matread(J_file)
J_up2 = vec(data["J_up"])
J_down2 = vec(data["J_down"])
IeEup2 = vec(data["IeEup"])
IeEdown2 = vec(data["IeEdown"])
t2 = vec(data["t"])

# Load the J data (parallel current)
full_path_to_directory = REVONTULI_MOUNT * "/mnt/data/etienne/Julia/AURORA.jl/data/" * directory_to_process3
J_file = joinpath(full_path_to_directory, "J.mat")
data = matread(J_file)
J_up3 = vec(data["J_up"])
J_down3 = vec(data["J_down"])
IeEup3 = vec(data["IeEup"])
IeEdown3 = vec(data["IeEdown"])
t3 = vec(data["t"])


##
f = Figure()
ax = Axis(f[1, 1], xlabel = "time (s)", ylabel = "")
sc1 = scatterlines!(ax, t1, IeEup1)
sc2 = scatterlines!(ax, t1, IeEup2)
# sc3 = scatterlines!(ax, t1, IeEup3)
axislegend(ax, [sc1, sc2, sc3], ["128", "32", "256"], "CFL", position=:rb)
DataInspector(f)
f
