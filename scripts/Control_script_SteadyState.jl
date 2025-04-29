using MAT
using Pkg
AURORA_dir = "/mnt/data/oliver/AURORA.jl/"
Pkg.activate(AURORA_dir)
using AURORA



## Setting parameters
altitude_lims = [60, 600];     # (km) altitude limits of the ionosphere
θ_lims = 180:-10:0;            # (°) angle-limits for the electron beams
#E_max = 3000;                   # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith


h_atmosphere_models = altitude_lims[1]-9:1:altitude_lims[2]+200
msis_file = find_nrlmsis_file(
    year=2006, month=12, day=12, hour=19, minute=30, lat=70, lon=19, height=h_atmosphere_models
    );
iri_file = find_iri_file(
    year=2006, month=12, day=12, hour=19, minute=30, lat=70, lon=19, height=h_atmosphere_models
    );


## Define where to save the results
root_savedir = "ionprod_v6"  # name of the root folder



## Define input parameters
input_type = "constant_onset"
IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
z₀ = altitude_lims[2];          # (km) altitude of the source
#E_min = E_max - 100;         # (eV) bottom energy of the FAB
Beams = 1:length(θ_lims)-1;   # beam numbers for the precipitation, starting with field aligned down
t0 = 0;                     # (s) time of start for smooth transition
t1 = 0;                     # (s) time of end for smooth transition


e_grid = 10 .^(range(1,stop=5,length=400))


for i in length(e_grid)-1:-1:1
    E_max = e_grid[i+1]
    E_min = e_grid[i]
    name_savedir = string(E_min)*'-'*string(E_max) * "eV"
    savedir = make_savedir(root_savedir, name_savedir)

    INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);

    ## Run the simulation
    calculate_e_transport_steady_state(altitude_lims, θ_lims, E_max, B_angle_to_zenith,
        msis_file, iri_file, savedir, INPUT_OPTIONS);
end


## Run the analysis
directory_to_process = joinpath(AURORA_dir, "data", root_savedir)
subdirs = readdir(directory_to_process, join=true)
for dir in subdirs
    make_volume_excitation_file(dir)
end



dE = diff(e_grid)                                   
E = e_grid[1:end-1] + dE/2                          # already centre energy bins

h_atm = []
QN2i      = []
QO1D      = []
QO1S      = []
QO2i      = []
QOi       = []
Q_from_O  = []
Q_from_O2 = []
Q_from_N2 = []

for i in length(e_grid)-1:-1:353
    E_max = e_grid[i+1]
    E_min = e_grid[i]
   
    name_savedir = string(E_min)*'-'*string(E_max) * "eV"
    f = matopen(joinpath(directory_to_process, name_savedir, "Qzt_all_L.mat"))
        push!(h_atm, read(f, "h_atm"))
        push!(QN2i      , read(f, "QN2i"     ))
        push!(QO1D      , read(f, "QO1D"     ))
        push!(QO1S      , read(f, "QO1S"     ))
        push!(QO2i      , read(f, "QO2i"     ))
        push!(QOi       , read(f, "QOi"      ))
        push!(Q_from_O  , read(f, "Q_from_O" ))
        push!(Q_from_O2 , read(f, "Q_from_O2"))
        push!(Q_from_N2 , read(f, "Q_from_N2"))
    close(f)
end

h_atm = h_atm[1]
QN2i      = reduce(hcat, QN2i     )
QO1D      = reduce(hcat, QO1D     )
QO1S      = reduce(hcat, QO1S     )
QO2i      = reduce(hcat, QO2i     )
QOi       = reduce(hcat, QOi      )
Q_from_O  = reduce(hcat, Q_from_O )
Q_from_O2 = reduce(hcat, Q_from_O2)
Q_from_N2 = reduce(hcat, Q_from_N2)


fig = Figure();
sc = display(fig);
ax = Axis(fig[1, 1], title = "Figure 1", xlabel = "x", ylabel = "y")
lines!(Q_from_O[:, 1], h_atm/1e3)
current_figure()
save(joinpath(directory_to_process, name_savedir, "figure.pdf"), fig, pdf_version="1.4")
p = lines(QN2i, h_atm/1e3, label = e_grid[i])


#to do: plot production profiles, assemble production in matrix