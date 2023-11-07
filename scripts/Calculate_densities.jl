# Will converts e- flux `Ie` (#e⁻/cm²/s) into number density `n_e` (#e⁻/m³).
# - Ie : electron flux (#e⁻/cm²/s), 3D array [n_beam * nz, nt, nE]
# - n_e: electron density (#e⁻/m³), 3D array [nz, nt, nE]

using AURORA

directory_to_process = "Visions2/Alfven_536s_correct_msis_and_scattering"
make_density_file(directory_to_process)
