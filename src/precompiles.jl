using PrecompileTools

let
    @setup_workload begin
        E, dE = make_energy_grid(3000)
        θ_lims = 180:-45:0
        @compile_workload begin
            AURORA.load_cross_sections(E, dE)
            AURORA.load_scattering_matrices(θ_lims, 10; verbose = false)
            AURORA.cosd.(θ_lims);
            AURORA.mu_avg(θ_lims);
            AURORA.beam_weight(θ_lims);
        end
    end
end
