using PrecompileTools: @setup_workload, @compile_workload

let
    @setup_workload begin
        E, dE = make_energy_grid(3000)
        θ_lims = 180:-45:0
        @compile_workload begin
            load_cross_sections(E, dE)
            load_scattering_matrices(θ_lims, 10; verbose = false)
            cosd.(θ_lims);
            mu_avg(θ_lims);
            beam_weight(θ_lims);
        end
    end
end
