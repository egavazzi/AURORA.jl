using PrecompileTools

let
    @setup_workload begin
        E, dE = make_energy_grid(3000)
        @compile_workload begin
            AURORA.load_cross_sections(E, dE)
        end
    end
end
