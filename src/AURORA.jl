module AURORA

include("constants.jl")

include("grids/abstract_grid.jl")
include("grids/altitude_grid.jl")
include("grids/energy_grid.jl")
include("grids/pitch_angle_grid.jl")
export AbstractGrid, AltitudeGrid, EnergyGrid, PitchAngleGrid

include("ionosphere/ionosphere.jl")
export Ionosphere
include("ionosphere/iri/iri.jl")
include("ionosphere/msis/msis.jl")
export find_msis_file, find_nrlmsis_file
export find_iri_file

include("physics/cross_sections/e_N2_cross_sections.jl")
include("physics/cross_sections/e_O2_cross_sections.jl")
include("physics/cross_sections/e_O_cross_sections.jl")
include("physics/cross_sections/emission_cross_sections.jl")
include("physics/cross_sections/cross_sections.jl")

include("physics/cache_policy.jl")
export CachePolicy
include("physics/scattering.jl")
include("physics/scattering_cache.jl")
export ScatteringData, clear_scattering_cache!
include("physics/phase_functions.jl")
export phase_fcn_N2, phase_fcn_O2, phase_fcn_O, convert_phase_fcn_to_3D

include("physics/laws.jl")
export ExprLaw, @law

include("physics/cascading.jl")
include("physics/cascading_cache.jl")
export clear_cascading_cache!

include("physics/species.jl")
export NeutralSpecies, MSISDensity, VectorDensity
export N2Species, O2Species, OSpecies

include("model.jl")
export AuroraModel, make_altitude_grid, make_energy_grid

include("input/spectra.jl")
include("input/modulations.jl")
include("input/input_flux.jl")
export AbstractSpectrum, FlatSpectrum, GaussianSpectrum, MaxwellianSpectrum, FileSpectrum
export AbstractModulation, ConstantModulation, SinusoidalFlickering, SquareFlickering, SmoothOnset
export InputFlux, evaluate_spectrum, apply_modulation, compute_flux
export Ie_top_from_file

include("solvers/transport_matrices.jl")
include("solvers/matrix_building.jl")

include("solvers/energy_degradation.jl")

include("solvers/sparse_indexing.jl")
include("solvers/steady_state.jl")
include("solvers/crank_nicolson.jl")

include("simulation/workspace.jl")
include("output/output_manager.jl")
export AuroraOutputManager
include("simulation/loop_planning.jl")
include("simulation/types.jl")
export AuroraSimulation
export AbstractMode, SteadyStateMode, TimeDependentMode, SteadyState, TimeDependent
export AbstractTimeConfig, SingleStepConfig, UniformTimeGrid, RefinedTimeGrid
include("simulation/initialize.jl")
export initialize!
include("output/write.jl")
include("output/read.jl")
export SimulationResult, load_results, load_coordinates
include("simulation/run.jl")
export run!

include("utilities.jl")
export v_of_E, mu_avg, beam_weight

include("analysis/analysis_types.jl")
export VolumeExcitationResult, ColumnExcitationResult, IeTopResult,
       load_volume_excitation, load_column_excitation, load_Ie_top
include("analysis/psd.jl")
include("analysis/emissions.jl")
include("analysis/fluxes.jl")
include("analysis/heating.jl")
export make_volume_excitation_file, make_column_excitation_file,
       make_Ie_top_file, make_current_file,
       make_heating_rate_file, make_psd_file

include("viz_interface.jl")

# Precompile selected functions
include("precompiles.jl")

end
