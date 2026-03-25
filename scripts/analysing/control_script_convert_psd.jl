# =============================================================================
# Batch conversion of AURORA IeFlickering-NN.mat files into phase-space density.
# =============================================================================

using AURORA

directories_to_process = ["data/Visions2/Alfven_train_335s"]

for current_dir in directories_to_process
    println("Processing $(basename(current_dir))...")
    make_psd_file(
        current_dir;
        compute = :F_only,        # :f_only, :F_only, or :both
        vpar_edges = nothing,
        output_prefix = "psd",  # will generate output files (`<output_prefix>-NN.mat`)
    )
end

println("Done.")
