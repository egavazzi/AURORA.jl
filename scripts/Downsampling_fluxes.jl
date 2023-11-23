# This script calls the function `downsampling_fluxes(directory_to_process, downsampling_factor)`
# which extracts Ie from simulation results in `directory_to_process`, downsample it in time
# according to `downsampling_factor`, and save the downsampled version in a new subdirectory
# `downsampled_Xx`. 'X' being the downsampling factor.
#
# - directory_to_process: relative path to the directory to process. Should start in `data`
# directory. Example : `directory_to_process = "backup/20231123-0912"`
#
# - downsampling_factor: downsampling factor in time.

using AURORA

directory_to_process = "backup/20231123-0912"
downsampling_factor = 10

# calling the function that does the work
downsampling_fluxes(directory_to_process, downsampling_factor)
