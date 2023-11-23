# This script extract Ie from simulation results, downsample it in time, and save the
# downsampled version in a new file
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
