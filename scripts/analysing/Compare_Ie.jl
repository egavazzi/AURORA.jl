using MAT
using AURORA


directory_to_process1 = "Benchmark/500km_7000eV_0-01s_megaturbo"
directory_to_process2 = "Benchmark/500km_3000eV_0-01s(3)"


##
data1 = AURORA.load_results(directory_to_process1)
data2 = AURORA.load_results(directory_to_process2)

##
data1.Ie ≈ data2.Ie
