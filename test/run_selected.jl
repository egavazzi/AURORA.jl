using AURORA
using TestItemRunner

@run_package_tests verbose=true filter=ti->(startswith(ti.filename, joinpath(@__DIR__, "atmosphere")))
