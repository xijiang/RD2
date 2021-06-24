module RD2

using JLD2, DataFrames, SparseArrays, LinearAlgebra, Statistics, StatsPlots,
    Random, CodecZlib

dat_dir = "/home/xijiang/Music/workspace/data/RD2"

include("workflow.jl")
include("milk.jl")
include("compare.jl")
include("cv-stage-1.jl")
include("refresh-data.jl")

# Test procedures that can serve as examples
# and can be ignored in package of release version
#include("tst/cv.jl")
#include("tst/cv23.jl")
#include("tst/cv-210422.jl")

export workflow

end # module
