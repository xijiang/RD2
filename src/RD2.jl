module RD2

using JLD2, DataFrames, SparseArrays, LinearAlgebra, Statistics, StatsPlots,
    Random

dat_dir = "/home/xijiang/Music/workspace/data/RD2"

include("workflow.jl")
include("milk.jl")
include("compare.jl")
include("cv-stage-1.jl")
include("refresh-data.jl")
include("tst/cv.jl")

export workflow

end # module
