module RD2

using JLD2, DataFrames
@info join(["Loading data",
            "May take a few minutes"],
           "\n")
dat_dir = "/home/xijiang/Music/workspace/data/RD2"

@load "$dat_dir/cv-setup.jld"     dcs gcs ncs
@load "$dat_dir/drp-training.jld" dts gts nts
@load "$dat_dir/ebv-all.jld"      dbv gbv nbv
@load "$dat_dir/GAID.jld"         G A ID;

include("workflow.jl")
include("milk.jl")
# include("compare.jl")

end # module
