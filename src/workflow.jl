"""
Records the test procedure I performed for data of stage II.
Notes,
- it has been demonstrated that DRP, and EBV are not good for training.
- The source codes for above can be found in the data directory.
- This package starts from using deregressed EBV (GSE, 2009, 41:55).
"""
function workflow()
    @load "$dat_dir/cv-setup.jld"     # => dcs gcs ncs
    @load "$dat_dir/drp-training.jld" # => dts gts nts
    @load "$dat_dir/ebv-all.jld"      # => dbv gbv nbv
    @load "$dat_dir/GAID.jld"         # => G A ID
end
