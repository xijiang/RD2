"""
    function cv_milk(ts, vs, G)
---
Cross validation with only G matrix
"""
function cv_milk(ts, vs, G)
    id = ts[:, :ix]
    yt = ts[:, :Milk]
    wt = ts[:, :Milk_ERC]
    nt = length(id)

    yv = vs[:, :Milk]
    append!(id, vs.ix)
    nv = length(yv)

    gi = inv(G[id, id])
    #ccwcv(yt, wt, yv, gi)
    cccv(yt, yv, gi)
end

"""
Test all the combinations to predict a 3rd country with 
- 3 x 1→1
- 3 x 2→2
- 1 x 3→1
"""
function cv_21_04_08()
    @load "$dat_dir/jld/cv-setup.jld" dcs gcs ncs
    @load "$dat_dir/jld/drp-training-2021-03-11.jld" dts gts nts
    @load "$dat_dir/jld/ebv-all.jld"  dbv gbv nbv
    @load "$dat_dir/jld/GAID.jld"     G A ID

    dts, dvs = mk_milk_dt(dts, dcs, dbv)
    gts, gvs = mk_milk_dt(gts, gcs, gbv)
    nts, nvs = mk_milk_dt(nts, ncs, nbv)
    
    @info "Prediction of one country with one of 3 coutries"
    # D → D
    cv_milk(dts, dvs, G)
    # G → D
    cv_milk(gts, dvs, G)
    # N → D
    cv_milk(nts, dvs, G)
    # D → G
    cv_milk(dts, gvs, G)
    # G → G
    cv_milk(gts, gvs, G)
    # N → G
    cv_milk(nts, gvs, G)
    # D → G
    cv_milk(dts, nvs, G)
    # G → G
    cv_milk(gts, nvs, G)
    # N → G
    cv_milk(nts, nvs, G)

    @info "Prediction of one country with two of 3 coutries"
    dgts = vcat(dts, gts)
    dnts = vcat(dts, nts)
    gnts = vcat(gts, nts)
    # D + G → D
    cv_milk(dgts, dvs, G)
    # D + N → D
    cv_milk(dnts, dvs, G)
    # G + N → D
    cv_milk(gnts, dvs, G)
    # D + G → G
    cv_milk(dgts, gvs, G)
    # D + N → G
    cv_milk(dnts, gvs, G)
    # G + N → G
    cv_milk(gnts, gvs, G)
    # D + G → N
    cv_milk(dgts, nvs, G)
    # D + N → N
    cv_milk(dnts, nvs, G)
    # G + N → N
    cv_milk(gnts, nvs, G)
    @info "Prediction of one country with all 3 coutries"
    dgnts = vcat(dts, gts, nts)
    # D + G + N → D
    cv_milk(dgnts, dvs, G)
    # D + G + N → G
    cv_milk(dgnts, gvs, G)
    # D + G + N → N
    cv_milk(dgnts, nvs, G)
end
