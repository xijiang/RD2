#=
Running order:
1. refresh_norsk_dt_210423()
2. cv_210510()
=#
using DataFrames, JLD2, CSV, GZip, LinearAlgebra, SparseArrays, Statistics, CodecZlib
import RD2:dat_dir
"""
    function read_n_drp(drp)
---
This is to read data sent by Mario on 2021-04-23.
More traits were added.
"""
function read_n_drp(drp)
    @info "  - Norwegian data training set"
    np = DataFrame(ID          = String[],
                   Milk        = Float64[],
                   Milk_ERC    = Float64[],
                   Fat         = Float64[],
                   Fat_ERC     = Float64[],
                   Protein     = Float64[],
                   Protein_ERC = Float64[],
                   SCS         = Float64[],
                   SCS_ERC     = Float64[],
                   x10         = Float64[], # NoInsHeifer DRP     
                   x11         = Float64[], # NoInsHeifer ERC     
                   x12         = Float64[], # NoInsCow DRP        
                   x13         = Float64[], # NoInsCow ERC        
                   x14         = Float64[], # calving interval DRP
                   x15         = Float64[], # calving interval ERC
                   N           = Int[])
    for line in eachline(drp)
        v = split(line)
        length(v) != 16 && continue # because two lines have only 2 columns in trnn set
        id = v[1]
        y = parse.(Float64, v[2:end-1])
        n = parse(Int, v[end])
        push!(np, [id; y; n])
    end
    np
end

function refresh_norsk_dt_210423()
    @info "Read phnotypes, EBV and CV sets"
    @load "$dat_dir/jld/drp-training-2021-03-11.jld" dts gts nts
    nts = read_n_drp(joinpath(dat_dir, "mario/20210423/norway_DRPs.txt"))
    @save joinpath(dat_dir, "jld/drp-training-2021-04-23.jld") {compress=true} dts gts nts

    @info "Insert new breeding values of Norwegian data"
    df = CSV.File(joinpath(dat_dir, "data-2020-07-02/norsk_RDC_REL_RBV.csv")) |> DataFrame
    @load "$dat_dir/jld/ebv-all.jld" dbv gbv nbv
    df.ID = string.(df.GenoId)
    df1 = select(df, :ID,
                 :r_Idx_KgProtein305d1_3=>:Protein,
                 :r_Idx_KgFat305d1_3=>:Fat)
    df1[df1.Fat .== ".", :Fat] .= "999"
    df1[df1.Protein .== ".", :Protein] .= "999"
    df1.Fat = parse.(Float64, df1.Fat)
    df1.Protein = parse.(Float64, df1.Protein)
    nbv = innerjoin(nbv, df1, on=:ID)
    @save joinpath(dat_dir, "jld/ebv_202105.jld") {compress=true} dbv gbv nbv
end

function ccwcv(yt, wt, yv, Gi, X)
    nt = length(yt)
    nv = length(yv)
    nx = size(X)[2] + 1

    #Ri = diagm(wt)
    XR = X .* wt
    Z  = sparse([I zeros(nt, nv)])
    ZR = Z .* wt
    lhs = [ XR'X XR'Z
            ZR'X ZR'Z + Gi ] # as va = ve = 1, h2=0.5.
    rhs = [ XR'yt
            ZR'yt ]
    bv = lhs \ rhs
    cc = cor(yv, bv[nt+nx:end])
    return cc, bv[nt+nx:end]
end

function cccv(yt, yv, Gi, X)
    nt = length(yt)
    nv = length(yv)
    nx = size(X)[2] + 1

    Z = sparse([I zeros(nt, nv)])
    lhs = [ X'X X'Z
            Z'X Z'Z + Gi ]      # as va = ve = 1, h2=0.5.
    rhs = [ X'yt
            Z'yt ]
    bv = lhs \ rhs
    cc = cor(yv, bv[nt+nx:end])
    return cc, bv[nt+nx:end]
end

function cv_G(ts, vs, G, trait, erc, X)
    id = ts[:, :ix]
    yt = ts[:, trait]
    wt = ts[:, erc]
    nt = length(id)

    yv = vs[:, trait]
    append!(id, vs.ix)
    nv = length(yv)

    gi = inv(G[id, id])
    c1, _ = ccwcv(yt, wt, yv, gi, X)
    c2, _ = cccv(yt, yv, gi, X)
    return c1, c2
end

function make_X(len...)
    nc = length(len)
    nr = sum(len)
    X = zeros(nr, nc)
    ac = 0
    for i in 1:nc
        X[ac+1:ac+len[i], i] .=1
        ac += len[i]
    end
    X
end

function mk_cv_tbl(dts, dvs, gts, gvs, nts, nvs, G, trait, erc)
    # target is each country in turn
    rst = zeros(6, 7)           # 6 = w/nw by (DGN); 7 = D G N DG DN GN DGN

    @info "Prediction of one country with one of 3 coutries"
    # D → D
    X = make_X(nrow(dts))
    rst[1, 1], rst[4, 1] = cv_G(dts, dvs, G, trait, erc, X)
    # G → D
    X = make_X(nrow(gts))
    rst[1, 2], rst[4, 2] = cv_G(gts, dvs, G, trait, erc, X)
    # N → D
    X = make_X(nrow(nts))
    rst[1, 3], rst[4, 3] = cv_G(nts, dvs, G, trait, erc, X)
    # D → G
    X = make_X(nrow(dts))
    rst[2, 1], rst[5, 1] = cv_G(dts, gvs, G, trait, erc, X)
    # G → G
    X = make_X(nrow(gts))
    rst[2, 2], rst[5, 2] = cv_G(gts, gvs, G, trait, erc, X)
    # N → G
    X = make_X(nrow(nts))
    rst[2, 3], rst[5, 3] = cv_G(nts, gvs, G, trait, erc, X)
    # D → G
    X = make_X(nrow(dts))
    rst[3, 1], rst[6, 1] = cv_G(dts, nvs, G, trait, erc, X)
    # G → G
    X = make_X(nrow(gts))
    rst[3, 2], rst[6, 2] = cv_G(gts, nvs, G, trait, erc, X)
    # N → G
    X = make_X(nrow(nts))
    rst[3, 3], rst[6, 3] = cv_G(nts, nvs, G, trait, erc, X)

    @info "Prediction of one country with two of 3 coutries"
    dgts = vcat(dts, gts)
    dnts = vcat(dts, nts)
    gnts = vcat(gts, nts)
    # D + G → D
    X = make_X(nrow(dts), nrow(gts))
    rst[1, 4], rst[4, 4] = cv_G(dgts, dvs, G, trait, erc, X)
    # D + N → D
    X = make_X(nrow(dts), nrow(nts))
    rst[1, 5], rst[4, 5] = cv_G(dnts, dvs, G, trait, erc, X)
    # G + N → D
    X = make_X(nrow(gts), nrow(nts))
    rst[1, 6], rst[4, 6] = cv_G(gnts, dvs, G, trait, erc, X)
    # D + G → G
    X = make_X(nrow(dts), nrow(gts))
    rst[2, 4], rst[5, 4] = cv_G(dgts, gvs, G, trait, erc, X)
    # D + N → G
    X = make_X(nrow(dts), nrow(nts))
    rst[2, 5], rst[5, 5] = cv_G(dnts, gvs, G, trait, erc, X)
    # G + N → G
    X = make_X(nrow(gts), nrow(nts))
    rst[2, 6], rst[5, 6] = cv_G(gnts, gvs, G, trait, erc, X)
    # D + G → N
    X = make_X(nrow(dts), nrow(gts))
    rst[3, 4], rst[6, 4] = cv_G(dgts, nvs, G, trait, erc, X)
    # D + N → N
    X = make_X(nrow(dts), nrow(nts))
    rst[3, 5], rst[6, 5] = cv_G(dnts, nvs, G, trait, erc, X)
    # G + N → N
    X = make_X(nrow(gts), nrow(nts))
    rst[3, 6], rst[6, 6] = cv_G(gnts, nvs, G, trait, erc, X)

    @info "Prediction of one country with all 3 coutries"
    dgnts = vcat(dts, gts, nts)
    # D + G + N → D
    X = make_X(nrow(dts), nrow(gts), nrow(nts))
    rst[1, 7], rst[4, 7] = cv_G(dgnts, dvs, G, trait, erc, X)
    # D + G + N → G
    X = make_X(nrow(dts), nrow(gts), nrow(nts))
    rst[2, 7], rst[5, 7] = cv_G(dgnts, gvs, G, trait, erc, X)
    # D + G + N → N
    X = make_X(nrow(dts), nrow(gts), nrow(nts))
    rst[3, 7], rst[6, 7] = cv_G(dgnts, nvs, G, trait, erc, X)

    round.(rst, sigdigits = 2)
end

function make_tv(ts, cs, bv, drp, erc)
    t1 = filter(row ->row.st == 't', cs)
    df1 = select(innerjoin(t1, ts, on = :ID), :ix, drp, erc)
    dft = filter(row -> row[drp] < 999, df1)
    t2  = filter(row -> row.st == 'v', cs)
    dfv = select(innerjoin(t2, bv, on = :ID), :ix, drp)
    return dft, dfv
end

function cv_210510()
    @load "$dat_dir/jld/cv-setup.jld" dcs gcs ncs
    @load "$dat_dir/jld/drp-training-2021-04-23.jld" dts gts nts
    @load "$dat_dir/jld/ebv_202105.jld"  dbv gbv nbv
    @load "$dat_dir/jld/GAID.jld"     G A ID
    
    begin                       # The milk trait
        dt, dv = make_tv(dts, dcs, dbv, :Milk, :Milk_ERC)
        gt, gv = make_tv(gts, gcs, gbv, :Milk, :Milk_ERC)
        nt, nv = make_tv(nts, ncs, nbv, :Milk, :Milk_ERC)
        rst = mk_cv_tbl(dt, dv, gt, gv, nt, nv, G, :Milk, :Milk_ERC)
        for row in eachrow(rst)
            println("| ", join(row, " | "), " |")
        end
    end
    begin                       # The protein trait
        dt, dv = make_tv(dts, dcs, dbv, :Protein, :Protein_ERC)
        gt, gv = make_tv(gts, gcs, gbv, :Protein, :Protein_ERC)
        nt, nv = make_tv(nts, ncs, nbv, :Protein, :Protein_ERC) 
        rst = mk_cv_tbl(dt, dv, gt, gv, nt, nv, G, :Protein, :Protein_ERC)
        for row in eachrow(rst)
            println("| ", join(row, " | "), " |")
        end
    end
    begin                       # The fat trait
        dt, dv = make_tv(dts, dcs, dbv, :Fat, :Fat_ERC)
        gt, gv = make_tv(gts, gcs, gbv, :Fat, :Fat_ERC)
        nt, nv = make_tv(nts, ncs, nbv, :Fat, :Fat_ERC)
        rst = mk_cv_tbl(dt, dv, gt, gv, nt, nv, G, :Fat, :Fat_ERC)
        for row in eachrow(rst)
            println("| ", join(row, " | "), " |")
        end
    end
    begin                       # The SCS trait
        dt, dv = make_tv(dts, dcs, dbv, :SCS, :SCS_ERC)
        gt, gv = make_tv(gts, gcs, gbv, :SCS, :SCS_ERC)
        nt, nv = make_tv(nts, ncs, nbv, :SCS, :SCS_ERC)
        rst = mk_cv_tbl(dt, dv, gt, gv, nt, nv, G, :SCS, :SCS_ERC)
        for row in eachrow(rst)
            println("| ", join(row, " | "), " |")
        end
    end
end
