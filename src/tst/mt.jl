#=
Running order:
1. refresh_norsk_dt_210423()
2. cv_210510()
=#
using DataFrames, JLD2, CSV, GZip, LinearAlgebra, SparseArrays, Statistics, CodecZlib
import RD2:dat_dir

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

"""
    function mtwcv(y1, y2, yv, wt, ct, Gi, X)
---
## `M`ultiple `t`rait `w`eighted `c`ross `v`alidation.

- `y1`: the 1st training set
- `y2`: the 2nd training set
- `yv`: the validation set
- `wt`: weight for training values
- `ct`: correlation between the traits
- `Gi`: inverse of GRM
- `X`: fixed effects matrix
"""
function mtwcv(y1, y2, yv, wt, ct, Gi, X, grp)
    XR = X .* wt
    yt = [y1; y2]

    n1 = length(y1)
    n2 = length(y2)
    nt = n1 + n2
    nv = length(yv)
    nn = nt + nv
    nx = size(X)[2]

    Z = sparse(zeros(nt, 2nn))
    Z[1:n1, 1:n1] = I(n1)
    Z[n1+1:n1+n2, nn+n1+1:nn+n1+n2] = I(n2)
    ZR = Z .* wt

    T = inv([1 ct; ct 1])
    lhs = [ XR'X XR'Z
            ZR'X ZR'Z + kron(T, Gi) ]
    rhs = [ XR'yt
            ZR'yt ]
    
    ebv = lhs \ rhs             # all ebv
    fra = (grp == 1) ?
        nx + n1 + n2 :
        nx + n1 + n2 + nn

    vbv = ebv[fra+1:fra+nv]     # ebv for the validation set
    cc = cor(yv, vbv)
    return cc, vbv
end

function mtcv(y1, y2, yv, ct, Gi, X, grp)
    yt = [y1; y2]

    n1 = length(y1)
    n2 = length(y2)
    nt = n1 + n2
    nv = length(yv)
    nn = nt + nv
    nx = size(X)[2]

    Z = sparse(zeros(nt, 2nn))
    Z[1:n1, 1:n1] = I(n1)
    Z[n1+1:n1+n2, nn+n1+1:nn+n1+n2] = I(n2)

    T = inv([1 ct; ct 1])
    lhs = [ X'X X'Z
            Z'X Z'Z + kron(T, Gi) ]
    rhs = [ X'yt
            Z'yt ]
    
    ebv = lhs \ rhs             # all ebv
    fra = (grp == 1) ?
        nx + n1 + n2 :
        nx + n1 + n2 + nn
    vbv = ebv[fra+1:fra+nv]     # ebv for the validation set
    cc = cor(yv, vbv)
    return cc, vbv
end


function cv_G(d1, d2, ct, vs, G, trait, erc, X; grp = 1)
    y1 = d1[:, trait]
    y2 = d2[:, trait]
    
    w1 = d1[:, erc]
    w2 = d2[:, erc]
    wt = [w1; w2]
    
    yv = vs[:, trait]
    
    i1 = d1[:, :ix]
    i2 = d2[:, :ix]
    i3 = vs[:, :ix]
    id = [i1; i2; i3]

    gi = inv(G[id, id])
    
    c1, _ = mtwcv(y1, y2, yv, wt, ct, gi, X, grp)
    c2, _ = mtcv(y1, y2, yv, ct, gi, X, grp)
    return c1, c2
end

function cv_G(ts, vs, G, trait, erc, X)
    yt = ts[:, trait]
    wt = ts[:, erc]

    yv = vs[:, trait]
    
    it = ts[:, :ix]
    iv = vs[:, :ix]
    id = [it; iv]
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
    ct = 0.9
    println("correlation = $ct")

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
    rst[1, 5], rst[4, 5] = cv_G(dts, nts, ct, dvs, G, trait, erc, X)
    # G + N → D
    X = make_X(nrow(gts), nrow(nts))
    rst[1, 6], rst[4, 6] = cv_G(gts, nts, ct, dvs, G, trait, erc, X)
    # D + G → G
    X = make_X(nrow(dts), nrow(gts))
    rst[2, 4], rst[5, 4] = cv_G(dgts, gvs, G, trait, erc, X)
    # D + N → G
    X = make_X(nrow(dts), nrow(nts))
    rst[2, 5], rst[5, 5] = cv_G(dts, nts, ct, gvs, G, trait, erc, X)
    # G + N → G
    X = make_X(nrow(gts), nrow(nts))
    rst[2, 6], rst[5, 6] = cv_G(gts, nts, ct, gvs, G, trait, erc, X)
    # D + G → N
    X = make_X(nrow(dts), nrow(gts))
    rst[3, 4], rst[6, 4] = cv_G(dts, gts, ct, nvs, G, trait, erc, X, grp = 2)
    # N + D → N
    X = make_X(nrow(nts), nrow(dts))
    rst[3, 5], rst[6, 5] = cv_G(nts, dts, ct, nvs, G, trait, erc, X)
    # N + G → N
    X = make_X(nrow(nts), nrow(gts))
    rst[3, 6], rst[6, 6] = cv_G(nts, gts, ct, nvs, G, trait, erc, X)

    @info "Prediction of one country with all 3 coutries"
    #dgnts = vcat(dts, gts, nts)
    # D + G + N → D
    X = make_X(nrow(dts), nrow(gts), nrow(nts))
    rst[1, 7], rst[4, 7] = cv_G(dgts, nts, ct, dvs, G, trait, erc, X)
    # D + G + N → G
    X = make_X(nrow(dts), nrow(gts), nrow(nts))
    rst[2, 7], rst[5, 7] = cv_G(dgts, nts, ct, gvs, G, trait, erc, X)
    # D + G + N → N
    X = make_X(nrow(dts), nrow(gts), nrow(nts))
    rst[3, 7], rst[6, 7] = cv_G(dgts, nts, ct, nvs, G, trait, erc, X, grp = 2)

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

function cv_210621()
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

#=
proto type from Theo.

function mt_am(y, RI, ifix, itr, id, GRMI, GI)
    # y = n*1 vector of records across traits
    # RI = n*1 vector of weights of records (inverse of residual variances)
    # ifix = n*1 indicates fixed effect class variables across traits (i.e. if ifix=1 records belong to the same fixed effect class irrespective trait)
    # itr = n*1 indicates trait of record
    # id = n*1 indicates ID of animal (agrees with row/column in GRM)
    # GRMI = q*q inverse of matrix of genomic relationships
    # GI = ntrait*ntrait inverse of matrix of genetic (co)variances of the traits
    ntrait = size(GI, 1)
    nanim = size(GRMI, 1)
    n = size(y, 1)
    X = zeros(n, maximum(ifix))
    XR = zero(n, maximum(ffix))
    for i = 1:n
        X[i,ifix[i]] = 1
        XR[i,ifix[i]] = RI[i]
    end
    Z = zeros(n, nanim * ntrait)
    ZR = zeros(n, nanim*ntrait)
    for i = 1:n
        iid = (itr[i] - 1) * nanim + id[i]
        Z[i, iid] = 1
        ZR[i, iid] = RI[i]
    end
    LHS = [XR'X XR'Z
           ZR'X ZR'Z + kron(GI, GRMI)]
    RHS = [XR'y
           ZR'y]
    EBV = LHS\RHS
    return EBV #size(EBV) = nanim*ntrait; EBV' = [ebv_tr1 ebv_tr2 ebv_tr3 ....] 
end
=#
