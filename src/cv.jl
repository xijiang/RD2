#=
# Function list:
- compare(yt, wt, yv, Gi, Ai) # using weighted MME
- compare(yt, yv, Gi, Ai) # not weighed MME
- cv_n()
- cv_d()
- cv_g()
- function CV_ebv(training, brding_v, cv_setup, G, A, ID)
=#

# @load "cv-setup.jld"     dcs gcs ncs
# @load "drp-training.jld" dts gts nts
# @load "ebv-all.jld"      dbv gbv nbv
# @load "GAID.jld"         G A ID;

function cv_n()
    ############################################################################
    # Norwegian data
    # note: variables end with '_n' are for Norwegian data
    id_n = copy(nts.ID)
    yt_n = copy(nts.Milk)
    wt_n = copy(nts.Milk_ERC)
    tmp  = filter(row -> row.st == 'v', ncs)
    bv_n = select(innerjoin(tmp, nbv, on = :ID), :ID, :Milk, :Milk_R)
    nv_n = length(bv_n.ID)
    append!(id_n, bv_n.ID)
    yv_n = copy(bv_n.Milk)
    ix_n  = zeros(Int, length(id_n))
    for i in 1:length(id_n)
        ix_n[i] = ID[id_n[i]]
    end
    gi_n = inv(G[ix_n, ix_n])
    ai_n = inv(A[ix_n, ix_n])
    bvg_wt_n, bva_wt_n = compare(yt_n, wt_n, yv_n, gi_n, ai_n)
    bvg_nw_n, bva_nw_n = compare(yt_n, yv_n, gi_n, ai_n);

    # Visualization of above data
    ## G vs. A
    tmp = randsubseq(ix_n, .023)  # only take a subset of ix_n
    scatter(vec(A[tmp, tmp]), vec(G[tmp, tmp]),
            leg=false, xlabel="A matrix sub", ylabel="G matrix sub", ms=1)
    savefig("GvA_n.png")

    ## DRP vs EBV, in the training set
    milk_drp_n = select(nts, :ID, r"Milk")
    milk_ebv_n = select(nbv, :ID, r"Milk")
    tmp = innerjoin(milk_drp_n, milk_ebv_n, on=:ID, renamecols = "_drp" => "_ebv")
    @df tmp scatter(:Milk_drp, :Milk_ebv,
                    ms = 1, xlabel = "Milk DRP", ylabel= "Milk_ebv", legend = false)
    savefig("DRPvEBV_n.png")

    ## EBV vs current GEBV, in the validation set
    p1 = scatter(yv_n, bvg_wt_n, leg=false, ms=1, ylabel="weighted")
    p2 = scatter(yv_n, bva_wt_n, leg=false, ms=1)
    p3 = scatter(yv_n, bvg_nw_n, leg=false, ms=1, ylabel="no weight", xlabel = "G")
    p4 = scatter(yv_n, bva_nw_n, leg=false, ms=1, xlabel="A")
    plot(p1, p2, p3, p4, layout=(2,2))
    savefig("EBVvsNEW_n.png")
    # Note x-axis are all EBV, y-axis are EBV from current analysis.
end

function cv_d()
    ################################################################################
    # Norwegian data
    # note: variables end with '_n' are for Norwegian data
    id_d = copy(dts.ID)
    yt_d = copy(dts.Milk)
    wt_d = copy(dts.Milk_ERC)
    tmp  = filter(row -> row.st == 'v', dcs)
    bv_d = select(innerjoin(tmp, dbv, on = :ID), :ID, :Milk, :Milk_R)
    nv_d = length(bv_d.ID)
    append!(id_d, bv_d.ID)
    yv_d = copy(bv_d.Milk)
    ix_d  = zeros(Int, length(id_d))
    for i in 1:length(id_d)
        ix_d[i] = ID[id_d[i]]
    end
    gi_d = inv(G[ix_d, ix_d])
    ai_d = inv(A[ix_d, ix_d])
    bvg_wt_d, bva_wt_d = compare(yt_d, wt_d, yv_d, gi_d, ai_d)
    bvg_nw_d, bva_nw_d = compare(yt_d, yv_d, gi_d, ai_d);

    # Visualization of above data
    ## G vs. A
    tmp = randsubseq(ix_d, .05) # only take a subset of ix_d
    scatter(vec(A[tmp, tmp]), vec(G[tmp, tmp]),
            leg=false, xlabel="A matrix sub", ylabel="G matrix sub", ms=1)
    savefig("GvA_d.png")

    ## DRP vs EBV, in the training set
    milk_drp_d = select(dts, :ID, r"Milk")
    milk_ebv_d = select(dbv, :ID, r"Milk")
    tmp = innerjoin(milk_drp_d, milk_ebv_d, on=:ID, renamecols = "_drp" => "_ebv")
    @df tmp scatter(:Milk_drp, :Milk_ebv,
                    ms = 1, xlabel = "Milk DRP", ylabel= "Milk_ebv", legend = false)
    savefig("DRPvEBV_d.png")

    ## EBV vs current GEBV, in the validation set
    p1 = scatter(yv_d, bvg_wt_d, leg=false, ms=1, ylabel="weighted")
    p2 = scatter(yv_d, bva_wt_d, leg=false, ms=1)
    p3 = scatter(yv_d, bvg_nw_d, leg=false, ms=1, ylabel="no weight", xlabel = "G")
    p4 = scatter(yv_d, bva_nw_d, leg=false, ms=1, xlabel="A")
    plot(p1, p2, p3, p4, layout=(2,2))
    savefig("EBVvsNEW_d.png")
    # Note x-axis are all EBV, y-axis are EBV from current analysis.
end

function cv_g()
    ################################################################################
    # Norwegian data
    # note: variables end with '_n' are for Norwegian data
    id_g = copy(gts.ID)
    yt_g = copy(gts.Milk)
    wt_g = copy(gts.Milk_ERC)
    tmp  = filter(row -> row.st == 'v', gcs)
    bv_g = select(innerjoin(tmp, gbv, on = :ID), :ID, :Milk, :Milk_R)
    nv_g = length(bv_g.ID)
    append!(id_g, bv_g.ID)
    yv_g = copy(bv_g.Milk)
    ix_g  = zeros(Int, length(id_g))
    for i in 1:length(id_g)
        ix_g[i] = ID[id_g[i]]
    end
    gi_g = inv(G[ix_g, ix_g])
    ai_g = inv(A[ix_g, ix_g])
    bvg_wt_g, bva_wt_g = compare(yt_g, wt_g, yv_g, gi_g, ai_g)
    bvg_nw_g, bva_nw_g = compare(yt_g, yv_g, gi_g, ai_g);

    # Visualization of above data
    ## G vs. A
    tmp = randsubseq(ix_g, .2)  # only take a subset of ix_g
    scatter(vec(A[tmp, tmp]), vec(G[tmp, tmp]),
            leg=false, xlabel="A matrix sub", ylabel="G matrix sub", ms=1)
    savefig("GvA_g.png")

    ## DRP vs EBV, in the training set
    milk_grp_g = select(gts, :ID, r"Milk")
    milk_ebv_g = select(gbv, :ID, r"Milk")
    tmp = innerjoin(milk_grp_g, milk_ebv_g, on=:ID, renamecols = "_grp" => "_ebv")
    @df tmp scatter(:Milk_grp, :Milk_ebv,
                    ms = 1, xlabel = "Milk DRP", ylabel= "Milk_ebv", legend = false)
    savefig("DRPvEBV_g.png")

    ## EBV vs current GEBV, in the validation set
    p1 = scatter(yv_g, bvg_wt_g, leg=false, ms=1, ylabel="weighted")
    p2 = scatter(yv_g, bva_wt_g, leg=false, ms=1)
    p3 = scatter(yv_g, bvg_nw_g, leg=false, ms=1, ylabel="no weight", xlabel = "G")
    p4 = scatter(yv_g, bva_nw_g, leg=false, ms=1, xlabel="A")
    plot(p1, p2, p3, p4, layout=(2,2))
    savefig("EBVvsNEW_g.png")
    # Note x-axis are all EBV, y-axis are EBV from current analysis.
end

function CV_ebv(training, brding_v, cv_setup, G, A, ID)
    ############################################################################
    # Norwegian data
    # note: variables end with '_n' are for Norwegian data

    # Training set
    tst = select(innerjoin(select(training, :ID), brding_v, on = :ID),
                 :ID, r"Milk")
    id = tst[:, :ID]
    yt = tst[:, :Milk]
    wt = tst[:, :Milk_R]
    nt = length(id)

    # validation set
    vst = select(innerjoin(select(filter(row -> row.st == 'v', cv_setup), :ID),
                           brding_v, on = :ID),
                 :ID, r"Milk")
    yv = vst[:, :Milk]
    append!(id, vst.ID)
    nv = length(vst.ID)
    ix = zeros(Int, length(id))
    for i in 1:length(id)
        ix[i] = ID[id[i]]
    end

    gi = inv(G[ix, ix])
    ai = inv(A[ix, ix])
    _, _ = compare(yt, wt, yv, gi, ai)
    _, _ = compare(yt, yv, gi, ai)
end

#function hist(ts, cs, df, ID)
#    @load "cv-setup.jld"     dcs gcs ncs
#    @load "drp-training.jld" dts gts nts
#    @load "ebv-all.jld"      dbv gbv nbv
#    @load "GAID.jld"         G A ID;
#    id = select(ts, :ID)
#    tt = select(inner
#end
