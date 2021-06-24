function compare(yt, wt, yv, Gi, Ai) # using weighted MME
    nt = length(yt)
    nv = length(yv)

    #Ri = diagm(wt)
    XR = X .* wt
    X  = ones(nt)
    Z  = sparse([I zeros(nt, nv)])
    ZR = Z .* wt
    lhs = [ XR'X XR'Z
            ZR'X ZR'Z + Gi ]    # as va = ve = 1, h2=0.5.
    rhs = [ XR'yt
            ZR'yt ]
    bvg = lhs \ rhs
    c1 = cor(yv, bvg[nt+2:end])

    lhs = [ XR'X XR'Z
            ZR'X ZR'Z + Ai ]    # as va = ve = 1, h2=0.5.
    bva = lhs \ rhs
    c2 = cor(yv, bva[nt+2:end])
    println("Weighted\n- using G: $c1;\n- using A: $c2")
    return bvg[nt+2:end], bva[nt+2:end]
end

function compare(yt, yv, Gi, Ai) # not weighed MME
    nt = length(yt)
    nv = length(yv)

    X = ones(nt)
    Z = sparse([I zeros(nt, nv)])
    lhs = [ X'X X'Z
            Z'X Z'Z + Gi ]      # as va = ve = 1, h2=0.5.
    rhs = [ X'yt
            Z'yt ]
    bvg = lhs \ rhs
    c1 = cor(yv, bvg[nt+2:end])

    lhs = [ X'X X'Z
            Z'X Z'Z + Ai ]      # as va = ve = 1, h2=0.5.
    bva = lhs \ rhs
    c2 = cor(yv, bva[nt+2:end])
    println("Not weighted\n- using G: $c1;\n- using A: $c2")
    return bvg[nt+2:end], bva[nt+2:end]
end

"""
    function ccwcv(yt, wt, yv, Gi)
---
This is a results of the meeting on 2021-03-24 for
`c`ross `c`ountry `w`eighted `c`ross `v`alidation.

The conclusion was that to include weigt is a must.
"""
function ccwcv(yt, wt, yv, Gi)
    nt = length(yt)
    nv = length(yv)

    Ri = diagm(wt)
    X = ones(nt)
    Z = sparse([I zeros(nt, nv)])
    lhs = [ X'Ri*X X'Ri*Z
            Z'Ri*X Z'Ri*Z + Gi ] # as va = ve = 1, h2=0.5.
    rhs = [ X'Ri*yt
            Z'Ri*yt ]
    bv = lhs \ rhs
    cc = cor(yv, bv[nt+2:end])
    return cc, bv[nt+2:end]
end

"""
    function cccv(yt, yv, Gi)
---
This is a results of the meeting on 2021-03-24 for
`c`ross `c`ountry `c`ross `v`alidation.

The conclusion was that to include weigt is a must.
"""
function cccv(yt, yv, Gi)
    nt = length(yt)
    nv = length(yv)

    X = ones(nt)
    Z = sparse([I zeros(nt, nv)])
    lhs = [ X'X X'Z
            Z'X Z'Z + Gi ]      # as va = ve = 1, h2=0.5.
    rhs = [ X'yt
            Z'yt ]
    bv = lhs \ rhs
    cc = cor(yv, bv[nt+2:end])
    return cc, bv[nt+2:end]
end

"""
    function make_tv(ts, cs, bv, drp, erc)
---
This function returns training set `:ix` and `:Milk` phenotype data in `ts`,
returns `:ix` and `:Milk` EBV in `cs`.
InterBull ID names, or other real names are not necessary,
as `ix` are indices in `G` and `A`.
"""
function make_tv(ts, cs, bv, drp, erc)
    t1 = filter(row ->row.st == 't', cs)
    df1 = select(innerjoin(t1, ts, on = :ID), :ix, drp, erc)
    dft = filter(row -> row.Milk < 999, df1)
    t2  = filter(row -> row.st == 'v', cs)
    dfv = select(innerjoin(t2, bv, on = :ID), :ix, drp)
    return dft, dfv
end

"""
    function mk_milk_dt(ts, cs, bv)
---
This function returns training set `:ix` and `:Milk` phenotype data in `ts`,
returns `:ix` and `:Milk` EBV in `cs`.
InterBull ID names, or other real names are not necessary,
as `ix` are indices in `G` and `A`.
"""
function mk_milk_dt(ts, cs, bv)
    t1  = filter(row -> row.st == 't', cs)
    df1 = select(innerjoin(t1, ts, on = :ID), :ix, :Milk, :Milk_ERC)
    dft = filter(row -> row.Milk < 999, df1)
    t2  = filter(row -> row.st == 'v', cs)
    dfv = select(innerjoin(t2, bv, on = :ID), :ix, :Milk)
    return dft, dfv
end


"""
    function mk_scs_dt(ts, cs, bv)
---
This function returns training set `:ix`, `:SCS`, `:SCS_ERC` phenotype data in `ts`,
returns `:ix` and `:SCS` EBV in `cs`.
InterBull ID names, or other real names are not necessary,
as `ix` are indices in `G` and `A`.
"""
function mk_scs_dt(ts, cs, bv)
    t1  = filter(row -> row.st == 't', cs)
    df1 = select(innerjoin(t1, ts, on = :ID), :ix, :SCS, :SCS_ERC)
    dft = filter(row -> row.SCS < 999, df1)
    t2  = filter(row -> row.st == 'v', cs)
    dfv = select(innerjoin(t2, bv, on = :ID), :ix, :SCS)
    return dft, dfv
end

"""
    function cv_G(ts, vs, G)
---
Cross validation with only G matrix.
Returns correlation of weighted and not weighted on ERC.
"""
function cv_G(ts, vs, G, trait, erc)
    id = ts[:, :ix]
    yt = ts[:, trait]
    wt = ts[:, erc]
    nt = length(id)

    yv = vs[:, trait]
    append!(id, vs.ix)
    nv = length(yv)

    gi = inv(G[id, id])
    c1, _ = ccwcv(yt, wt, yv, gi)
    c2, _ = cccv(yt, yv, gi)
    return c1, c2
end
