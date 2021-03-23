
"""
    function cv_milk(ts, cs, bv, G, A, ID)
---
Using data 2021-03-11 for cross validation again.
"""
function cv_milk(ts, vs, G, A)
    id = ts[:, :ix]
    yt = ts[:, :Milk]
    nt = length(id)

    yv = vs[:, :Milk]
    append!(id, vs.ix)
    nv = length(yv)

    gi = inv(G[id, id])
    ai = inv(A[id, id])
    _, _ = compare(yt, yv, gi, ai);
end
"""
    function mkcvdt(ts, cs, bv)
---
This function returns training set `:ix` and `:Milk` phenotype data in `ts`,
returns `:ix` and `:Milk` EBV in `cs`.
InterBull ID names are not necessary, as `ix` are indices in `G` and `A`.
"""
function mkcvdt(ts, cs, bv)
    t1  = filter(row -> row.st == 't', cs)
    df1 = select(innerjoin(t1, ts, on = :ID), :ix, :Milk)
    dft = filter(row -> row.Milk<10, df1)
    t2  = filter(row -> row.st == 'v', cs)
    dfv = select(innerjoin(t2, bv, on = :ID), :ix, :Milk)
    return dft, dfv
end

function cv_2021_03_11()
    @load "$dat_dir/jld/cv-setup.jld"     dcs gcs ncs
    @load "$dat_dir/jld/drp-training.jld" dts gts nts
    @load "$dat_dir/jld/ebv-all.jld"      dbv gbv nbv
    @load "$dat_dir/jld/GAID.jld"         G A ID
    ts, vs = mkcvdt(dts, dcs, dbv)
    cv_milk(ts, vs, G, A);
    ts, vs = mkcvdt(gts, gcs, gbv)
    cv_milk(ts, vs, G, A);
    ts, vs = mkcvdt(nts, ncs, nbv)
    cv_milk(ts, vs, G, A);
end

function cv_2020_07_02()
    @load "$dat_dir/jld/cv-setup.jld"     dcs gcs ncs
    @load "$dat_dir/jld/drp-training-2020-07-02.jld" dts gts nts
    @load "$dat_dir/jld/ebv-all.jld"      dbv gbv nbv
    @load "$dat_dir/jld/GAID.jld"         G A ID
    ts, vs = mkcvdt(dts, dcs, dbv)
    cv_milk(ts, vs, G, A);
    ts, vs = mkcvdt(gts, gcs, gbv)
    cv_milk(ts, vs, G, A);
    ts, vs = mkcvdt(nts, ncs, nbv)
    cv_milk(ts, vs, G, A);
end

function compare_data_2021_3()
    @load "$dat_dir/jld/cv-setup.jld"     dcs gcs ncs
    @load "$dat_dir/jld/ebv-all.jld"      dbv gbv nbv
    pre = pst = nothing
    #test 1
    begin
        @load "$dat_dir/jld/drp-training.jld" dts gts nts
        ts, vs = mkcvdt(dts, dcs, dbv)
        pre = copy(ts)
        println(first(ts, 10))
        println(first(vs, 10))
        println(describe(ts))
        println(describe(vs))
    end
    #test 2
    begin
        @load "$dat_dir/jld/drp-training-2020-07-02.jld" dts gts nts
        ts, vs = mkcvdt(dts, dcs, dbv)
        pst = copy(ts)
        println(first(ts, 10))
        println(first(vs, 10))
        println(describe(ts))
        println(describe(vs))
    end
    #test 3
    pre, pst
end
