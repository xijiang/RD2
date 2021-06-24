
"""
    function cv_milk(ts, cs, bv, G, A, ID)
---
Using data 2021-03-11 for cross validation again.
"""
function cv_milk(ts, vs, G, A)
    id = ts[:, :ix]
    yt = ts[:, :Milk]
    wt = ts[:, :Milk_ERC]
    nt = length(id)

    yv = vs[:, :Milk]
    append!(id, vs.ix)
    nv = length(yv)

    gi = inv(G[id, id])
    ai = inv(A[id, id])
    _, _ = compare(yt, yv, gi, ai)
    _, _ = compare(yt, wt, yv, gi, ai)
end

function cv_2021_03_11()
    @load "$dat_dir/jld/cv-setup.jld"     dcs gcs ncs
    @load "$dat_dir/jld/drp-training.jld" dts gts nts
    @load "$dat_dir/jld/ebv-all.jld"      dbv gbv nbv
    @load "$dat_dir/jld/GAID.jld"         G A ID
    ts, vs = mk_milk_dt(dts, dcs, dbv)
    cv_milk(ts, vs, G, A);
    ts, vs = mk_milk_dt(gts, gcs, gbv)
    cv_milk(ts, vs, G, A);
    ts, vs = mk_milk_dt(nts, ncs, nbv)
    println(size(ts), ' ', size(vs))
    cv_milk(ts, vs, G, A);
end

function cv_2020_07_02()
    @load "$dat_dir/jld/cv-setup.jld"     dcs gcs ncs
    @load "$dat_dir/jld/drp-training-2020-07-02.jld" dts gts nts
    @load "$dat_dir/jld/ebv-all.jld"      dbv gbv nbv
    @load "$dat_dir/jld/GAID.jld"         G A ID
    ts, vs = mk_milk_dt(dts, dcs, dbv)
    cv_milk(ts, vs, G, A);
    ts, vs = mk_milk_dt(gts, gcs, gbv)
    cv_milk(ts, vs, G, A);
    ts, vs = mk_milk_dt(nts, ncs, nbv)
    cv_milk(ts, vs, G, A);
end

"""
This function is important.
It found that I mistakenly included the missing data,
i.e., 999.0, for the calculation,
which led to very wrong results.
Below line was added after.
```julia
    dft = filter(row -> row.Milk<990, df1)
```
"""
function compare_data_2021_3()
    @load "$dat_dir/jld/cv-setup.jld"     dcs gcs ncs
    @load "$dat_dir/jld/ebv-all.jld"      dbv gbv nbv
    pre = pst = nothing
    #test 1
    begin
        @load "$dat_dir/jld/drp-training.jld" dts gts nts
        ts, vs = mk_milk_dt(dts, dcs, dbv)
        pre = copy(ts)
        println(first(ts, 10))
        println(first(vs, 10))
        println(describe(ts))
        println(describe(vs))
    end
    #test 2
    begin
        @load "$dat_dir/jld/drp-training-2020-07-02.jld" dts gts nts
        ts, vs = mk_milk_dt(dts, dcs, dbv)
        pst = copy(ts)
        println(first(ts, 10))
        println(first(vs, 10))
        println(describe(ts))
        println(describe(vs))
    end
    #test 3
    pre, pst
end
