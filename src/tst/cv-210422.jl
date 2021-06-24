function mk_cv_tbl(dts, dvs, gts, gvs, nts, nvs, G, trait, erc)
    # target is each country in turn
    rst = zeros(6, 7)           # 6 = w/nw by (DGN); 7 = D G N DG DN GN DGN

    @info "Prediction of one country with one of 3 coutries"
    # D → D
    rst[1, 1], rst[4, 1] = cv_G(dts, dvs, G, trait, erc)
    # G → D
    rst[1, 2], rst[4, 2] = cv_G(gts, dvs, G, trait, erc)
    # N → D
    rst[1, 3], rst[4, 3] = cv_G(nts, dvs, G, trait, erc)
    # D → G
    rst[2, 1], rst[5, 1] = cv_G(dts, gvs, G, trait, erc)
    # G → G
    rst[2, 2], rst[5, 2] = cv_G(gts, gvs, G, trait, erc)
    # N → G
    rst[2, 3], rst[5, 3] = cv_G(nts, gvs, G, trait, erc)
    # D → G
    rst[3, 1], rst[6, 1] = cv_G(dts, nvs, G, trait, erc)
    # G → G
    rst[3, 2], rst[6, 2] = cv_G(gts, nvs, G, trait, erc)
    # N → G
    rst[3, 3], rst[6, 3] = cv_G(nts, nvs, G, trait, erc)

    @info "Prediction of one country with two of 3 coutries"
    dgts = vcat(dts, gts)
    dnts = vcat(dts, nts)
    gnts = vcat(gts, nts)
    # D + G → D
    rst[1, 4], rst[4, 4] = cv_G(dgts, dvs, G, trait, erc)
    # D + N → D
    rst[1, 5], rst[4, 5] = cv_G(dnts, dvs, G, trait, erc)
    # G + N → D
    rst[1, 6], rst[4, 6] = cv_G(gnts, dvs, G, trait, erc)
    # D + G → G
    rst[2, 4], rst[5, 4] = cv_G(dgts, gvs, G, trait, erc)
    # D + N → G
    rst[2, 5], rst[5, 5] = cv_G(dnts, gvs, G, trait, erc)
    # G + N → G
    rst[2, 6], rst[5, 6] = cv_G(gnts, gvs, G, trait, erc)
    # D + G → N
    rst[3, 4], rst[6, 4] = cv_G(dgts, nvs, G, trait, erc)
    # D + N → N
    rst[3, 5], rst[6, 5] = cv_G(dnts, nvs, G, trait, erc)
    # G + N → N
    rst[3, 6], rst[6, 6] = cv_G(gnts, nvs, G, trait, erc)
    @info "Prediction of one country with all 3 coutries"
    dgnts = vcat(dts, gts, nts)
    # D + G + N → D
    rst[1, 7], rst[4, 7] = cv_G(dgnts, dvs, G, trait, erc)
    # D + G + N → G
    rst[2, 7], rst[5, 7] = cv_G(dgnts, gvs, G, trait, erc)
    # D + G + N → N
    rst[3, 7], rst[6, 7] = cv_G(dgnts, nvs, G, trait, erc)

    round.(rst, sigdigits = 2)
end

"""
Test all the combinations to predict a 3rd country with 
- 3 x 1→1
- 3 x 2→2
- 1 x 3→1
## DRP
- drp-training-2021-04-08.jld: with Norwegian DRP updated
- drp-training-2021-04-21.jld: with german data updated
"""
function cv_21_04_22(drp)
    @load "$dat_dir/jld/cv-setup.jld" dcs gcs ncs
    @load "$dat_dir/jld/$drp"         dts gts nts
    @load "$dat_dir/jld/ebv-all.jld"  dbv gbv nbv
    @load "$dat_dir/jld/GAID.jld"     G A ID

    begin                       # The milk trait
        dt, dv = mk_milk_dt(dts, dcs, dbv)
        gt, gv = mk_milk_dt(gts, gcs, gbv)
        nt, nv = mk_milk_dt(nts, ncs, nbv)
        rst = mk_cv_tbl(dt, dv, gt, gv, nt, nv, G, :Milk, :Milk_ERC)
        for row in eachrow(rst)
            println("| ", join(row, " | "), " |")
        end
    end

    begin                       # the SCS trait
        dt, dv = mk_scs_dt(dts, dcs, dbv)
        gt, gv = mk_scs_dt(gts, gcs, gbv)
        nt, nv = mk_scs_dt(nts, ncs, nbv)
        rst = mk_cv_tbl(dt, dv, gt, gv, nt, nv, G, :SCS, :SCS_ERC)
        for row in eachrow(rst)
            println("| ", join(row, " | "), " |")
        end
    end
end

"""
    function cv_21_04_22_2(drp)
---
Include the validation set also
"""
function cv_21_04_22_2(drp)
    @load "$dat_dir/jld/cv-setup.jld" dcs gcs ncs
    @load "$dat_dir/jld/$drp"         dts gts nts
    @load "$dat_dir/jld/ebv-all.jld"  dbv gbv nbv
    @load "$dat_dir/jld/GAID.jld"     G A ID
    G += 0.01I

    begin                       # The milk trait
        dt, dv = mk_milk_dt_2(dts, dcs, dbv)
        gt, gv = mk_milk_dt_2(gts, gcs, gbv)
        nt, nv = mk_milk_dt_2(nts, ncs, nbv)
        rst = mk_cv_tbl(dt, dv, gt, gv, nt, nv, G, :Milk, :Milk_ERC)
        for row in eachrow(rst)
            println("| ", join(row, " | "), " |")
        end
    end

    begin                       # the SCS trait
        dt, dv = mk_scs_dt(dts, dcs, dbv)
        gt, gv = mk_scs_dt(gts, gcs, gbv)
        nt, nv = mk_scs_dt(nts, ncs, nbv)
        rst = mk_cv_tbl(dt, dv, gt, gv, nt, nv, G, :SCS, :SCS_ERC)
        for row in eachrow(rst)
            println("| ", join(row, " | "), " |")
        end
    end
end

function mk_milk_dt_2(ts, cs, bv)
    df1 = select(innerjoin(cs, ts, on = :ID), :ix, :Milk, :Milk_ERC)
    dft = filter(row -> row.Milk < 999, df1)
    t2  = filter(row -> row.st == 'v', cs)
    dfv = select(innerjoin(t2, bv, on = :ID), :ix, :Milk)
    return dft, dfv
end

function mk_scs_dt_2(ts, cs, bv)
    df1 = select(innerjoin(cs, ts, on = :ID), :ix, :SCS, :SCS_ERC)
    dft = filter(row -> row.SCS < 999, df1)
    t2  = filter(row -> row.st == 'v', cs)
    dfv = select(innerjoin(t2, bv, on = :ID), :ix, :SCS)
    return dft, dfv
end
