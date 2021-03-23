
"""
    Function read_training()
---
## For bulls
I checked the DRP files: the weights in the files are ERC’s (Effective Record
Contributions) and not reliabilities, so that does seem to be correct.

The values do look somewhat funny, because I rescaled all heritabilities to be
0.5 (details are in the document attached, which I sent last year).

With `h² = 0.5, λ =(1 - h²)/h² = (1 - 0.5)/0.5 = 1`

With a reliability of 99.0 (on a 0-100% scale), the ERC becomes: 
`ERC = λ⋅rel/(1-rel) = 1 × 0.99/(1-0.99) = 99.0` ! 
(note that `rel` is the reliability on a 0-1 scale).

With a reliability of 98.0 (on a 0-100% scale), the ERC becomes: 
`ERC = λ⋅rel/(1-rel) = 1×0.98/(1-0.98) = 49.0`

## And for cows:

With a reliability of 63.0 (on a 0-100% scale), the ERC becomes: 
`ERC = λ⋅rel/(1-rel) = 1×0.63/(1-0.63) = 1.7

So, I think we need to revert to simpler deregression (as we discussed last 
time!) which simply is: `DRP = EBV/rel`, after first subtracting the mean EBV.

I have attached files with DRPs computed using this simpler deregression; the 
files are in the same format as before.  The weights are still ERC’s.

## Two additional points:

I think there may have been a mistake in the file with Norwegian DRPs, where I 
linked IDs to the wrong data.  That should be fixed now.

The difference between the Dutch bulls and cows (at least this is what I think 
we saw from one of Xijiang’s plots), may be due to a different base pertaining 
to the cow and bull EBVs.  So I wonder if I should do the deregression 
separately for cows and bulls.  Note that with the deregression the mean of the 
DRP is removed.  Doing it separately for cows and bulls would effectively give 
both the same mean. 

Anyway, I’m curious to see the results with these new DRPs.

-- Mario, 2021.3
"""
function read_training()
    trn_dir = joinpath(dat_dir, "mario/20210311")
    
    @info "  - Dutch data training set"
    dp = DataFrame(ID          = String[],  #  1. Interbull ID
                   Milk        = Float64[], #  2. milk yield DRP
                   Milk_ERC    = Float64[], #  3. milk yield ERC
                   Fat         = Float64[], #  4. fat yield DRP
                   Fat_ERC     = Float64[], #  5. fat yield ERC
                   Protein     = Float64[], #  6. protein yield DRP
                   Protein_ERC = Float64[], #  7. protein yield ERC
                   SCS         = Float64[], #  8. somatic cell score DRP
                   SCS_ERC     = Float64[], #  9. somatic cell score ERC
                   x10         = Float64[], # 10. NR56 DRP
                   x11         = Float64[], # 11. NR56 ERC
                   x12         = Float64[], # 12. calving to 1st ins DRP
                   x13         = Float64[], # 13. calving to 1st ins ERC
                   x14         = Float64[], # 14. 1st-Last ins. Cows DRP
                   x15         = Float64[], # 15. 1st-Last ins. Cows ERC
                   x16         = Float64[], # 16. calving interval DRP
                   x17         = Float64[], # 17. calving interval ERC
                   x18         = Int[]      # 18. Numerical ID (as in “all.ped”)
                   )
    for line in eachline(joinpath(trn_dir, "dutch_DRPs.txt"))
        v = split(line)
        id = v[1]
        y = parse.(Float64, v[2:end-1])
        n = parse(Int, v[end])
        push!(dp, [id; y; n])
    end

    @info "  - German data training set"
    gp = DataFrame(ID          = String[],  #  1. Interbull ID
                   Milk        = Float64[], #  2. milk yield DRP
                   Milk_ERC    = Float64[], #  3. milk yield ERC
                   Fat         = Float64[], #  4. fat yield DRP
                   Fat_ERC     = Float64[], #  5. fat yield ERC
                   Protein     = Float64[], #  6. protein yield DRP
                   Protein_ERC = Float64[], #  7. protein yield ERC
                   SCS         = Float64[], #  8. somatic cell score DRP
                   SCS_ERC     = Float64[], #  9. somatic cell score ERC
                   x10         = Float64[], # 10. NR56 heifers DRP
                   x11         = Float64[], # 11. NR56 heifers ERC
                   x12         = Float64[], # 12. calving to 1st ins DRP
                   x13         = Float64[], # 13. calving to 1st ins ERC
                   x14         = Float64[], # 14. NR56 cows DRP
                   x15         = Float64[], # 15. NR56 cows ERC
                   x16         = Float64[], # 16. 1st-Last ins. Cows DRP
                   x17         = Float64[], # 17. 1st-Last ins. Cows ERC
                   x18         = Float64[], # 18. Days open cows DRP
                   x19         = Float64[], # 19. Days open cows ERC
                   N           = Int[]      # 20. Numerical ID (as in “all.ped”)
                   )
    for line in eachline(joinpath(trn_dir, "german_DRPs.txt"))
        v = split(line)
        id = v[1]
        y = parse.(Float64, v[2:end-1])
        n = parse(Int, v[end])
        push!(gp, [id; y; n])
    end

    @info "  - Norwegian data training set"
    np = DataFrame(ID       = String[],
                   Milk     = Float64[],
                   Milk_ERC = Float64[],
                   SCS      = Float64[],
                   SCS_ERC  = Float64[],
                   N        = Int[])
    for line in eachline(joinpath(trn_dir, "norway_DRPs.txt"))
        v = split(line)
        length(v) != 6 && continue # because two lines have only 2 columns in trnn set
        id = v[1]
        y = parse.(Float64, v[2:end-1])
        n = parse(Int, v[end])
        push!(np, [id; y; n])
    end
    return (dp, gp, np)
end

function refresh_drp()
    @info "Read phenotypes, EBV and CV sets"
    dts, gts, nts = read_training() # ts = training set
    @save joinpath(dat_dir, "jld/drp-training.jld") {compress=true} dts gts nts
end
