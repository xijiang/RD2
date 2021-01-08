function compare(yt, wt, yv, Gi, Ai) # using weighted MME
    nt = length(yt)
    nv = length(yv)

    Ri = diagm(wt)
    X = ones(nt)
    Z = sparse([I zeros(nt, nv)])
    lhs = [ X'Ri*X X'Ri*Z
            Z'Ri*X Z'Ri*Z + Gi ] # as va = ve = 1, h2=0.5.
    rhs = [ X'Ri*yt
            Z'Ri*yt ]
    bvg = lhs \ rhs
    c1 = cor(yv, bvg[nt+2:end])

    lhs = [ X'Ri*X X'Ri*Z
            Z'Ri*X Z'Ri*Z + Ai ] # as va = ve = 1, h2=0.5.
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
