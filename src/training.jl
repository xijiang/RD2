################################################################################
# This part is from Calus et al, JDS, 2016, 99:6403-6419
################################################################################

"""
    function DRP()
---
# Some terms
- DYD: Daughter Yield Deviation, VanRaden and Wiggans, 1991
- DRP: DeRegressed proofs
- EDC: Effective Daughter contributions, Fikse and Banos, 2001.

This DRP is kind of matrix deregression.
"""
function DRP(W, A, ng, nu, y, λ; conv=1.)
    # A is the inverse of thee A matrix.
    # R is the weight of DRP
    # l = 1:ng
    # r = ng+1:ng+nu
    # 
    # lhs = [ng          ones(ng)'W   zeros(nu)'
    #        R*ones(ng)  W+A[l, l].*λ A[l, r].*λ
    #        zeros(r, r) A[r, l].*λ   A[r, r].*λ]
    # rhs = [ ones(ng)'W*y
    #         W*y
    #         zeros(nu) ]
    # b = lhs \ rhs
    μ = 0
end

function ERC(λ, Rel)
    λ * Rel / (1 - Rel)
end

function EDC(h², Rel)
    λₛ = (4 - h²) / h²
    λₛ * Rel / (1 - Rel)
end
