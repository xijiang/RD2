using SparseArrays, LinearAlgebra, Statistics, StatsPlots, Random

"""
From: Theodorus Meuwissen
Sent: mandag 28. desember 2020 21:26

Hi Xijiang

It seems the DRP don’t work well.  Even when using EBV, the use of weights does
not work well.  Please try to use deregressed EBV Garrick et al.  (Genetics
Selection Evolution 2009 41:55)

Basically:

DRP = EBV/R2

DRP = deregressed proof which is to be analysed; R2 = reliability of the EBV = EDC / (EDC+lambda)

EDC = effective daughter contributions

Lambda = (4-h2)/h2 if EDC reflects number of daughters; h2 is heritability being used.

If EDC reflects effective number of records; lamda is from animal model.

The weight of the deregressed proofs is: lambda*R2/(1-R2)

Correction which does not affect relative weights but actual weights so that they are on the correct scale.
"""
function dr_ebv()
    
end
