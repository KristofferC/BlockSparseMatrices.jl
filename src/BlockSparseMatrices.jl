module BlockSparseMatrices

using BlockArrays
using LinearAlgebra
using SparseArrays

#import BlockArrays: nblocks, getblock, getblock!, setblock!, BlockIndex

import Base: Order.Forward

export SparseMatrixBSC, nnzblocks, nblocks

include("SparseMatrixBSC.jl")

end # module
