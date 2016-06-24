module BlockSparseMatrices

using BlockArrays

import BlockArrays: nblocks, getblock, getblock!, setblock!, BlockIndex

import Base: Order.Forward, A_mul_B!, MultiplicativeInverses.SignedMultiplicativeInverse, SparseMatrixCSC

export SparseMatrixBSC, nnzblocks, nblocks

include("SparseMatrixBSC.jl")

end # module
