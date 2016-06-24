# BlockSparseMatrices

[![Build Status](https://travis-ci.org/KristofferC/BlockSparseMatrices.jl.svg?branch=master)](https://travis-ci.org/KristofferC/BlockSparseMatrices.jl)

```jl
A = sprand(10, 10, 0.5)

block_size_rows = 2
block_size_columns = 2

# Creating block sparse matrix 
B = SparseMatrixBSC(A, block_size_rows, block_size_columns)

# BlockMatrix times vector
B * rand(size(B,1))

# Convert back to CSC
A2 = SparseMatrixCSC(B)

# Check roundtrip
A == A2
```
