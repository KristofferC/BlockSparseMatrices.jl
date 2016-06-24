immutable SparseMatrixBSC{Tv, Ti <: Integer}  <: AbstractBlockMatrix{Tv, SparseMatrixBSC}
    R::Int                 # Block size in rows
    C::Int                 # Block size in columns
    m::Int                 # Number of rows
    n::Int                 # Number of columns
    colptr::Vector{Ti}     # Column i is in colptr[i]:(colptr[i+1]-1)
    rowval::Vector{Ti}     # Row values of blocks
    nzval::Array{Tv, 3}    # Nonzero values, one "matrix" per block, nzval[i, j, block]

  function SparseMatrixBSC(R::Integer, C::Integer,
                           m::Integer, n::Integer,
                           colptr::Vector{Ti}, rowval::Vector{Ti},  nzval::Array{Tv, 3})
      m < 0 && throw(ArgumentError("number of rows (m) must be ≥ 0, got $m"))
      n < 0 && throw(ArgumentError("number of columns (n) must be ≥ 0, got $n"))

      R < 0 && throw(ArgumentError("block size y must be ≥ 0, got $x_block_size"))
      C < 0 && throw(ArgumentError("block size x must be ≥ 0, got $y_block_size"))

      m % R != 0 && throw(ArgumentError("row block size: $(R) must evenly divide number of rows: $m"))
      n % C != 0 && throw(ArgumentError("column block size: $(C) must evenly divide number of rows: $n"))
      new(Int(R), Int(C), Int(m), Int(n), colptr, rowval, nzval)
  end
end

function SparseMatrixBSC{Tv, Ti <: Integer}(m::Integer, n::Integer, colptr::Vector{Ti}, rowval::Vector{Ti}, nzval::Array{Tv, 3})
    R, C = size(nzval, 1), size(nzval, 2)
    SparseMatrixBSC{Tv, Ti}(R, C, m, n, colptr, rowval, nzval)
end

#####################
# Utility functions #
#####################

nblocks(A::SparseMatrixBSC) = (length(A.colptr) - 1, A.n ÷ A.C)
nblocks(A::SparseMatrixBSC, i::Int) = nblocks(A)[i]

blocksize(A::SparseMatrixBSC) = A.R, A.C
blocksize(A::SparseMatrixBSC, i::Int) = blocksize(A)[i]

nnzblocks(A::SparseMatrixBSC) = size(A.nzval, 1) * size(A.nzval, 2)
nzblockrange(A::SparseMatrixBSC, col::Integer) =  Int(A.colptr[col]):Int(A.colptr[col + 1] - 1)

Base.LinearIndexing(::Type{SparseMatrixBSC}) = Base.LinearSlow()
Base.size(A::SparseMatrixBSC) = (A.m, A.n)
Base.nnz(A::SparseMatrixBSC) = length(A.nzval)

# Computes the blockindex which is the block I, J where the element with global index i, j
# would occupy
@inline blockindex(A::SparseMatrixBSC, i::Integer, j::Integer) = (i - 1) ÷ A.R + 1, (j - 1) ÷ A.C + 1

# Given a blockindex I, J and a global index i, j computes the offset α, β into the block such that
# A[i, j] = A.nzval[α, β, blockindex(I, J)]
@inline function blockoffsets(A::SparseMatrixBSC, I::Integer, J::Integer, i::Integer, j::Integer)
     α = i - A.R * (I - 1)
     β = j - A.C * (J - 1)
     return α, β
end

# Transforms from a global index i, j to a `BlockIndex`, see the BlockArrays package
function global2blockindex(A::SparseMatrixBSC, i::Integer, j::Integer)
    I, J = blockindex(A, i, j)
    α, β = blockoffsets(A, I, J, i, j)
    return BlockIndex((I, J), (α, β))
end

# Returns -1 if no block found else returns the block A[., ., block]
function findblock(block_array::SparseMatrixBSC, block_index::BlockIndex{2})
  nzrange = nzblockrange(A, block_index.I[2])
  # No block stored in this column
  if start(nzrange) > last(nzrange)
    return -1
  end

  block = searchsortedlast(A.rowval, block_index.I[1], start(nzrange), last(nzrange), Forward)
  if start(nzrange) > block
      return -1
  end

  if rowinblock(A, A.rowval[block], i)
      return block
  else
      return -1
  end
end

function Base.getindex(block_array::SparseMatrixBSC, block_index::BlockIndex{2})
  block = findblock(block_array, block_index)
  if block == -1
      return zero(T)
  else
    return A.nzval[block_index.α[1], block_index.α[2], block]
  end
end



 function Base.setindex!{T,N}(block_array::BlockArray{T, N}, v, block_index::BlockIndex{N})
      getblock(block_array, block_index.I...)[block_index.α...] = v
  end

function getblock{T}(A::SparseMatrixBSC{T}, I::Integer, J::Integer)
    @boundscheck checkbounds(A, I * blocksize(A, 1), J * blocksize(A, 2))
    getblock!(zeros(T, blocksize(A)), A, I, J)
end

function getblock!{T}(block::Matrix{T}, A::SparseMatrixBSC{T}, I::Integer, J::Integer)
    @boundscheck checkbounds(A, I * blocksize(A, 1), J * blocksize(A, 2))
    @boundscheck @assert size(block) == blocksize(A, 1)

    nzrange = nzblockrange(A, J)

    if start(nzrange) > last(nzrange)
       return fill!(block, 0.0)
    end
    block = searchsortedlast(A.rowval, I, start(nzrange), last(nzrange), Forward)

    if (start(nzrange) > block) || (A.rowval[block] != I)
         return fill!(block, 0.0)
    else
        @inbounds for i in blocksize(A, 1), j in blocksize(A, 2)
            block[i, j] = A.nzval[i, j, block]
        end
        return block
    end
end


############
# Printing #
############
function Base.showarray(io::IO, A::SparseMatrixBSC;
                   header::Bool=true, repr=false)
    print(io, A.m, "×", A.n, " sparse block matrix with ", nnzblocks(A), " ", eltype(A), " block entries of size ",
        blocksize(A, 1),"×",blocksize(A, 2))
end

function Base.getindex{T}(A::SparseMatrixBSC{T}, i::Integer, j::Integer)
    @boundscheck checkbounds(A, i, j)
    J = blockcol(A, j)
    I = blockrow(A, i)
    nzrange = nzblockrange(A, J)
    # No block stored in this column
    if start(nzrange) > last(nzrange)
      return zero(T)
    end

    block = searchsortedlast(A.rowval, I, start(nzrange), last(nzrange), Forward)
    if start(nzrange) > block
        return zero(T)
    end

    @inbounds if rowinblock(A, A.rowval[block], i)
        i_blk, j_blk = blockoffsets(A, I, J, i, j)
        return A.nzval[i_blk, j_blk, block]
    else
        return zero(T)
    end
end


###########################
# Conversions CSC <-> BSC #
###########################

# Strategy:
function SparseMatrixCSC{Tv, Ti <: Integer}(A::SparseMatrixBSC{Tv, Ti})
    if blocksize(A) == (1,1)
        return SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, vec(A.nzval))
    end
    rowval = zeros(Ti, nnz(A))
    colptr = zeros(Ti, size(A, 2) + 1)
    nzval = zeros(Tv, length(A.nzval))

    count_row = 1
    count_col = 2
    colptr[1] = 1
    @inbounds for col in 1:nblocks(A, 2)
        blockrange = nzblockrange(A, col)
        n_blocks_col = length(blockrange)
        nnz_values_col = n_blocks_col * blocksize(A, 1)
        for j_blk in 1:blocksize(A, 2)
            # The new colptr is the previous one plus the number of nonzero elements in this column.
            colptr[count_col] = colptr[count_col - 1] + nnz_values_col
            count_col += 1
            for block in blockrange
                i_offset = (A.rowval[block] - 1) * blocksize(A, 1)
                for i_blk in 1:blocksize(A, 1)
                    nzval[count_row] = A.nzval[i_blk, j_blk, block]
                    rowval[count_row] = i_blk + i_offset
                    count_row += 1
                end
            end
        end
    end
    return SparseMatrixCSC(A.m, A.n, colptr, rowval, nzval)
end


function SparseMatrixBSC{Tv, Ti <: Integer}(A::SparseMatrixCSC{Tv, Ti}, R::Integer, C::Integer)
    if (R, C) == (1,1)
        return SparseMatrixBSC(1, 1, A.m, A.n, A.colptr, A.rowval, reshape(A.nzval, length(A.nzval), 1, 1))
    end

    A.m % R != 0 && throw(ArgumentError("row block size: $(R) must evenly divide number of rows: $(A.m)"))
    A.n % C != 0 && throw(ArgumentError("column block size: $(C) must evenly divide number of rows: $(A.n)"))

    Anzval = A.nzval
    Arowval = A.rowval
    Acolptr = A.colptr

    # Upper bound of number of nonzero blocks is nnz(A).
    rowval = zeros(Ti, nnz(A))
    colptr = zeros(Ti, A.n ÷ C + 1)

    n_colblocks = div(A.n, C)

    colptr[1] = 1
    row_counter = 1
    rows = Int[]

    # Strategy to compute rowval and colptr:
    # For each column block accumulate all the rowvalues in the CSC matrix for that block.
    # Convert these to what rowblock they represent
    # Each unique rowblock should now be entered in order into rowval.
    # Colptr for this column block is incremented by the number of unique rowblock values.
    @inbounds for colblock in 1:n_colblocks
        j_offset = (colblock - 1) * C
        row_block_counter = 1

        # Count the number of non zero values for the columns in this column block
        nzvals_block = Acolptr[j_offset + C + 1] - Acolptr[j_offset + 1]
        if nzvals_block == 0 # No nz values in this column block, exit early
            colptr[colblock + 1] = row_counter
            continue
        end

        # Accumulate rowvals for this block
        resize!(rows, nzvals_block)
        for j_blk in 1:C
            col = j_offset + j_blk
            nz_range = Acolptr[col]:Acolptr[col + 1] - 1
            for r in nz_range
                rows[row_block_counter] = Arowval[r]
                row_block_counter += 1
            end
        end

        # Convert from row values -> block rows
        @simd for i in 1:length(rows)
            rows[i] = (rows[i] - 1) ÷ R + 1
        end

        # Pick out the unique block rows and put them into rowval
        sort!(rows) # <- A bit of a bottle enck, takes about 30% of tot time.
        rowval[row_counter] = rows[1] # We have at least one value in rows so this is ok
        row_counter += 1
        for i in 2:length(rows)
            if rows[i] > rows[i-1]
                rowval[row_counter] = rows[i]
                row_counter += 1
            end
        end
        colptr[colblock + 1] = row_counter
    end

    # We now know the true number of non zero blocks so we reshape rowval
    # and allocate the exact space we need for nzval
    deleteat!(rowval, row_counter:length(rowval))
    nzval = zeros(Tv, R, C, length(rowval))


    @inbounds for colblock in 1:n_colblocks
        j_offset = (colblock - 1) * C
        for j_blk in 1:C
            current_block = colptr[colblock]
            col = j_offset + j_blk
            for r in Acolptr[col]:Acolptr[col + 1] - 1
                row = Arowval[r]
                # Looking for the correct block for this column
                while row > rowval[current_block] * R
                    current_block += 1
                end
                i_blk = row - (rowval[current_block] - 1) * R
                nzval[i_blk, j_blk, current_block] = Anzval[r]
            end
        end
    end

    SparseMatrixBSC(A.m, A.n, colptr, rowval, nzval)
end



##########
# LinAlg #
##########

# We do dynamic dispatch here so that the size of the blocks are known at compile time
function A_mul_B!{Tv}(b::Vector{Tv}, A::SparseMatrixBSC{Tv}, x::Vector{Tv})
    if blocksize(A) == (1,1)
        A_mul_B!(b, SparseMatrixCSC(A), x)
    else
        _A_mul_B!(b, A, x, Val{A.C}, Val{A.R})
    end
end

function _A_mul_B!{Tv, C, R}(b::Vector{Tv}, A::SparseMatrixBSC{Tv}, x::Vector{Tv},
                             ::Type{Val{C}}, ::Type{Val{R}})
    fill!(b, 0.0)
    n_cb, n_rb = nblocks(A)
    for J in 1:n_cb
        j_offset = (J - 1)  * blocksize(A, 2)
        for r in nzblockrange(A, J)
            @inbounds i_offset = (A.rowval[r] - 1) * blocksize(A, 1)
            matvec_kernel!(b, A, x, r, i_offset, j_offset, Val{C}, Val{R})
        end
    end
    return b
end

# TODO: Possibly SIMD.jl and do a bunch of simd magic coolness for small mat * vec
@inline function matvec_kernel!{C, R, T}(b::Vector{T}, A::SparseMatrixBSC{T}, x::Vector{T}, r,
                                              i_offset, j_offset, ::Type{Val{C}}, ::Type{Val{R}})
    @inbounds for j in 1:R
        for i in 1:C
            b[i_offset + i] += A.nzval[i, j, r] * x[j_offset + j]
        end
    end
end

function Base.:*{Tv}(A::SparseMatrixBSC{Tv}, x::Vector{Tv})
    A_mul_B!(similar(x, A.m), A, x)
end
