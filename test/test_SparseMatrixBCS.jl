srand(1)
@time for i = 1:8
    n = prod(1:i)
    a = sprand(n, n, 0.1 / n)
    for j in 1:i
        @test sumabs2(SparseMatrixCSC(SparseMatrixBSC(a, j, i)) - a) == 0
        v = rand(size(a, 2))
        @test SparseMatrixBSC(a, j, i) * v == a * v
    end
end
