function swap = swap(dimA, dimB)

    testbasis = 0;
    basis = zeros(dimA*dimB);
    for i = 1:1:dimA
        for j = 1:1:dimB
            basis((j-1)*dimA + i, (i-1)*dimB + j) = 1;
            testbasis = testbasis + kron(zket(dimB,j)*zket(dimA,i)',zket(dimA,i)*zket(dimB,j)');
        end
    end
    flag = isequal(testbasis-basis, zeros(dimA*dimB));
    swap = basis;
end