function temp = GetPhotonCombo(N,pulseNum)
    if pulseNum == 1
        temp = {[N]};
        return;
    end
    if N == 1
        temp = {};
        for i = 1:pulseNum
            comb = zeros(1,pulseNum);
            comb(i) = 1;
            temp{end+1} = comb;
        end
        return
    end
    comb = zeros(1,pulseNum);
    comb(1) = N;
    temp{1} = comb;
    for i = 1:N
        comb(1) = N-i;
        subCombSet = GetPhotonCombo(i,pulseNum-1);
        for j = 1:numel(subCombSet)
            subComb = subCombSet{j};
            for k = 1:(pulseNum-1)
                comb(k+1) = subComb(k);
            end
            temp{end+1} = comb;
        end
    end
end
