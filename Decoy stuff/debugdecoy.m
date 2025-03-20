for i = 1:1:length(rho)
    d = eig(rho{i});
    if all(d>=-eps) && isequal(rho{i},rho{i}')
        continue;
    else
        i
    end
end

for i = 1:1:length(sigma)
    d = eig(sigma{i});
    if all(d>=-eps) && isequal(sigma{i},sigma{i}')
        continue;
    else
        i
    end
end