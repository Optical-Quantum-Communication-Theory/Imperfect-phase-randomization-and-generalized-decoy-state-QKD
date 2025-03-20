%% FUNCTION NAME: primalDf
% This function calculates the gradient of primal problem objective
% function.
%%
function dfval = primalDf(rho,keyMap,krausOperators)

    if nargin == 2 || isempty(krausOperators)
        % if there is no post-selection map
        eigMin = lambda_min(rho);
        dim = size(rho,1);
        if eigMin <= 0
            epsilon = (1e-14-eigMin)*dim;
            rho = (1-epsilon)*rho + epsilon*eye(dim)/dim;
        end
        zRho = 0;
        for j = 1:numel(keyMap)
            zRho = zRho + keyMap{j}*rho*keyMap{j};
        end
        eigMin = lambda_min(zRho);
        if eigMin <= 0
            epsilon = (1e-14-eigMin)*dim;
            zRho = (1-epsilon)*zRho + epsilon*eye(dim)/dim;
        end
        dfval = transpose(logm(rho))-transpose(logm(zRho));
    else
        % if there is a post-selection map.
        gRho = krausFunc(rho,krausOperators);
        dim = size(gRho,1);
        eigMin = lambda_min(gRho);
        if eigMin <= 0
            epsilon = (1e-14-eigMin)*dim;
            gRho = (1-epsilon)*gRho + epsilon*eye(dim)/dim;
        end
    
        zRho = 0;
        for j = 1:numel(keyMap)
            zRho = zRho + keyMap{j}*gRho*keyMap{j};
        end
        eigMin = lambda_min(zRho);
        if eigMin <= 0
            epsilon = (1e-14-eigMin)*dim;
            zRho = (1-epsilon)*zRho + epsilon*eye(dim)/dim;
        end
        dfval = transpose(krausFunc(logm(gRho)-logm(zRho),krausOperators,'transpose'));
    end

end