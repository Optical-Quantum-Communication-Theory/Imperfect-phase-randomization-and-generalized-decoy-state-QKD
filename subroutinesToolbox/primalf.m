%%  FUNCTION NAME: primalObjective
% This file contains the primal problem objective function.
%
% $f(\rho) := D(\mathcal{G}(\rho)||\mathcal{Z}(\mathcal{G}(\rho)))$
%
% Syntax:  fval = primalf(rho,keyMap,krausOperators)
%
% Input: 
%
%  * rho  - density matrix shared by Alice and Bob
% 
%  * keyMap - Alice's key map PVM (If Alice's key map POVM is not projective, use Naimark's extension)
%
%  * krausOperators - The Kraus operators for the post-selection map of
%  Alice and Bob.
%
% Output:
%  
%  * fval - the objective function value. 
%%

function fval = primalf(rho,keyMap,krausOperators)

% if there is no Kraus operator, then proceed the calculation without the G
% map.
if nargin == 2 || isempty(krausOperators)
    
    eigMin = lambda_min(rho);  % check the minimum eigenvalue of this density matrix
    dim = size(rho,1); % get the dimension of the density matrix.
    if eigMin <= 0
        % if the eigenvalue is less than one, do a perturbation by using a
        % pinching channel.
        epsilon = (1e-14-eigMin)*dim; % choose the epsilon value for perturbation
        rho = (1-epsilon)*rho + epsilon*eye(dim)/dim;
    end
    
    zRho = 0;
    for jMapElement = 1:numel(keyMap)
        zRho = zRho + keyMap{jMapElement}*rho*keyMap{jMapElement};
    end % calculate the Z(\rho)
    
    eigMin = lambda_min(zRho); % check the minimum eigenvalue of Z(G(\rho))
    if eigMin <= 0
       epsilon = (1e-14-eigMin)*dim;
        zRho = (1-epsilon)*zRho + epsilon*eye(dim)/dim;
    end

    fval = real(trace(rho*(logm(rho)-logm(zRho)))); % calculate the quantum relative entropy
else
    % for the case there is a post-selection map.
    
    gRho = krausFunc(rho,krausOperators); % calculate G(\rho).
    dim = size(gRho,1); % get the dimension of G(\rho).
    eigMin = lambda_min(gRho); % check the minimum eigenvalue of this density matrix
    if eigMin <= 0
       % if the eigenvalue is less than one, do a perturbation by using a
        % pinching channel.
       epsilon = (1e-14-eigMin)*dim;
       gRho = (1-epsilon)*gRho + epsilon*eye(dim)/dim;
    end

    zRho = 0;
    for jMapElement = 1:numel(keyMap)
        zRho = zRho + keyMap{jMapElement}*gRho*keyMap{jMapElement};
    end %calculate the Z(G(\rho))
    
    eigMin = lambda_min(zRho); % check the minimum eigenvalue of Z(G(\rho))
    if eigMin <= 0
       epsilon = (1e-14-eigMin)*dim;
        zRho = (1-epsilon)*zRho + epsilon*eye(dim)/dim;
    end
    fval = real(trace(gRho*(logm(gRho)-logm(zRho)))); % calculate the quantum relative entropy
end

end