%%  FUNCTION NAME: decoyBoundsSDP
% This file contains the function to compute decoy bounds without assuming
% that all the states are diagonal by solving an SDP.
% Check OneNote for more detailed notes.
%
% Syntax:  [YL,YU] = decoyBoundsSDP(Gamma, rho, sigma, decoy_expectations)
%
% Input: 
%
% * rho - Cell which contains the projected states. So rho{i} = P rho_i P
%         the ith state.
%
% * sigma - Cell which contains the states whose statistics we would like
%           to estimate. Note that these are assumed to lie within the
%           projected space.
%
% * decoy_expectations - Matrix with all measurement outcomes for different
%                        intensities. decoy_expectations(i,j) is the
%                        detection statistics of the ith measurement for
%                        rho_j.
%
% * Wb - The weight outside the projected subspace where the projection is
%        made on the state Bob receives (dimension reduction on the POVMs).
%
% * a - Measure of the degree to which rho commutes with the projection
%
% Output:
%  
% * YL - Lower bounds on yields
%
% * YU - Upper bounds on yields
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [YL, YU] = decoyBoundsSDP(ConstraintPOVM, ObjPOVM, rho, sigma, decoy_expectations, Wb, a)
    
    linearconstrainttolerance = 10^-10;
    
    dimM = length(ConstraintPOVM{1});
    dimS = length(rho{1});
    dim = dimM*dimS;
    for k = 1:numel(sigma)
        for i =  1:1:numel(ObjPOVM)
            cvx_begin sdp quiet
                cvx_precision default
                cvx_solver Mosek
                variable J(dim,dim) hermitian semidefinite
                minimize real(trace(kron(ObjPOVM{i},sigma{k})*J))
                for j = 1:numel(rho)
                    W = 1-trace(rho{j});
                    if W<= eps
                        W = 0;
                    end
                    epsilon = (1-a)*sqrt(W);
                    for l = 1:numel(ConstraintPOVM)
                        real(trace((kron(ConstraintPOVM{l},rho{j})*J))) <= min(decoy_expectations(l,j)+epsilon+linearconstrainttolerance,1);
                        real(trace((kron(ConstraintPOVM{l},rho{j})*J))) >= max(decoy_expectations(l,j)-2*epsilon-2*W-Wb(j)-linearconstrainttolerance,0);
                    end
                end
                1-W-Wb(j)-epsilon <= real(trace(kron(eye(dimM),rho{j})*J)) <= 1;
                if Wb == 0
                    norm(PartialTrace(J,1,[dimM,dimS])-eye(dimS)) <= linearconstrainttolerance;
                else
                    PartialTrace(J,1,[dimM,dimS]) <= eye(dimS)
                end
            cvx_end
            YL(i,k) = real(trace(kron(ObjPOVM{i},sigma{k})*J));
        
            cvx_begin sdp quiet
                cvx_precision default
                cvx_solver Mosek
                variable J(dim,dim) hermitian semidefinite
                maximize real(trace(kron((ObjPOVM{i}),sigma{k})*J))
                for j = 1:numel(rho)
                    W = 1-trace(rho{j});
                    if W<= eps
                        W = 0;
                    end
                    epsilon = (1-a)*sqrt(W);
                    for l = 1:numel(ConstraintPOVM)
                        real(trace((kron(ConstraintPOVM{l},rho{j})*J))) <= min(decoy_expectations(l,j)+epsilon+linearconstrainttolerance,1);
                        real(trace((kron(ConstraintPOVM{l},rho{j})*J))) >= max(decoy_expectations(l,j)-2*epsilon-2*W-Wb(j)-linearconstrainttolerance,0);
                    end
                end
                1-W-Wb(j)-epsilon <= real(trace(kron(eye(dimM),rho{j})*J)) <= 1;
                if Wb == 0
                    norm(PartialTrace(J,1,[dimM,dimS])-eye(dimS)) <= linearconstrainttolerance;
                else
                    PartialTrace(J,1,[dimM,dimS]) <= eye(dimS)
                end
            cvx_end
            YU(i,k) = real(trace(kron((ObjPOVM{i}),sigma{k})*J));     
        end
    end 
end
