%% FUNCTION NAME:    step2Solver
% Step 2 solver for dimension reduction method, modified for infinite 
% shield system specifically. Refer to Jan 25 notes for details.
%% Copyrights => will become open source
% Author: Twesh Upadhyaya
%
% Created: January 25, 2021
% 
%% Inputs
% Accepts a POVM Gamma, with corresponding expectations gamma, and W

function [lowerbound, debugging, flag] = step2SolverAB(rho, keyMap, sigmaAAs,...
    epsilonsource, Gamma, gammaL, gammaU, W, krausOperators, options)
    
    warning('off','MATLAB:logm:nonPosRealEig');
    defaultoptions.epsilon = 0; % 0<epsilon<=1/(e(d'-1)), epsilon is related to perturbed channel
    defaultoptions.epsilonprime = 1e-12; % epsilonprime is related to constraint tolerance
    
%     if nargin == 8
%         if ~isfield(options,'epsilon')
%             options.epsilon = defaultoptions.epsilon;
%         end
%         if ~isfield(options,'epsilonprime')
%             options.epsilonprime = defaultoptions.epsilonprime;
%         end
%     else 
%         options = defaultoptions;
%     end
    %epsilon = max(0,-lambda_min(rho));
    epsilonprime = options.epsilonprime;
    epsilonprimeprime = options.epsilonprimeprime;
    [fval, epsilon1] = primalfep(rho, keyMap, krausOperators, options);
    [gradf, epsilon2] = primalDfep(rho, keyMap, krausOperators, options); % epsilon1 == epsilon2 should hold

    gradf = (gradf + gradf')/2;
    fval = real(fval);
    if isempty(krausOperators)
        dprime = size(rho, 1);
    else
        dprime = size(krausFunc(rho, krausOperators), 1);
    end
    
    epsilon = max(epsilon1, epsilon2);
    %epsilon = 0;
    if epsilon> 1/(exp(1)*(dprime-1))
      ME = MException('step2Solver:epsilon too large','Theorem cannot be applied. Please have a better rho to start with');
      throw(ME);
    
    end
    debugging = {};
    debugging.lambdamin = lambda_min(rho);
    debugging.epsilon1 = epsilon1;
    debugging.epsilon2 = epsilon2;
    

    Lepsilon = real(fval - trace(rho*gradf.'));

    debugging.Lepsilon = Lepsilon;
    debugging.rho = rho;
    debugging.keyMap = keyMap;
    debugging.krausOp = krausOperators;
    debugging.gradf = gradf;
    debugging.fval = fval;    
        
    dim = size(rho,1);
    dimAAs=size(sigmaAAs,1);
    dimB=dim/dimAAs;
    d=size(gammaL,1);
    
    gamma=[gammaU;-gammaL];        
    
    cvx_begin sdp quiet
%         cvx_precision default
        cvx_solver mosek
        variable y(2*d) nonnegative
        variable y2d1 nonnegative
        variable Yb(dimAAs,dimAAs) hermitian semidefinite
        maximize real(-transpose(y)*(gamma+epsilonprime)-y2d1*(2*epsilonsource+epsilonprimeprime)-trace(Yb*sigmaAAs))
        sdpCondition(y,Gamma)+kron(Yb,eye(dimB))>=-transpose(gradf);
        y2d1*eye(dimAAs)>=Yb;
    cvx_end
    
    flag=strcat(cvx_status,'dual');

    Mepsilonprime = cvx_optval;
    debugging.Mepsilonprime = Mepsilonprime;
    if epsilon == 0
        zetaEp = 0;
    else 
       	zetaEp = 2 * epsilon * (dprime-1) * log(dprime/(epsilon*(dprime-1)));
    end
    lowerbound = Lepsilon + real(Mepsilonprime) - zetaEp; 
    debugging.epsilonprime1 = options.epsilonprime;
    debugging.zetaEp = zetaEp;
    %debugging.zetaEp2 =  2 * epsilon_per * (dprime-1) * log(dprime/(epsilon*(dprime-1)));
    debugging.lowerbound = lowerbound;
    % note: we use natural log in all the calculation
    % one needs to convert natural log to log2. 
end

function result = sdpCondition(dualY, GammaVector)
   	result = 0;
    for iConstraint = 1 : length(dualY)/2
        result = result + (dualY(iConstraint)-dualY(iConstraint+length(dualY)/2)) * GammaVector{iConstraint};
    end        
end