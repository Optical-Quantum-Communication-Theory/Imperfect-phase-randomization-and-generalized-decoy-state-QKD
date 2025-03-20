%% FUNCTION NAME: primalSolver
% Step 1 solver for dimension reduction method, modified for infinite 
% shield system specifically. Refer to Jan 25 notes for details.
%% Copyrights => will become open source
% Author: Twesh Upadhyaya
%
% Created: January 25, 2021
% 
%% Inputs
%Accepts a POVM Gamma, with corresponding expectations gamma, and W

function [fvalvec,rho,fval,gap,flag] = primalSolverAB(rho0, keyMap, sigmaAAs,...
    epsilonsource, Gamma, gammaL,gammaU, W, krausOperators)
    warning('off','CVX:UnsymmetricLMI');
    % options
    defaultOptions.maxiter = 30; % Maximum number of iterations
    defaultOptions.maxgap = 1e-10; % Maximum gap as specified by the Frank-Wolfe algorithm
    defaultOptions.linesearchprecision = 1e-20;
    defaultOptions.linesearchminstep = 1e-20;
    defaultOptions.linearconstrainttolerance = 10^-10;%1e-10;
    defaultOptions.initmethod = 1; % 1 for closest to rho0, 2 for maximize minimum eigenvalue
    defaultOptions.verbose = 0;
    
    options = defaultOptions; %options not passed in as argument
    
    dimAAs=size(sigmaAAs,1); %could also pass this in as a parameter
    dim = size(rho0,1);
    %dimB = dim/dimAAs;
    fval = NaN;
    gap = Inf;
    flag = 0; % success flag

    %commented out for now
    %[observables,independentCols] = removeLinearDependence(observables);
    %expectations = expectations(independentCols);

    % project rho0 onto the set of density matrices consistent with observations
    rho = closestDensityMatrix(rho0,sigmaAAs,dimAAs,epsilonsource,Gamma,gammaL,gammaU,W,options);
    rho = full(rho);
    prev_rho = rho;
    
%     if lambda_min(rho) < 0
%         flag = 1;
%         if options.verbose
%             display('Error: minimium eigenvalue less than 0.');
%             display('Try increasing the constraint tolerance.');
%         end
%         return;
%     end

    fvalvec=[inf inf inf inf];
    % Main optimization loop
    for i = 1:options.maxiter
        
        gradf = primalDf(rho,keyMap,krausOperators);
        
        [deltaRho,status] = subproblem(rho,sigmaAAs,dimAAs,epsilonsource,Gamma,gammaL,gammaU,W,gradf,options);
        
        if(strcmp(status, "Infeasible") || strcmp(status, 'Failed'))
            if(options.verbose)
                fprintf('Error on iteration %d; using rho from last iteration\n', i)
                if i==1
                    fprintf('Note: No optimization was successful. Initial guess for rho is being used.\n')
                end
            end
            rho = prev_rho;
            break;
        end
        
        % perform an exact line search
        optimoptions = optimset('TolX',options.linesearchprecision);
        stepSize = fminbnd(@(t)primalf(rho+t*deltaRho,keyMap,krausOperators),options.linesearchminstep,1,optimoptions);
        gap = abs(trace(rho*gradf.')-trace((rho+deltaRho)*gradf.'));
        
        %stopping condition if gap is small
        if gap < options.maxgap
            rho = rho + stepSize*deltaRho;
            break;
        end
               
        prev_rho=rho;
        rho = rho + stepSize*deltaRho;
        
        fval=primalf(rho,keyMap,krausOperators);
        fvalvec(i+4)=fval; %offset by 4 to allow indexing up to 4 previous values
        
        if ((fvalvec(i+4) > fvalvec(i+3)) && (fvalvec(i+3) > fvalvec(i+2)))... %stopping condition if fval is not decreasing last two times
           || ((fvalvec(i+4) > fvalvec(i+3)) && (fvalvec(i+3) < fvalvec(i+2)) &&... %stopping condition if oscillatory behaviour 
           (fvalvec(i+2) > fvalvec(i+1)) && (fvalvec(i+1) < fvalvec(i)))
            
            rho=prev_rho; %return the previous rho which has lower fval
            break;
        end

        disp('iter = '+string(i));
        disp('fval = '+string(fval));
        
        if i == options.maxiter
            flag = 1;
            display('Warning: Maximum iteration limit reached.');
        end
    end
    
    if options.verbose
        display(['Current gap: ',num2str(gap),'  Num iters: ',num2str(i)]);
    end
    
    if lambda_min(rho) < 0
        flag = 1;
        if options.verbose
            display('Warning: minimum eigenvalue less than 0.');
        end
    end
    
    fvalvec=fvalvec(5:end);
    fval = primalf(rho,keyMap,krausOperators);
end

function [deltaRho,status] = subproblem(rho,sigmaAAs,dimAAs,epsilonsource,Gamma,gammaL,gammaU,W,gradf,options)
    dim = size(rho,1);
    cvx_begin sdp quiet
        cvx_precision default
        cvx_solver mosek
        variable deltaRho(dim,dim) hermitian
        variable S(dimAAs,dimAAs) hermitian semidefinite
        minimize real(trace(transpose(gradf)*deltaRho))
        for j = 1:numel(Gamma)
            %expectations(j)-W-options.linearconstrainttolerance <= real(trace((rho+deltaRho)*observables{j})) <= expectations(j)+options.linearconstrainttolerance
            gammaL(j)-options.linearconstrainttolerance <= real(trace((rho+deltaRho)*Gamma{j})) <= gammaU(j)+options.linearconstrainttolerance;
        end
%         1-W <= real(trace(rho+deltaRho)) <= 1;
        real(trace(S)) <= 2*epsilonsource + options.linearconstrainttolerance;
        PartialTrace(rho+deltaRho,2,dimAAs) <= sigmaAAs+S;
        0 <= rho + deltaRho;
    cvx_end
    status = cvx_status;
end


function rho = closestDensityMatrix(rho0,sigmaAAs,dimAAs,epsilonsource,Gamma,gammaL,gammaU,W,options)
    dim = size(rho0,1);
    cvx_begin sdp quiet
        cvx_precision default
        cvx_solver mosek
        variable rho(dim,dim) hermitian semidefinite
        variable S(dimAAs,dimAAs) hermitian semidefinite
        if options.initmethod == 1
            minimize norm(rho0-rho)
        elseif options.initmethod == 2
            minimize -lambda_min(rho)
        end
        for j = 1:numel(Gamma)
            gammaL(j)-options.linearconstrainttolerance <= real(trace((rho)*Gamma{j})) <= gammaU(j)+options.linearconstrainttolerance;
        end
%         1-W <= real(trace(rho)) <= 1;
        real(trace(S)) <= 2*epsilonsource + options.linearconstrainttolerance;
        PartialTrace(rho,2,dimAAs) <= sigmaAAs+S;
    cvx_end
    cvx_status
end
