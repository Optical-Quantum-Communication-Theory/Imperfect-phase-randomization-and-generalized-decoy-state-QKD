%% FUNCTION NAME: ThreeState_PrivAmp_decoy
% This function allows you to calculate the privacy amplification (PA) term
% in the key rate formula of the 3-state protocol with decoy states.
%
% Input:
% * N : preserving up to N-photon subspace on Bob's side
%
% * nA : cut-off for eigenvectors on Alice's shield system
%
% * decoySDPcutoff: cut-off for photon number projection on Alice's signal
%                   states in the decoy-state SDP
%
% * a : probability of Alice's sending a completely phase randomised state
%
% * amp : Amplitude of the coherent state produced by the laser
%
% * eta : Channel loss
%
% * prob : probability of producing a 0/1 state (as opposed to the decoy
%          state)
%
% * t : fraction of photons going into the monitoring line
%
% * iter: maximum number of iterations in primalSolver
%
% * rho0: initial guess for optimisation in primalSolver
%  
% Output:
% * Rho: a cell of sub-optimal solutions (density matrices) to the primal SDP problem conditioned on Alice sending out n-photon signal
% 
% * primalfvals: an array of upperbounds on the optimal PA term from the primal SDP conditioned on Alice sending out n-photon signal
%
% * dualfvals: an array of lowerbounds on the optimal PA term from the dual SDP conditioned on Alice sending out n-photon signal
% 
% * gaps: [gaps(i) = abs(trace(dRho{i} * gradf.'))] where i labels which n-photon signal state Alice sent
% 
% * flagprimal: flags in primal optimisation
% 
% * flagdual: flags in dual optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Shlok Nahar                Last updated: 16th July 2020
function [Rho, primalfvals, dualfvals, keyrate, gaps, flagprimal, flagdual] = ThreeState_PrivAmp_decoy(N, nA, decoySDPcutoff, a, amp, eta, prob, t, iter, rho0)
    
    % The completely phase random part is N+1 dimensional. The flag system
    % and the coherent part add one dimension each since those are already
    % pure.
    dimAshield = 1; % Dimension of Alice's shield system
    
    dimA = 3*dimAshield; % dimension of Alice's system in the PM scheme(3 possible states she can send)
%     dimAshield = 1;
%     dimA = 2;
    nondiagdimB = (N+1)*(N+2)/2;  % <=N-photon dimension
    diagdimB = 0;   % flag-state dimension (number of coarse-grained 'click' events)
    dimB = nondiagdimB + diagdimB; % Bob's total dimension
    dim = dimA * dimB;
    
    vecZero = zket(2,1);
    vecOne = zket(2,2);
    
    if nargin == 8 % if no input rho0
        rho0 = eye(dim)/dim;
        options.initmethod = 2;
        
    elseif nargin == 9 % specified initial guess for rho0
        options.initmethod = 1;
    end
    
    
    BPOVMs = COWPOVM(N, 1, t);
    BPOVM = CoarseGrainedPOVMs(BPOVMs);
    [BobPOVM, BobPOVMNoflag] = MultiCoarseGrainedPOVMs(BPOVM);
    KraussBob = MultiBobKrauss(BobPOVMNoflag);
    KraussAlice = AliceKrauss(dimAshield,1);
  
    krausOp = {kron(KraussAlice, KraussBob)};
    keyProj1 = kron(vecZero*vecZero', eye(dimAshield*dimB)); 
    keyProj2 = kron(vecOne*vecOne', eye(dimAshield*dimB));
    keyMap = {keyProj1, keyProj2};
    
    %% Primal CVX Settings
    options.verbose = 1;
    options.linearconstrainttolerance = 1e-10;

    options.maxgap = 1e-7;
    options.maxiter = iter;

    AlicePOVM = {kron(diag(zket(3,1)),eye(dimAshield)), ...
                kron(diag(zket(3,2)),eye(dimAshield)), ...
                kron(diag(zket(3,3)),eye(dimAshield)), ...
                eye(dimA)};
    pA = [prob/2; prob/2; 1-prob; 1];
    q1 = [];
    
%     n = 17;
    n = decoySDPcutoff;
    for i = 0:1:n
        if(poisspdf(i,amp(1)^2) <= eps)
            n = i-1;
            break;
        end
    end
    
    for k = 1:1:length(amp)
    [p{4*k-3}, p{4*k-2}, p{4*k-1}, p{4*k}] = ChannelSimulationStatistics(sqrt(eta)*amp(k), prob, t);
    q{k} = MultiChannelSimulationStatistics(p{4*k-3}, p{4*k-2}, p{4*k-1}, p{4*k});
    q1 = [q1;MultiSimulationStatistics(sqrt(eta)*amp(k), t)];
    
    alpha = amp(k);
    % Projected State
    
    
    rho{k} = zeros(n+1);
    for i = 0:1:n
        for j = 0:1:n
            if poisspdf(i,alpha^2) > eps
                if i == j
                    rho{k}(i+1,j+1) = poisspdf(i,alpha^2);
                else
                    rho{k}(i+1,j+1) = (1-a)*sqrt(poisspdf(i,alpha^2)*poisspdf(j,alpha^2));
                end
            end
        end
    end
    end
    
    Fproj = trace(rho{1}); % Fidelity between approx. rho and infinite rho
    if Fproj>=1
        Fproj = 1-eps;
    end
    wt = 1-Fproj;
    epsilonproj = (1-a)*sqrt(wt);
    
    [psi,lambda,epsilonvec] = closestEigenvectors(rho{1}, epsilonproj, nA+1); % Eigenvector of the projected rho with the maximum distance between this eigenvector and the infinite dimensional states eigenvector

    V{1} = kron(zket(n+1,1),eye(n+1));
    V{2} = kron(eye(n+1),zket(n+1,1));
    V{3} = 0;
    for i = 0:1:n
        for j = 0:1:i
            V{3} = V{3} + sqrt(nchoosek(i,j)/2^i)*kron(zket(n+1,j+1),zket(n+1,i-j+1))*zket(n+1,i+1)';
        end
    end
    knownStates = {};
    for k = 1:1:length(amp)
        knownStates{3*k-2} = V{1}*rho{k}*V{1}';
        knownStates{3*k-1} = V{2}*rho{k}*V{2}';
        knownStates{3*k} = V{3}*rho{k}*V{3}';
    end 
    
    fullexpect = [];
    for i = 1:1:3
        expect = [];
        for k = 1:1:length(amp)
            PgeqNx = (1-ThreeStateleqN (N, t, p{4*k-4+i})); % Upper bound on greater than equal to N subspace
%             Wb(3*k-3+i) = PgeqNx; % single intensity
            Wb(k) = PgeqNx; % multiple intensities
            consStat = q1(i+3*k-3,:);
            consPOVM = BobPOVMNoflag;
            expect = [expect, consStat'];
        end
        [YL{i}, YU{i}] = decoyBoundsSDP(consPOVM,consPOVM, rho, {psi*psi'}, expect, Wb, a);
%         fullexpect(:,i) = expect; % For single intensity
    end
        
%     for i = 1:1:3 For single intensity
%         [YL{i}, YU{i}] = decoyBoundsSDP(consPOVM,consPOVM, knownStates, {V{i}*psi*psi'*V{i}'}, fullexpect, Wb, a);
%     end

    %   All equality constraints become inequality constraints in the projection method
    Gamma = {};
    gammaL = [];
    gammaU = [];
    
    
    for i = 1:1:3
        for j = 1:length(consPOVM)
            Gamma = [Gamma; kron(AlicePOVM{i}, consPOVM{j})];
            gammaL = [gammaL; max(pA(i)*YL{i}(j)-epsilonvec,0)];
            gammaU = [gammaU; min(pA(i)*YU{i}(j)+epsilonvec,pA(i))];
        end
    end
    
    W=0;
    rhoA = zeros(3);
    for i = 1:1:3
        for j = 1:1:3
            rhoA(i,j) = sqrt(pA(i)*pA(j))*psi'*V{i}'*V{j}*psi;
        end
    end
    
    %% Primal CVX
    [fvalvec,rho,primalfvals,gaps,flagprimal] = primalSolverAB(rho0, keyMap, rhoA,...
                        epsilonvec, Gamma, gammaL, gammaU, W, krausOp);

    % Perturb the suboptimal rho if its eigenvalues are not all
    % positive to avoid the 'zetaEp' in 'step2Solver.m' being to big
    eigMin = lambda_min(rho);
    if eigMin < 0
        epsilon = -eigMin*dim;
        if epsilon < 1
            rho = (1-epsilon)*rho + epsilon *eye(dim)/dim;
            disp('lambda_min(Rho) with lambda_min before perturbation = '+string(eigMin));
            disp('lambda_min(Rho) with lambda_min = '+string(lambda_min(rho)));
        else
           disp('Warning: epsilon >= 1');
        end
    end
    Rho = rho;

    %% Dual CVX
    disp('running dual cvx...');
    tic
    
    cons = zeros(1, numel(Gamma));  % +numel(IneqOp)
    for iBasisElm = 1:numel(Gamma)
            if real(trace(rho * Gamma{iBasisElm})) > gammaU(iBasisElm)
                cons(iBasisElm) = abs(real(trace(rho * Gamma{iBasisElm}) - gammaU(iBasisElm)));
            end
    end
    option2 = {};
    option2.epsilonprime = max(cons); % maximum constraint violation
    disp('max cons = '+string(max(cons)));
    option2.epsilonprimeprime = 1e-10
    
    %% Dubugging
    
%     testA = PartialTrace(Rho, 2, dimA, dimB);
%     testB = PartialTrace(Rho, 1, dimA, dimB);
%     for i = 1:dimA
%         for j =1 :dimA
%             if(abs(testA(i,j))<1e-8)
%                 testA(i,j) = 0;
%             else
%                 testA(i,j) = testA(i,j);
%             end
%         end
%     end
    
%     for i = 1:dimB
%         for j =1 :dimB
%             if(abs(testB(i,j))<1e-8)
%                 testB(i,j) = 0;
%             else
%                 testB(i,j) = testB(i,j);
%             end
%         end
%     end
    
    [dualfvals, debugging, flagdual] = step2SolverAB(Rho, keyMap, rhoA,...
                                epsilonvec, Gamma, gammaL, gammaU, sum(W), krausOp, option2);
%     disp('dualfval = '+string(dualfvals(nA))+' at nA = '+string(nA+nA_min-1));
    if dualfvals <= 0 || (lambda-epsilonproj) <= 0
        keyrate = 0;
    else
        keyrate = dualfvals*(lambda-epsilonproj);
    end
end