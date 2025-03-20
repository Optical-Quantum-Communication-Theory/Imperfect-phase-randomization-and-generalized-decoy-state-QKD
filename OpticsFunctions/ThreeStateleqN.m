%% Calculate the minimum weight of the <=N-photon subspace of the signal state of Three State Protocol
% Input:
% * N  : number of photons 
% * t  : fraction of photons going into the monitoring line
% * pbx: Bob's outcome probability cell conditioned on Alice choice of
%        state where x = 2 indicates her sending the "+" state
%
% Output:
% * PleqN: probability of signal state in <=N-photon subspace
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Shlok Nahar                                     21st July 2020
function PleqN = ThreeStateleqN(N, t, pbx)
    
    % Refer to CoarseGrainedPOVMS.m for the numbers for pbx

    data = 0; % prob of click in data line or in middle time bin of minus detector
    minus = 0; % minus detector-only click prob
    
    data = pbx{4,2} + pbx{5,2} + pbx{2,2} + ...  % single clicks
           pbx{6,3} + pbx{7,3} + pbx{10,3} + ... % double clicks
           pbx{9,4};                             % triple clicks
    
    minus = pbx{1,2} + pbx{2,2} + pbx{3,2} + ... % single clicks
            pbx{1,3} + pbx{2,3} + pbx{5,3} + ... % double clicks
            pbx{1,4};                            % triple clicks
    
    % Cross-click prob
    cc = 1 - data - minus - pbx{1,1} + ...
         pbx{2,2};                         % to avoid double-counting
    if cc < 0    % cc maybe slightly negative e.g. -1e-16 due to numerical precision
        cc = 0;
    end
    
    % Refer to Eq. (C26) from PhysRevApplied.20.064031 for this
    cc_min_0 = 0;
    ccMinNPlusOne = 1 - t^(N+1)-(1-t/4)^(N+1) + (3*t/4)^(N+1);

    PleqN = 1 - (cc - cc_min_0)/(ccMinNPlusOne - cc_min_0);
    
    if PleqN > 1
        PleqN = 1;
    elseif PleqN < 0
        PleqN = 0;
    end

end
