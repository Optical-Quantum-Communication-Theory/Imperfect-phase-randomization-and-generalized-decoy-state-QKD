%% Bob's Krauss operators  for the 3-State protocol
% The announcements correspond to the basis in which Bob got an outcome, or
% if the detection event should be discarded. The chosen post-processing here 
% is as follows:
% 1) The 0/1 basis corresponds to all clicks in the data line and outer clicks
% in the "minus" detector. 
% 2) The +/- basis corresponds to all clicks in the middle time slot of the 
% "minus" detector.
% 3) All multi-clicks are discarded.
% Input:
% * POVM   : POVMs as outputted by the MultiCoarseGrainedPOVM function
% 
% Output:
% * BKrauss : a cell of the Krauss operators
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Shlok Nahar                       Last Updated: 11th June 2020


function BKrauss = MultiBobKrauss(CPOVM)
    
    dim = size(CPOVM{1,1},1); % Dimension of the flag-state POVMs
    
    % Coarse grained POVM operators
    PZero = zeros(dim); % Measurement result 0. 
    POne = zeros(dim); % Measurement result 1.
      
    PZero = CPOVM{1,5} + ...   % Single-click in just the 0/1 detector. (7)
            CPOVM{1,2};        % Single-click in the "-" detector in the first time slot (1)
    POne = CPOVM{1,6} + ...    % Single-click in just the 0/1 detector. (8)
           CPOVM{1,4};         % Single-click in the "-" detector in the last time slot (3)
    
    % Using Nicky's sqrtFlagStatePOVM function to compute the sqrt of the
    % flag-state POVMs. He said the in-built function gave some errors
    % sometimes so using this just to be safe.
%     sqrtPZero = sqrtFlagStatePOVM(PZero, dim, dim-7);
%     sqrtPOne = sqrtFlagStatePOVM(POne, dim, dim-7);
    
    vecZero = zket(2,1);
    vecOne = zket(2,2);

    % Krauss operator for the 0/1 basis
%     KraussData = sqrtFlagStatePOVM(PZero + POne, dim, dim-7);
    KraussData = sqrt(PZero+POne);
    % Krauss operator for the 0/1 basis
%    KraussData = kron(kron(sqrtPZero,vecZero),vecZero) + kron(kron(sqrtPOne,vecZero),vecOne);
    % Krauss operator for the +/- basis
%    KraussMonitoring = kron(kron(sqrtPMinus,vecOne),vecOne); Don't really
%    need this.
    
    BKrauss = KraussData;
%    BKrauss{2} = KraussMonitoring;
end