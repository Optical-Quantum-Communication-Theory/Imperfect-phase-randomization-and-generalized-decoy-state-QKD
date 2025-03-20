%% Alice's Krauss operators  for the 3-State protocol
% The announcements correspond to the basis in which Alice's prepared state
% is.
% Input:
% * dimAshield : dimension of the shield system
% 
% * nd : Number of decoy settings
%
% Output:
% * AKrauss : a cell of the Krauss operators
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Shlok Nahar                        Last Updated: 29th May 2020


function AKrauss = AliceKrauss(dimAshield, nd)
    
    % Dimension of Alice's shield system
    
    dim = dimAshield;

    vecZero = zket(2,1);
    vecOne = zket(2,2);

    % Projecting onto Alice's system and adding the key map register corresponding to her measurement outcome with Identity 	on the shield system
    PZero = kron(vecZero*zket(3*nd,1)',eye(dim)); % Measurement result 0.
    POne = kron(vecOne*zket(3*nd,2)',eye(dim)); % Measurement result 1.
%    PTwo = kron(zket(3,3)*zket(3,3)',eye(dim)); % Measurement result 2.
    
    
    % Krauss operator for the 0/1 basis
    KraussBit = PZero + POne;

%    KraussBit = kron(vecZero,sqrtm(PZero)) + kron(vecOne,sqrtm(POne));

    % Krauss operator for the 0/1 basis
%    KraussBit = kron(kron(sqrtm(PZero),vecZero),vecZero) + kron(kron(sqrtm(POne),vecZero),vecOne);
    % Krauss operator for the +/- basis
 %   KraussDecoy = kron(kron(sqrtm(PTwo),vecOne),vecZero); Don't really
 %   need this.
    
    AKrauss = KraussBit;
%    AKrauss{2} = KraussDecoy;
end