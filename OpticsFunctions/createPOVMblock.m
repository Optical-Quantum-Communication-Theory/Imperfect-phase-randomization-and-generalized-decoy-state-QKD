 %% Forming the POVM from a specified set of annihilation operators assuming Nphoton photons
% 

% We are given as input the following : 
% Input:
% * indices   : Vector of numbers that signify the indices for the previous
%         summations. This stores the number of photons you have already
%         put in the first few detectors. We start with {0}. 
%         

% * annOp   : Cell containing all the annihilation operators for the different
%         detectors and time slots.

% * NPhoton   : Total number of photons that are detected.
%
% * vac : The vacuum state |0>.


% We have to sum over all  photon numbers for the different modes with the 
% constraint that the sum of all the photon numbers in
% the modes is NPhoton. Thus, we have multiple summations that are dependent
% on each other through the constraint. This function is recursive.
%
% We will first populate the first modes with some number of photons.
% Now we have to populate the remaining modes with a smaller number of
% photons. Thus, this function can be called recursively, until we only
% have a single modes left. In this case, ALL the remaining photons must go
% into that mode.
% 
%
% This code can be used to find multi-click POVM elements for any protocol
% if you know the annihilation operators for each detector.



% Output:
% * P : a (2n,1) cell which will be the POVMs after adding the contribution
%       due to the detector and time-slot represented by the annihilation
%       operator m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Shlok Nahar          Last Updated: 13th Jan 2020

function POVMblock = createPOVMblock(indices, annOp, NPhoton, vac)

    if NPhoton >= length(annOp) % Need to have more photons than click events

        remainingPhotons = NPhoton-sum([indices{:}])-(length(annOp)-1); % The total number of photons left for the rest of the summations
        % we substract the number of photons we have already put in
        % detectors.
        % we then subtract the number of remaining detectors (because there
        % should be atleast 1 photon in each of the REMAINING detectors).

        
        if(length(annOp) > 1)
            POVMblock = 0; % initialising sum as 0


            for i = 1:1:remainingPhotons % sum over all photon numbers from 1 to the number of photons left
                newannOp = annOp;
                newannOp(1) = []; % Deleting the annihilation operator just used in this sum
                
                newindices = indices;
                newindices{end + 1} = i; % Adding the number of photons detected in this operator to the vector
                POVMblock = POVMblock + annOp{1}'^i * createPOVMblock(newindices, newannOp, NPhoton, vac) * annOp{1}^i/factorial(i);
            end
        else % Base case when we have just one annihilation operator
            POVMblock = annOp{1}'^remainingPhotons * (vac) * (vac') * annOp{1}^remainingPhotons/factorial(remainingPhotons);
        end
    else
        POVMblock = zeros(size(vac,2)); % Zero matrix of whatever dimension space we have.
    end
end