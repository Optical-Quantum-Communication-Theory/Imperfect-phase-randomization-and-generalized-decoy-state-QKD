%% Bob's Flag-State Squashed POVM for COW protocol in <=N subspace
% See for e.g.- Fig. 1 from 10.1088/1367-2630/ac1e41 for a schematic of the optical detection setup.
% Input:
% * N   : photon number corresponding to unsquashed <=N-photon subspace
% 
% * n   : number of bits in a block. For 3-State Protocol, this can be set to 1.
% 
% * t   : fraction of photons going into the monitoring line. For DPS, this
% can be set to 1.
% 
% Output:
% * POVM : a cell of the POVM elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Shlok Nahar          Last Updated: 27th May 2020
%
% This code is not written in the most readable manner. Feel free to email
% sanahar@uwaterloo.ca for any questions about it (or the functions called
% inside).


function POVM = COWPOVM(N, n, t)

    %% Notice we have 2n input modes, and we are interested in the space corresponding to <= N photons in ALL
    %% the modes together. This is important, because restricting to the <= N photon in each mode will give you 
    %% a bigger space than you need. However this is easier to construct.
    % Our first task is to construct a projector that converts the latter to the former.


    dim = (N+1) * nchoosek(2*n+N, 2*n-1)/(2*n); % sum (2*n+i-1) choose (2*n-1) from i = 0 to i = N
    
    projector = zeros(dim, (N+1)^(2*n)); % Projector onto the <= N photon subspace with block diagonal structure 
    % This projector takes a vector in the bigger space (of <= N photon in
    % EACH mode). If there are > N photon total, it produces zero. If there
    % are <=N photon total, it leaves it unchanged. We will construct this
    % project as follows.


    projector(1,1) = 1;   %the vacuum state in ALL input modes goes to the first dimension             

                                  
    for photonnumber = 1:1:N % Considering all photon numbers till the photon number cut-off

        %% Forming vector of all the ways that $x_1+....+x_{2n} (sum of photons in modes) = photonnumber$
        temp = GetPhotonCombo(photonnumber,2*n);
        
        l = size(temp,2); % Count number of such combinations that are there.
        for i = 1:1:l % Loop over all combinations
            combtemp = temp{i};
            index = 1; % Column that will be non-zero in our projection matrix
            for j = 1:1:2*n % Loop over all the elements in the combination
                index = index + combtemp(j)*(N+1)^(j-1); % |j2n....,j3,j2,j1> has only the 1+j1+(N+1)*j2+(N+1)^2*j3+... element non-zero
                                                        % Thus, this is the column of the projector that will be non-zero
            end
            usedDims = nchoosek(2*n+photonnumber-1,2*n); % Number of dims you have used for photon numbers < photonnumber
            % this can be computed using straightforward combinatorics.
            projector(usedDims+i, index) = 1;
        end
    end

    %% Now we create generic operators for the input modes to Bob's setup. 
    % We create a generic annihilation operator restricted to <=N-photon subspace
    a = zeros(N+1);
    
    for k = 1:1:N
        a(k,k+1) = sqrt(k);
    end
    
    % Annihilation operators u_i where the different spaces in the
    % tensor product represent different input time-bins 
    u = cell(2*n,1);
    for i = 1:1:2*n
        u{i} = a; % u_i=a
        for j=1:1:i-1                   % u_i='a' tensored with identity (i-1) times
            u{i}=kron(u{i},eye(N+1));
        end
        for j=1:1:(2*n-i)               % u_i=(2n-i) times identity tensor 'a' tensor identity tensored (i-1) times
            u{i}=kron(eye(N+1),u{i});
        end
    end

   
    % Project these operators into the <= N total photons subspace
    for i = 1:1:2*n
        u{i} = projector * u{i} * projector';
    end
    
    %% we have now constructed u{i}, which is the annihilation operator on the ith input time bin mode,
    %% restricted to the subspace of <=N TOTAL photons in all the time bins..

     %% For the next section, we specialize to Bob's COW detector setup.


    % We now define the vacuum vector |0,0,0,.. 0 2n times> in the <= N total photon subspace
    vac = projector * zket((N+1)^(2*n),1);

   
    m = cell(2*n+1,1); %Annihilation operators that represent the destructive interference of the adjacent time-bins    
    % These correspond to the "-" detector clicking in various time bins.

    p = cell(2*n+1,1); %Annihilation operators that represent the constructive interference of the adjacent time-bins
    % These correspond to the "+" detector clicking in various time bins.

    % Note that while writing these annihilation operators we have ignored
    % the terms that would come from the vacuum components of the beam
    % splitter as when writing the POVMs and imposing the condition that
    % said state was the vacuum state, these annihilation operators would
    % not play a role. See Varun Narasimhachar's masters thesis, equations
    % A.5 and A.6 for a more detailed explanation.
    
    % Annihilation operators for the first time slot
    m{1} = -0.5*sqrt(t)*u{1};
    p{1} = 0.5*sqrt(t)*u{1};
    
    % Annihilation operators for the last time slot
    m{2*n+1} = 0.5*sqrt(t)*u{2*n};
    p{2*n+1} = 0.5*sqrt(t)*u{2*n};
    
    % Annihilation operators for middle-slots
    for i = 2:1:2*n
        m{i} = 0.5*sqrt(t)*(u{i-1}-u{i});
        p{i} = 0.5*sqrt(t)*(u{i-1}+u{i});
    end

    % In the similar way, we construct the operators in the other basis,
    % and add it to a cell array which stores ALL the operators.

    numberofdetectors = 6*n+2; % Number of detectors in the protocol. Here we think of clicks in different time
                              % slots as different detectors.
    w = cell(numberofdetectors,1);
    for i = 1:1:2*n
        w{i} = m{i}; %first the destructive interference operators in the monitoring line
        w{2*n+1+i} = p{i}; %then the constructive interference operators in the monitoring line
        w{4*n+2+i} = sqrt(1-t)*u{i}; %then the operators for the detection line
    end
    w{2*n+1} = m{2*n+1};
    w{4*n+2} = p{2*n+1};
    
    
    %% We now construct the POVM elements corresponding to every possible click pattern
    %% We do NOT add flags, and we assume zero loss within the detectors. 

    % we construct a cell array POVM, where the second index is equal to
    % the total number of detectors clicking + 1. The first index tells you
    % the specific click pattern when that many detectors click.

    POVM{1,1} = vac*vac'; % no-click event (index2-1 = 0) 
    % no detectors click, and there is only one click pattern corresponding
    % to this.

    % In general, if we have k detectors clicks, then there are 6n+2 choose
    % k possible click patterns. The first index (for a given value of k , i.e index2 = k+1)
    % iterates over all these click patterns according to the MATLAB nchoosek
    % function.
    

    % To know which POVM corresponds to which detection event, look at how
    % the function nchoosek orders it.

   
    for numClicks = 1:1:numberofdetectors % This is the number of detectors that click in a particular detection event
        
        combo = nchoosek(w, numClicks); % Creates an array of all possible 
                                % combinations of i annihilation operators
                                % for the i click event
        l = size(combo, 1); % Find number of combinations of i annihilation operators that exist.
        % Taking all possible combinations of multi-clicks and adding them to
        % the set of POVMs.
        for j = 1:1:l % Loop over all combinations
            temp = combo(j,:); % Stores the annihilation operators corresponding to the jth click pattern.
            currentPOVM = zeros(dim); % Dummy variable to store the POVM as we loop over all photon numbers till the cut-off
            
          
            
            for photonnumber = 1:1:N % Considering all photon numbers till the photon number cut-off
                
                currentPOVM = currentPOVM + createPOVMblock({0}, temp, photonnumber, vac);  
                % Forms the POVMblock for a specific photonnumber,
                % for the combination stored in temp. We do this by
                % acting the operators (in temp) on the vacuum state in an
                % appropriate way.
                                                         
            end
            POVM{j,numClicks+1} = currentPOVM; % The second index is the number of detectors that click +1
              
        end
    end
    
%%%% Debugging

    x = POVM(~cellfun('isempty',POVM)); % To remove all the empty elements
    nPovmElm = numel(x); % Number of elements in the array
    
    % Checking if they sum to the identity
    sumElms = 0;
    for iElm = 1 : nPovmElm 
        sumElms = sumElms + x{iElm};
    end
    
    % Checking for positive-semidefiniteness
    for iElm = 1 : nPovmElm
        if eigs(-x{iElm},1) > 0 % Checking if the maximum eigenvalue of the
                                % negative of the POVM is positive
            error = eigs(-x{iElm},1);
            break
        end
    end
end
