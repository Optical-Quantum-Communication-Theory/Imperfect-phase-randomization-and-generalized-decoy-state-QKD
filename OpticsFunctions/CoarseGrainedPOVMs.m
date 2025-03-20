%% Bob's coarse-grained POVM operators  for the 3-State protocol
% We do not actually have the "plus" detector in the experimental
% implementation of the protocol. So any detection events with clicks in
% the "plus" detector and clicks in any other detector would be coarse
% grained to clicks in just the other detector and all clicks in just the 
% "plus" detector would be coarse grained to the no-click event.
%
% Input:
% * POVM   : POVMs as outputted by the COWPOVM function with n=1
% 
% Output:
% * CPOVM : Coarse grained POVMs
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Shlok Nahar                     Last Updated: 27th August 2020

%% Note that there is clearly a better way to do this, by writing functions that 
%% switch between click patterns, and the (index1,index2) style of storing POVMs. 
%% Unfortunately, this code does not follow that logic, and is more handwritten.  

function [CPOVM] = CoarseGrainedPOVMs(POVM)


    CPOVM = cell(10, 6); % The second index i is for the number of detectors
                         % that click and the first will tell us exactly
                         % what (i-1)-click event has occured.
    
    %% The coarse grained no-click event
    
    CPOVM{1,1} = POVM{1,1} + ...                             % Vacuum event
                 POVM{4,2} + POVM{5,2} + POVM{6,2} + ...     % Single-click events
                 POVM{19,3} + POVM{20,3} + POVM{23,3} + ...  % Double-click events
                 POVM{47,4};                                 % Triple-click event
                 
             
    
    %% Coarse grained single-click events
    
    % The "-" detector in the first time slot (1)
    CPOVM{1,2} = POVM{1,2} + ...                             % Single-click event
                 POVM{3,3} + POVM{4,3} + POVM{5,3} + ...     % Double-click events
                 POVM{12,4} + POVM{13,4} + POVM{16,4} + ...  % Triple-click events
                 POVM{26,5};                                 % Quadruple-click event

    % The "-" detector in the middle time slot (2)
    CPOVM{2,2} = POVM{2,2} + ...                             % Single-click event
                 POVM{9,3} + POVM{10,3} + POVM{11,3} + ...   % Double-click events
                 POVM{27,4} + POVM{28,4} + POVM{31,4} + ...  % Triple-click events
                 POVM{46,5};                                 % Quadruple-click event

    % The "-" detector in the last time slot (3)
    CPOVM{3,2} = POVM{3,2} + ...                             % Single-click event
                 POVM{14,3} + POVM{15,3} + POVM{16,3} + ...  % Double-click events
                 POVM{37,4} + POVM{38,4} + POVM{41,4} + ...  % Triple-click events
                 POVM{56,5};                                 % Quadruple-click event

    % The "0/1" detector in the first time slot (7)
    CPOVM{4,2} = POVM{7,2} + ...                             % Single-click event
                 POVM{21,3} + POVM{24,3} + POVM{26,3} + ...  % Double-click events
                 POVM{48,4} + POVM{50,4} + POVM{53,4} + ...  % Triple-click events
                 POVM{66,5};                                 % Quadruple-click event
    
    % The "0/1" detector in the last time slot (8)
    CPOVM{5,2} = POVM{8,2} + ...                             % Single-click event
                 POVM{22,3} + POVM{25,3} + POVM{27,3} + ...  % Double-click events
                 POVM{49,4} + POVM{51,4} + POVM{54,4} + ...  % Triple-click events
                 POVM{67,5};                                 % Quadruple-click event
             
    %% Coarse-grained double-click events
    
    % Double-click in 1 and 2
    CPOVM{1,3} = POVM{1,3} + ...                             % Double-click event
                 POVM{2,4} + POVM{3,4} + POVM{4,4} + ...     % Triple-click events
                 POVM{6,5} + POVM{7,5} + POVM{10,5} + ...    % Quadruple-click events
                 POVM{11,6};                                 % Quintuple-click event
    
    % Double-click in 1 and 3
    CPOVM{2,3} = POVM{2,3} + ...                             % Double-click event
                 POVM{7,4} + POVM{8,4} + POVM{9,4} + ...     % Triple-click events
                 POVM{16,5} + POVM{17,5} + POVM{20,5} + ...  % Quadruple-click events
                 POVM{21,6};                                 % Quintuple-click event
             
    % Double-click in 1 and 7
    CPOVM{3,3} = POVM{6,3} + ...                             % Double-click event
                 POVM{14,4} + POVM{17,4} + POVM{19,4} + ...  % Triple-click events
                 POVM{27,5} + POVM{29,5} + POVM{32,5} + ...  % Quadruple-click events
                 POVM{31,6};                                 % Quintuple-click event
    
    % Double-click in 1 and 8
    CPOVM{4,3} = POVM{7,3} + ...                             % Double-click event    
                 POVM{15,4} + POVM{18,4} + POVM{20,4} + ...  % Triple-click events
                 POVM{28,5} + POVM{30,5} + POVM{33,5} + ...  % Quadruple-click events
                 POVM{32,6};                                 % Quintuple-click event
    
    % Double-click in 2 and 3
    CPOVM{5,3} = POVM{8,3} + ...                             % Double-click event
                 POVM{22,4} + POVM{23,4} + POVM{24,4} + ...  % Triple-click events
                 POVM{36,5} + POVM{37,5} + POVM{40,5} + ...  % Quadruple-click events
                 POVM{36,6};                                 % Quintuple-click event
    
    % Double-click in 2 and 7
    CPOVM{6,3} = POVM{12,3} + ...                            % Double-click event
                 POVM{29,4} + POVM{32,4} + POVM{34,4} + ...  % Triple-click events
                 POVM{47,5} + POVM{49,5} + POVM{52,5} + ...  % Quadruple-click events
                 POVM{46,6};                                 % Quintuple-click event
    
    % Double-click in 2 and 8
    CPOVM{7,3} = POVM{13,3} + ...                            % Double-click event
                 POVM{30,4} + POVM{33,4} + POVM{35,4} + ...  % Triple-click events
                 POVM{48,5} + POVM{50,5} + POVM{53,5} + ...  % Quadruple-click events
                 POVM{47,6};                                 % Quintuple-click event
    
    % Double-click in 3 and 7
    CPOVM{8,3} = POVM{17,3} + ...                            % Double-click event
                 POVM{39,4} + POVM{42,4} + POVM{44,4} + ...  % Triple-click events
                 POVM{57,5} + POVM{59,5} + POVM{62,5} + ...  % Quadruple-click events
                 POVM{51,6};                                 % Quintuple-click event
    
    % Double-click in 3 and 8
    CPOVM{9,3} = POVM{18,3} + ...                            % Double-click event
                 POVM{40,4} + POVM{43,4} + POVM{45,4} + ...  % Triple-click events
                 POVM{58,5} + POVM{60,5} + POVM{63,5} + ...  % Quadruple-click events
                 POVM{52,6};                                 % Quintuple-click event
    
    % Double-click in 7 and 8
    CPOVM{10,3} = POVM{28,3} + ...                           % Double-click event
                  POVM{52,4} + POVM{55,4} + POVM{56,4} + ... % Triple-click events
                  POVM{68,5} + POVM{69,5} + POVM{70,5} + ... % Quadruple-click events
                  POVM{56,6};                                % Quintuple-click event

    %% Coarse-grained triple-click events
    
    % Triple-click in 1,2 and 3
    CPOVM{1,4} = POVM{1,4} + ...                             % Triple-click event
                 POVM{1,5} + POVM{2,5} + POVM{3,5} + ...     % Quadruple-click events
                 POVM{1,6} + POVM{2,6} + POVM{5,6} + ...     % Quintuple-click events
                 POVM{1,7};                                  % Sextuple-click event
    
    % Triple-click in 1,2 and 7
    CPOVM{2,4} = POVM{5,4} + ...                             % Triple-click event
                 POVM{8,5} + POVM{11,5} + POVM{13,5} + ...   % Quadruple-click events
                 POVM{12,6} + POVM{14,6} + POVM{17,6} + ...  % Quintuple-click events
                 POVM{11,7};                                 % Sextuple-click event
    
    % Triple-click in 1,2 and 8
    CPOVM{3,4} = POVM{6,4} + ...                             % Triple-click event
                 POVM{9,5} + POVM{12,5} + POVM{14,5} + ...   % Quadruple-click events
                 POVM{13,6} + POVM{15,6} + POVM{18,6} + ...  % Quintuple-click events
                 POVM{12,7};                                 % Sextuple-click event
    
    % Triple-click in 1,3 and 7
    CPOVM{4,4} = POVM{10,4} + ...                            % Triple-click event
                 POVM{18,5} + POVM{21,5} + POVM{23,5} + ...  % Quadruple-click events
                 POVM{22,6} + POVM{24,6} + POVM{27,6} + ...  % Quintuple-click events
                 POVM{16,7};                                 % Sextuple-click event
    
    % Triple-click in 1,3 and 8
    CPOVM{5,4} = POVM{11,4} + ...                            % Triple-click event
                 POVM{19,5} + POVM{22,5} + POVM{24,5} + ...  % Quadruple-click events
                 POVM{23,6} + POVM{25,6} + POVM{28,6} + ...  % Quintuple-click events
                 POVM{17,7};                                 % Sextuple-click event
    
    % Triple-click in 1,7 and 8
    CPOVM{6,4} = POVM{21,4} + ...                            % Triple-click event
                 POVM{31,5} + POVM{34,5} + POVM{35,5} + ...  % Quadruple-click events
                 POVM{33,6} + POVM{34,6} + POVM{35,6} + ...  % Quintuple-click events
                 POVM{21,7};                                 % Sextuple-click event
    
    % Triple-click in 2,3 and 7
    CPOVM{7,4} = POVM{25,4} + ...                            % Triple-click event
                 POVM{38,5} + POVM{41,5} + POVM{43,5} + ...  % Quadruple-click events
                 POVM{37,6} + POVM{39,6} + POVM{41,6} + ...  % Quintuple-click events
                 POVM{22,7};                                 % Sextuple-click event
    
    % Triple-click in 2,3 and 8
    CPOVM{8,4} = POVM{26,4} + ...                            % Triple-click event
                 POVM{39,5} + POVM{42,5} + POVM{44,5} + ...  % Quadruple-click events
                 POVM{38,6} + POVM{40,6} + POVM{42,6} + ...  % Quintuple-click events
                 POVM{23,7};                                 % Sextuple-click event
    
    % Triple-click in 2,7 and 8
    CPOVM{9,4} = POVM{36,4} + ...                            % Triple-click event
                 POVM{51,5} + POVM{54,5} + POVM{55,5} + ...  % Quadruple-click events
                 POVM{48,6} + POVM{49,6} + POVM{50,6} + ...  % Quintuple-click events
                 POVM{27,7};                                 % Sextuple-click event
    
    % Triple-click in 3,7 and 8
    CPOVM{10,4} = POVM{46,4} + ...                           % Triple-click event
                  POVM{61,5} + POVM{64,5} + POVM{65,5} + ... % Quadruple-click events
                  POVM{53,6} + POVM{54,6} + POVM{56,6} + ...  % Quintuple-click events
                  POVM{28,7};                                 % Sextuple-click event
    
    %% Coarse-grained quadruple-click events

    % Quadraple-click in 1,2,3 and 7
    CPOVM{1,5} = POVM{4,5} + ...                             % Quadruple-click event
                 POVM{3,6} + POVM{6,6} + POVM{8,6} + ...     % Quintuple-click events
                 POVM{2,7} + POVM{4,7} + POVM{7,6} + ...     % Sextuple-click events
                 POVM{1,8};                                  % Septuple-click event
    
    % Quadraple-click in 1,2,3 and 8
    CPOVM{2,5} = POVM{5,5} + ...                             % Quadruple-click event
                 POVM{4,6} + POVM{7,6} + POVM{9,6} + ...     % Quintuple-click events
                 POVM{3,7} + POVM{5,7} + POVM{8,6} + ...     % Sextuple-click events
                 POVM{2,8};                                  % Septuple-click event
                 
    % Quadraple-click in 1,2,7 and 8
    CPOVM{3,5} = POVM{15,5} + ...                            % Quadruple-click event
                 POVM{16,6} + POVM{19,6} + POVM{20,6} + ...  % Quintuple-click events
                 POVM{13,7} + POVM{14,7} + POVM{15,6} + ...  % Sextuple-click events
                 POVM{6,8};                                  % Septuple-click event
    
    % Quadraple-click in 1,3,7 and 8
    CPOVM{4,5} = POVM{25,5} + ...                            % Quadruple-click event
                 POVM{26,6} + POVM{29,6} + POVM{30,6} + ...  % Quintuple-click events
                 POVM{18,7} + POVM{19,7} + POVM{20,6} + ...  % Sextuple-click events
                 POVM{7,8};                                  % Septuple-click event
    
    % Quadraple-click in 2,3,7 and 8
    CPOVM{5,5} = POVM{45,5} + ...                            % Quadruple-click event
                 POVM{41,6} + POVM{44,6} + POVM{45,6} + ...  % Quintuple-click events
                 POVM{24,7} + POVM{25,7} + POVM{26,6} + ...  % Sextuple-click events
                 POVM{8,8};                                  % Septuple-click event
    
    %% Coarse-grained all click event

    CPOVM{1,6} = POVM{10,6} + ...                            % Quintuple-click event
                 POVM{6,7} + POVM{9,7} + POVM{10,6} + ...    % Sextuple-click events
                 POVM{3,8} + POVM{4,8} + POVM{5,8} + ...     % Septuple-click events
                 POVM{1,9};                                  % All-click event
    
%% Debugging

    x = CPOVM(~cellfun('isempty',CPOVM)); % To remove all the empty elements
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





