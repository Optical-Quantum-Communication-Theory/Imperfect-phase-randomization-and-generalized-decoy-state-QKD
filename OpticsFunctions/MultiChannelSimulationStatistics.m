%% Bob's observed detection statistics for the 3-state protocol with a loss-only channel
%
% We assume that we have a loss-only channel and so whenever Alice sends
% coherent states, Bob receives coherent states with a reduced amplitude.
% Alice sends the state |mu, 0> with probability q/2, |0, mu> with
% probability q/2, |mu, mu> with probability 1-q. After channel loss, the
% states going into Bob's apparatus will be essentially the same as what
% Alice sends except with mu replaced with alpha, the amplitude of the
% coherent state after loss.
%
% Input:
% * p : Cell containing all the probabilities corresponding to the
%       different POVM elements as outputted by ChannelSimulationStatistics.
%
% * pzeroalpha : Cell containing the probabilities corresponding to the
%                different POVM elements conditioned on Alice having sent 
%                the zeroalpha state as outputted by ChannelSimulationStatistics.
%
% * palphazero : Cell containing the probabilities corresponding to the
%                different POVM elements conditioned on Alice having sent 
%                the zeroalpha state as outputted by ChannelSimulationStatistics.
%
% * palphaalpha : Cell containing the probabilities corresponding to the
%                different POVM elements conditioned on Alice having sent 
%                the decoy (alphaalpha) state as outputted by ChannelSimulationStatistics.
%
% Output:
% * P : Vector containing all the probabilities corresponding to the
%       different coarse grained POVM elements.
%
% * Pzeroalpha : Vector containing the probabilities corresponding to the
%                different coarse grained POVM elements conditioned on Alice having sent 
%                the zeroalpha state.
%
% * Palphazero : Vector containing the probabilities corresponding to the
%                different coarse grained POVM elements conditioned on Alice having sent 
%                the zeroalpha state.
%
% * Palphaalpha : Vector containing the probabilities corresponding to the
%                different coarse grained POVM elements conditioned on Alice having sent 
%                the decoy (alphaalpha) state.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Shlok Nahar                     Last Updated: 27th August 2020


function Q = MultiChannelSimulationStatistics(pzeroalpha, palphazero, palphaalpha, p)
    

    Q = zeros(4, 7);   % q(1,:) is Pzeroalpha, q(2,:) is Palphazero and so on...
    q = cell(4,1);
    q{1} = pzeroalpha;
    q{2} = palphazero;
    q{3} = palphaalpha;
    q{4} = p;
    
    for j = 1:1:4
        Q(j,1) = q{j}{1,1};         % No-click event
        probsum = q{j}{1,1};
        for i = 1:1:5
            Q(j,i+1) = q{j}{i,2};   % Single-click events
            probsum = probsum + Q(j,i+1);
        end
        Q(j,7) = 1 - probsum;           % Multi-click event
    end
    Pzeroalpha = Q(1,:); 
    Palphazero = Q(2,:);
    Palphaalpha = Q(3,:);
    P = Q(4,:);
end