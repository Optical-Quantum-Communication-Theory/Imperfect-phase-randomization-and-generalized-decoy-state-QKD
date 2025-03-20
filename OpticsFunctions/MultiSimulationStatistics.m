%% Bob's observed detection statistics for the 3-state protocol with a loss-only channel
%
% We assume that we have a loss-only channel and so whenever Alice sends
% coherent states, Bob receives coherent states with a reduced amplitude.
% Alice sends the state |mu, 0> with probability q/2, |0, mu> with
% probability q/2, |mu/sqrt(2), mu/sqrt(2)> with probability (1-q). After channel loss, the
% states going into Bob's apparatus will be essentially the same as what
% Alice sends except with mu replaced with alpha, the amplitude of the
% coherent state after loss.
%
% Input:
% * alpha   : Amplitude of the coherent states entering Bob's apparatus.
% 
% * t       : Fraction of photons going into the monitoring line.
%
% Output:
% * p : Cell containing all the probabilities corresponding to the
%       different POVM elements.
%
% * pzeroalpha : Cell containing the probabilities corresponding to the
%                different POVM elements conditioned on Alice having sent 
%                the zeroalpha state.
%
% * palphazero : Cell containing the probabilities corresponding to the
%                different POVM elements conditioned on Alice having sent 
%                the zeroalpha state.
%
% * palphaalpha : Cell containing the probabilities corresponding to the
%                different POVM elements conditioned on Alice having sent 
%                the decoy (alphaalpha) state.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Shlok Nahar                       Last Updated: 16th June 2020


function p = MultiSimulationStatistics(alpha, t)
    
    p = zeros(3,7);
    a = exp(-t*alpha^2/4);
    b = exp(-(1-t)*alpha^2);
    pzeroalpha(1) = a^2*b;
    pzeroalpha(2) = (1-a)*b*a;
    pzeroalpha(3) = (1-a)*b*a;
    pzeroalpha(4) = 0;
    pzeroalpha(5) = (1-b)*a^2;
    pzeroalpha(6) = 0;
    pzeroalpha(7) = 1-sum(pzeroalpha,'all');
    
    palphazero(1) = a^2*b;
    palphazero(2) = 0;
    palphazero(3) = (1-a)*b*a;
    palphazero(4) = (1-a)*b*a;
    palphazero(5) = 0;
    palphazero(6) = (1-b)*a^2;
    palphazero(7) = 1-sum(palphazero,'all');
    
    a = exp(-t*alpha^2/8);
    b = exp(-(1-t)*alpha^2/2);
    palphaalpha(1) = a^2*b^2;
    palphaalpha(2) = (1-a)*b^2*a;
    palphaalpha(3) = 0;
    palphaalpha(4) = (1-a)*b^2*a;
    palphaalpha(5) = (1-b)*a^2*b;
    palphaalpha(6) = (1-b)*a^2*b;
    palphaalpha(7) = 1-sum(palphaalpha,'all');
    
    p(1,:) = pzeroalpha;
    p(2,:) = palphazero;
    p(3,:) = palphaalpha;
    
end