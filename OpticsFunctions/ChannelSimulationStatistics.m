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
% * q       : Probability of Alice sending a state in the key map basis.
%             Further details given above.
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


function [pzeroalpha, palphazero, palphaalpha, p] = ChannelSimulationStatistics(alpha, q, t)
    
    alpha = abs(alpha);
    czeroalpha = cell(5,1);
    czeroalpha{1} = exp(t*alpha^2/4)-1;
    czeroalpha{2} = exp(t*alpha^2/4)-1;
    czeroalpha{3} = 0;
    czeroalpha{4} = exp((1-t)*alpha^2)-1;
    czeroalpha{5} = 0;
    
    calphazero = cell(5,1);
    calphazero{1} = 0;
    calphazero{2} = exp(t*alpha^2/4)-1;
    calphazero{3} = exp(t*alpha^2/4)-1;
    calphazero{4} = 0;
    calphazero{5} = exp((1-t)*alpha^2)-1;
    
    calphaalpha = cell(5,1);
    calphaalpha{1} = exp(t*alpha^2/8)-1;
    calphaalpha{2} = 0;
    calphaalpha{3} = exp(t*alpha^2/8)-1;
    calphaalpha{4} = exp((1-t)*alpha^2/2)-1;
    calphaalpha{5} = exp((1-t)*alpha^2/2)-1;
    
    a = exp(-alpha^2); % |<0|alpha>|^2
    p{1,1} = q*a*exp(t*alpha^2/2) + (1-q)*a^2*exp(3*t*alpha^2/2);
    pzeroalpha{1,1} = a*exp(t*alpha^2/2);
    palphazero{1,1} = a*exp(t*alpha^2/2);
    palphaalpha{1,1} = a*exp(3*t*alpha^2/4);
    
    for i = 1:1:5
        comboalphazero = nchoosek(calphazero, i);
        combozeroalpha = nchoosek(czeroalpha, i);
        comboalphaalpha = nchoosek(calphaalpha, i);
        l = size(comboalphazero, 1);
        for j = 1:1:l % Loop over all combinations
            tempalphazero = comboalphazero(j,:);
            tempzeroalpha = combozeroalpha(j,:);
            tempalphaalpha = comboalphaalpha(j,:);
            pzeroalpha{j, i+1} = pzeroalpha{1,1};
            palphazero{j, i+1} = palphazero{1,1};
            palphaalpha{j, i+1} = palphaalpha{1,1};
            for k = 1:1:i
                pzeroalpha{j, i+1} = pzeroalpha{j, i+1}*tempzeroalpha{k};
                palphazero{j, i+1} = palphazero{j, i+1}*tempalphazero{k};
                palphaalpha{j, i+1} = palphaalpha{j, i+1}*tempalphaalpha{k};
            end
            p{j, i+1} = 0.5*q*pzeroalpha{j, i+1} + 0.5*q*palphazero{j, i+1} + (1-q)*palphaalpha{j, i+1};
        end
    end
end