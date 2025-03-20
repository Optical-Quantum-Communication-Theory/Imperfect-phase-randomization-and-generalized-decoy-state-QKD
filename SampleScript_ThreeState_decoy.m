clear;
distance = [0:10:300];
% distance = 210;
eta = 10.^(-0.016*distance); % Loss in dB
E = length(eta);
alpha = [0;0.5]; % Decoy amplitudes

N = 1; % Bob's photon number cutoff
nA = 1; % Alice's eigenvector (analogous to photon number) that is being used for key rate
decoySDPcutoff = 17; % Cutoff for decoy-state SDP
dimAshield = 1;
dimA = 3*dimAshield;
dimB = (N+1)*(N+2)/2;
dimTot = dimA*dimB;

a = [0.9128, 0.9564, 1]; % Degree of phase randomisation
iter = 30;
F = length(a);
t = 0.1; % Bob's basis choice beam splitting ratio
prob = 2/3; % Alice's probability of sending Z vs X basis states

rho0 = eye(dimTot)/dimTot;  % replace with better starting points

primalfvals = cell(F,E);
dualfvals = cell(F,E);
gaps = cell(F,E);
flagprimal = cell(F,E);
flagdual = cell(F,E);
keyRate = zeros(F,E);

tstart = tic;

for j = 1:E
     for i = 1:F
        tpointstart = tic;
        optimoptions = optimset('TolX',1e-2,'MaxFunEvals',30, 'Display', 'Iter');  % options field "TolX" for func "fminbnd"
        amp = fminbnd_parallel(@(am)-PrivAmpOpt_decoy(N, nA, decoySDPcutoff, ...
            a(i), [am;alpha], eta(j), prob, t, iter, rho0),0.01,0.99,10,3);
        optamp(i,j) = amp;
        % amp = 0.1189+4*0.0242;

        [Rho{i,j}, primalfvals{i,j}, dualfvals{i,j}, keyRate(i,j), ...
            gaps{i,j}, flagprimal{i,j}, flagdual{i,j}] = ...
            ThreeState_PrivAmp_decoy(N, nA, decoySDPcutoff, a(i), [amp;alpha], eta(j), prob, t, iter, rho0);
        tpointend = toc(tpointstart)
    end
end
tend = toc(tstart);


disp("Time taken = " + string(tend));

save("results");
