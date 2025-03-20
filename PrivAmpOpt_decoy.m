function keyrate = PrivAmpOpt_decoy(N, nA, decoySDPcutoff, a, amp, eta, prob, t, iter, rho0)
        
    [~, ~,~, keyrate, ~, ~, ~] = ThreeState_PrivAmp_decoy(N, nA, decoySDPcutoff, a, amp, eta, prob, t, iter, rho0);
%     amp(1)
%     if(dualfvals<0)
%         dualfvals = 0;
%     end
end