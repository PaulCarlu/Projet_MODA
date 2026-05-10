function [w] = solve_omega_K(K,M,uL,ui,uR)
    samples = 500;
    k_range=linspace(-2*pi,2*pi,samples);
    for i=1:length(k_range)
        k = k_range(i);
        lbpos = exp(k*1i);
        lbneg = 1/lbpos;
        [KLL, KLR, KRL, KRR, ~, ~] = condensation_dyn_2(K, uL, ui, uR);
        [MLL, MLR, MRL, MRR, ~, ~] = condensation_dyn_2(M, uL, ui, uR);

        M_cond = lbneg*MRL+MRR+MLL+lbpos*MLR;
        K_cond = lbneg*KRL+KRR+KLL+lbpos*KLR;
        [~,omega]=eig(M_cond,K_cond);
        w(:,i) = diag(omega); 
    end
end

