function fixed_V1 = fixed_V1(Const_tau,params,N_regions)
    betta = params(1);
    nu = params(3);
    tol = 10^(-16);
    conv = 1;
    V_prev = zeros(N_regions,1);
    V_next = V_prev; 
    iter= 1;
    while conv>tol
        %V_max = max(V_prev);
        exp_beta_v_nu = exp(betta*V_prev/nu);
        Const_tau_sum = Const_tau*exp_beta_v_nu;
        V_next = nu * log(Const_tau_sum); %- V_max;
        conv = sum(abs(V_prev - V_next));
        V_prev = V_next;
        iter = iter +1;
    end
    fixed_V1 = V_next;
end