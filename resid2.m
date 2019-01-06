function res = resid2(x,params,N_regions)
    %Parameters:
    betta = params(1);
    betta_diff = params(2);
    nu = params(3);
    gamma = params(4);
    theta= params(5);
    l_bar = params(6);
    elasticity = params(7);

    gamma_H = params(8);
    
    sigma_eps=  params(9);
    rho= params(10);
    auto_corr = params(11);
    cross_autocorr = params(12);
    ppsi = params(13);
    rich = 	params(14);
    delta = params(15);
    mean_z = params(16);
    H = params(17);
    elasticity_w = params(18);
    elasticity_h = params(19);
    S_bar = params(20);
    
    gamma_mat = params(21:(21+N_regions - 1));
    %tau_util_mat = params((21+N_regions):(21+2*N_regions-1));
    pphi = params((21+N_regions):(21+2*N_regions-1));
    d =  params((21+2*N_regions):(21+2*N_regions + N_regions*N_regions - 1));
    d = reshape(d,[N_regions,N_regions]);
    tau_real = params((21+2*N_regions + N_regions*N_regions):(21+2*N_regions + 2* N_regions*N_regions -1));
    tau_real = reshape(tau_real,[N_regions,N_regions]);
    tau_util_mat = params((21+2*N_regions + 2*N_regions*N_regions):(21+2*N_regions + 3* N_regions*N_regions -1));
    tau_util_mat = reshape(tau_util_mat,[N_regions,N_regions]);
    lower_bound = 10^(-8);
	res = ones(2 * (N_regions-1),1);
    price_mat = zeros(N_regions,N_regions);
    l = ones(N_regions,1);
    %v = zeros(N_regions,1);
    price_mat(1,1) = 1;
    for ii = 2:(N_regions)
        price_mat(ii,ii) = x(ii-1);
        l(ii) = x(ii - 1 + N_regions -1);
        %v(ii) = x(ii - 1 + 2*N_regions - 2);
    end

    l(1) = l_bar - sum(l(2:N_regions));
    price_mat = d * price_mat;
	diag_price_mat = diag(price_mat);
    SDF_ss = betta;
    price_fin = (price_mat.^(1 - elasticity) *  gamma_mat.^(1 - elasticity)).^(1/(1 - elasticity));
    r =  price_fin*(1/betta - 1 + delta);
    K_Lh_ratio = (r./diag_price_mat/theta).^(1 / (theta - 1)); 
    w =  (1 - theta) * diag_price_mat.* K_Lh_ratio.^(theta);
    h = (w./(H* price_fin)).^(1/gamma_H);
    k = K_Lh_ratio.* h.*l;
    output = k.^(theta).*(l.* h).^( 1 - theta);
    fed_transf = S_bar *output;
    invest = delta * k;
    dividend = sum(r.* k - price_fin .* invest);
    transf = pphi .* dividend./l;
    
    inc = (repmat(w.* h+ transf - price_fin.* fed_transf./l,1,N_regions) - tau_real.*price_mat)./repmat(price_fin,1,N_regions);
    tau = - (inc - repmat(H* h.^(1 + gamma_H)/(1 + gamma_H),1,N_regions) ).^(1-gamma)/(1-gamma);
    tau = tau + tau_util_mat;
    Const_tau = exp(-tau/nu);
    V = fixed_V1(Const_tau,params,N_regions);
    exp_beta_v_nu = exp(betta*V/nu);
    exp_beta_v_nu_mat = repmat(exp_beta_v_nu,1,N_regions);
    Const_tau_id = Const_tau.*exp_beta_v_nu_mat;
    V_mat = repmat(exp(V/nu),1,N_regions);
    
    mu = Const_tau_id./V_mat;
    res(N_regions: 2*(N_regions - 1)) = l(2:N_regions) - mu(:,2:N_regions)' * l;
    
    cons = zeros(N_regions,N_regions,N_regions);
    invest_local = zeros(N_regions,N_regions);
    for ii = 1:N_regions
        for jj = 1:N_regions
            invest_local(ii,jj) = price_mat(ii,jj)^(-elasticity)*gamma_mat(jj)^(1-elasticity) * (invest(ii) + fed_transf(ii))* price_fin(ii)^(elasticity);
            for kk = 1:N_regions
                cons(ii,jj,kk) = inc(ii,kk)*price_fin(ii)^(elasticity)*price_mat(ii,jj)^(- elasticity) * gamma_mat(jj)^(1 - elasticity);
            end
        end
    end
    totalcons = zeros(N_regions,1);
    for jj = 1:N_regions
        for ii = 1:N_regions
            totalcons(jj) = totalcons(jj) + (invest_local(ii,jj)+ tau_real(ii,jj) * mu(ii,jj)*l(ii))*d(ii,jj);
            for kk = 1:N_regions
                totalcons(jj) = totalcons(jj) + ((cons(ii,jj,kk))*l(ii)*mu(ii,kk))*d(ii,jj);
            end
        end        
    end
    res(1:N_regions-1) = output(2:N_regions)-totalcons(2:N_regions);
    %res = mu;
%     price_11_ss = 1;
%     price_21_ss = d * price_11_ss;
%     price_12_ss = d * price_22_ss;
%     price_fin_1_ss = ((gamma_1 * price_11_ss)^(1 - elasticity) + (gamma_2 * price_12_ss)^(1 - elasticity))^(1/(1 - elasticity));
%     price_fin_2_ss = ((gamma_1 * price_21_ss)^(1 - elasticity) + (gamma_2 * price_22_ss)^(1 - elasticity))^(1/(1 - elasticity));
%     r_ss = price_fin_1_ss*(1/betta - 1 + delta);
%     r_f_ss = price_fin_2_ss*(1/betta - 1 + delta);
% 	K_Lh_ratio = (r_ss/rich/price_11_ss/theta)^(1 / (theta - 1)); 
% 	K_Lh_ratio_f = (r_f_ss/price_22_ss/theta)^(1 / (theta - 1)); 
% 	w_ss = price_11_ss *rich * (1 - theta) * K_Lh_ratio^(theta);
% 	w_f_ss = price_22_ss * (1 - theta) * K_Lh_ratio_f^(theta);
% 	h_1_ss = (w_ss/(H* price_fin_1_ss))^(1/gamma_H);
% 	h_2_ss = (w_f_ss/(H * price_fin_2_ss))^(1/gamma_H);
% 	k_ss = K_Lh_ratio * h_1_ss *(l_ss);
% 	k_f_ss = K_Lh_ratio_f * h_2_ss *(l_f_ss);
% 	invest_ss = delta * k_ss;
% 	invest_f_ss = delta * k_f_ss;
% 	dividend_ss = r_ss * k_ss + r_f_ss * k_f_ss - price_fin_1_ss * invest_ss - price_fin_2_ss * invest_f_ss;
% 	transf_1_ss = pphi * dividend_ss/l_ss;
% 	transf_2_ss = (1 - pphi) * dividend_ss/l_f_ss;
% 	inc_11_ss = max((w_ss * h_1_ss + transf_1_ss)/price_fin_1_ss,lower_bound);
% 	inc_12_ss = max((w_ss * h_1_ss + transf_1_ss - price_12_ss * tau_real)/price_fin_1_ss,lower_bound);
% 	inc_21_ss = max((w_f_ss * h_2_ss + transf_2_ss - price_21_ss * tau_real)/price_fin_2_ss,lower_bound);
% 	inc_22_ss = max((w_f_ss * h_2_ss + transf_2_ss)/price_fin_2_ss,lower_bound);
% 	tau_11_ss = 0;
% 	tau_12_ss = -inc_12_ss^(1 - gamma)/(1 - gamma) +inc_11_ss^(1 - gamma)/(1 - gamma)+tau_util_2;
% 	tau_21_ss = -inc_21_ss^(1 - gamma)/(1 - gamma) + inc_22_ss^(1 - gamma)/(1 - gamma)+tau_util_1;
% 	tau_22_ss = 0;
% 	Const_tau_11 = max(exp(-tau_11_ss/nu),lower_bound);
% 	Const_tau_12 = max(exp(-tau_12_ss/nu),lower_bound);
% 	Const_tau_21 = max(exp(-tau_21_ss/nu),lower_bound);
% 	Const_tau_22 = max(exp(-tau_22_ss/nu),lower_bound);
% 	val_x =max(((exp(v_f_ss/nu) - exp(betta * v_f_ss/nu)*Const_tau_22)/Const_tau_21),lower_bound);
% 	v_ss = nu * log(val_x^(1/betta));
% 	res(3) = val_x^(1/betta) - val_x * Const_tau_11 - exp(betta * v_f_ss/nu)* Const_tau_12;
% 	mu_ss = exp(betta * v_f_ss/nu)*Const_tau_12/ (exp(betta * v_f_ss/nu)*Const_tau_12 + exp(betta * v_ss/nu)*Const_tau_11);
% 	mu_f_ss = exp(betta * v_ss/nu)*Const_tau_21/ (exp(betta * v_ss/nu)*Const_tau_21 + exp(betta * v_f_ss/nu)*Const_tau_22);
% 	res(2) = mu_f_ss * l_f_ss - mu_ss * l_ss;
% 	cons_11_ss = inc_11_ss * price_11_ss^(- elasticity) * gamma_1^(1 - elasticity)* price_fin_1_ss^(elasticity);
% 	cons_12_ss = inc_11_ss * price_12_ss^(- elasticity) * gamma_2^(1 - elasticity)* price_fin_1_ss^(elasticity);
% 	cons_21_ss = inc_22_ss * price_21_ss^(- elasticity) * gamma_1^(1 - elasticity) * price_fin_2_ss^(elasticity);
% 	cons_22_ss = inc_22_ss * price_22_ss^(- elasticity) * gamma_2^(1 - elasticity) * price_fin_2_ss^(elasticity);
% 	cons_11_m_ss = inc_12_ss * price_11_ss^(- elasticity) * gamma_1^(1 - elasticity)* price_fin_1_ss^(elasticity);
% 	cons_12_m_ss = inc_12_ss * price_12_ss^(- elasticity) * gamma_2^(1 - elasticity)* price_fin_1_ss^(elasticity);
% 	cons_21_m_ss = inc_21_ss * price_21_ss^(- elasticity) * gamma_1^(1 - elasticity) * price_fin_2_ss^(elasticity);
% 	cons_22_m_ss = inc_21_ss * price_22_ss^(- elasticity) * gamma_2^(1 - elasticity) * price_fin_2_ss^(elasticity);
% 	invest_11_ss = price_11_ss^(-elasticity)*gamma_1^(1-elasticity)*invest_ss* price_fin_1_ss^(elasticity);
% 	invest_12_ss = price_12_ss^(-elasticity)*gamma_2^(1-elasticity)*invest_ss* price_fin_1_ss^(elasticity);
% 	invest_21_ss = price_21_ss^(-elasticity)*gamma_1^(1-elasticity)*invest_f_ss * price_fin_2_ss^elasticity;
% 	invest_22_ss = price_22_ss^(-elasticity)*gamma_2^(1-elasticity)*invest_f_ss * price_fin_2_ss^elasticity;
% 	totalcons2_ss = (cons_12_ss * l_ss * (1 - mu_ss) + invest_12_ss) * d  +  (cons_22_ss*(1 -  mu_f_ss)*l_f_ss+ invest_22_ss)  + (cons_12_m_ss + tau_real) * d * mu_ss * l_ss  +  cons_22_m_ss* mu_f_ss * l_f_ss;
%     res(1) = totalcons2_ss -  k_f_ss^(theta)*(l_f_ss * h_2_ss)^( 1 - theta);
end