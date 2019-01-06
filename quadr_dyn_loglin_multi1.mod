
@#define N_regions = 13

@#for co in 1 : N_regions
var w_@{co} r_@{co} l_@{co} z_@{co} k_@{co} invest_@{co} transf_@{co} fed_transf_@{co} fed_spending_@{co} price_fin_@{co} output_@{co} h_@{co} avg_cons_@{co} v_@{co};
@#for co2 in 1 : N_regions
var mu_@{co}_@{co2} invest_local_@{co}_@{co2} price_@{co}_@{co2} inc_@{co}_@{co2} tau_@{co}_@{co2};
@#for co3 in 1 : N_regions
var cons_@{co}_@{co2}_@{co3};
@#endfor
@#endfor
varexo eps_@{co};
@#endfor
var dividend SDF utilitarian_welfare utilitarian_welfare_pres;

parameters betta betta_diff nu gamma theta l_bar elasticity 
gamma_H sigma_eps rho auto_corr cross_autocorr ppsi 
rich delta  mean_z H S_bar elasticity_w elasticity_h;
@#for co in 1 : N_regions
parameters gamma_@{co}  pphi_@{co};
@#for co2 in 1 : N_regions
parameters d_@{co}_@{co2} tau_real_@{co}_@{co2} tau_util_@{co}_@{co2};
@#endfor
@#endfor
load input;
betta = params_multi(1);
betta_diff = params_multi(2);
nu = params_multi(3);
gamma = params_multi(4);
theta= params_multi(5);
l_bar = params_multi(6);
elasticity = params_multi(7);
gamma_H = params_multi(8);
sigma_eps=  params_multi(9);
rho= params_multi(10);
auto_corr = params_multi(11);
cross_autocorr = params_multi(12);
ppsi = params_multi(13);
rich = 	params_multi(14);
delta = params_multi(15);
mean_z = params_multi(16);
H = params_multi(17);
elasticity_w = params_multi(18);
elasticity_h = params_multi(19);
S_bar = params_multi(20);
@#for co in 1 : N_regions
gamma_@{co} = params_multi(20 + @{co});
pphi_@{co}= params_multi(20 + @{N_regions} + @{co});
@#for co2 in 1 : N_regions
d_@{co}_@{co2} = params_multi(20 + 2 * @{N_regions} + @{co} + (@{co2}-1)*(@{N_regions}));
tau_real_@{co}_@{co2} = params_multi(20 + 2 * @{N_regions} + @{N_regions}*@{N_regions} + @{co} + (@{co2}-1)*(@{N_regions}));
tau_util_@{co}_@{co2} = params_multi(20 + 2 * @{N_regions} + 2 * @{N_regions}*@{N_regions} + @{co} + (@{co2}-1)*(@{N_regions}));
@#endfor
@#endfor    

model;
@#for co in 1 : N_regions
[name = 'Utility']
exp(v_@{co}/nu) = (
@#for ii in 1:N_regions
+exp((betta * v_@{ii}(+1) - tau_@{co}_@{ii})/nu )
@#endfor
);
[name = 'Labor decision']
exp(h_@{co}) = (exp(w_@{co})/(H * exp(price_fin_@{co})))^(1/gamma_H);
[name = 'Price normalization']
@# if co == 1
price_@{co}_@{co} = 0;
@# endif
exp(price_fin_@{co}) = (
@#for ii in 1 : N_regions
+ (exp(price_@{co}_@{ii}) * gamma_@{ii})^(1 - elasticity)
@#endfor
)^(1 / (1 - elasticity));
[name = 'Resource constraint']
@# if co != 1
exp(output_@{co}) = (
@#for ii in 1 : N_regions
+ (exp(invest_local_@{ii}_@{co})+ tau_real_@{ii}_@{co} * exp(mu_@{ii}_@{co})*exp(l_@{ii})) * d_@{ii}_@{co}
+ (
@#for kk in 1 : N_regions
+ (exp(cons_@{ii}_@{co}_@{kk}) * exp(l_@{ii}) * exp(mu_@{ii}_@{kk}) )* d_@{ii}_@{co}
@#endfor
)
@#endfor
); 
@#endif
/*[name = 'Population change']*/
exp(l_@{co}) = (
@#for ii in 1 : N_regions
+ exp(mu_@{ii}_@{co}) * exp(l_@{ii}(-1))
@#endfor
);
[name = 'Competitive wages']
exp(w_@{co}) = (1 - theta)* exp(price_@{co}_@{co})* exp(z_@{co}) * (exp(k_@{co}(-1))/(exp(l_@{co}(-1)) * exp(h_@{co})))^(theta);
[name = 'Competitive rental rate']
r_@{co} = theta * exp(price_@{co}_@{co})* exp(z_@{co}) * (exp(k_@{co}(-1))/(exp(l_@{co}(-1)) * exp(h_@{co})))^(theta - 1);
[name = 'Average consumption']
exp(avg_cons_@{co}) = exp(l_@{co}(-1))/l_bar * exp(inc_@{co}_@{co});
[name = 'Euler equation']
exp(price_fin_@{co}) * (1 + ppsi * (exp(invest_@{co})/exp(k_@{co}(-1)) - delta)) = SDF(+1) * (r_@{co}(+1) + exp(price_fin_@{co}(+1)) * ( ppsi/2 * (exp(invest_@{co}(+1))/exp(k_@{co}) - delta)^2 + 1- delta
+ ppsi * (exp(invest_@{co}(+1))/exp(k_@{co}) - delta)) );
[name = 'Investment']
exp(k_@{co}) = exp(invest_@{co})+ (1-delta)*exp(k_@{co}(-1));
[name = 'Productivity']
z_@{co} = auto_corr*z_@{co}(-1) +eps_@{co} + mean_z* (1 - auto_corr);
[name = 'CapitalTransfers']
transf_@{co} = pphi_@{co} * dividend/exp(l_@{co}(-1));
[name = 'FedTransfers']
exp(fed_transf_@{co}) = S_bar  *(exp(elasticity_w * (w_@{co} - steady_state(w_@{co}))))*(exp(elasticity_h * (h_@{co} - steady_state(h_@{co}))))*
(
@#for ii in 1 : N_regions
+ pphi_@{ii} * exp(output_@{ii})
@#endfor
)
;
[name = 'FedSpending']
exp(fed_spending_@{co}) = pphi_@{co} *
(
@#for ii in 1 : N_regions
+ exp(fed_transf_@{ii})
@#endfor
);
[name = 'Output']
exp(output_@{co}) = exp(z_@{co}) * exp(k_@{co}(-1))^(theta) * (exp(h_@{co}) * exp(l_@{co}(-1)))^(1 - theta);
@#for co2 in 1 : N_regions
[name = 'Income']
exp(inc_@{co}_@{co2}) = (exp(w_@{co}) * exp(h_@{co}) + transf_@{co} - exp(price_fin_@{co}) * exp(fed_transf_@{co})/exp(l_@{co}(-1)) - exp(price_@{co}_@{co2}) *tau_real_@{co}_@{co2})/exp(price_fin_@{co});
[name = 'Resistance to migration']
tau_@{co}_@{co2} = - (exp(inc_@{co}_@{co2}) - H* exp(h_@{co})^(1 + gamma_H)/(1 + gamma_H) )^(1 - gamma)/(1 - gamma) + tau_util_@{co}_@{co2};
[name = 'Migration rates']
exp(mu_@{co}_@{co2}) = exp((betta * v_@{co2}(+1) - tau_@{co}_@{co2})/nu) / exp(v_@{co}/nu);
[name = 'Local investment']
exp(invest_local_@{co}_@{co2}) = exp(price_@{co}_@{co2})^(-elasticity)*gamma_@{co2}^(1-elasticity)*(exp(invest_@{co}) + exp(fed_spending_@{co}) + ppsi/2 * exp(k_@{co}(-1)) * (exp(invest_@{co})/exp(k_@{co}(-1)) - delta)^2)*exp(price_fin_@{co})^(elasticity);
[name = 'Local prices']
@# if co != co2
exp(price_@{co}_@{co2}) = d_@{co}_@{co2} *exp(price_@{co2}_@{co2});
@# endif
@#for co3 in 1 : N_regions
exp(cons_@{co}_@{co2}_@{co3}) =  exp(inc_@{co}_@{co3})* exp(price_@{co}_@{co2})^(- elasticity)* exp(price_fin_@{co})^(elasticity)* gamma_@{co}^(1 - elasticity);
@#endfor
@#endfor
@#endfor 
[name = 'Dividends']
dividend = (
@#for ii in 1 : N_regions
+(r_@{co} * exp(k_@{co}(-1)) - exp(price_fin_@{co}) * exp(invest_@{co}))
@#endfor
);
[name = 'Welfare']
utilitarian_welfare_pres = (
@#for ii in 1 : N_regions
@#for jj in 1 : N_regions
+ -exp(mu_@{ii}_@{jj}) * exp(l_@{ii}(-1)) * (tau_@{ii}_@{jj} - tau_util_@{ii}_@{jj})
@#endfor
@#endfor
);
utilitarian_welfare = (
@#for ii in 1 : N_regions
+  exp(l_@{ii}(-1)) *v_@{ii}
@#endfor
);
[name = 'SDF']
SDF = (betta + betta_diff) * (
@#for ii in 1 : N_regions
+(exp(avg_cons_@{co}) /exp(avg_cons_@{co}(-1)))^(-gamma) *pphi_@{co}
@#endfor
);
end;

initval;
@#for co in 1 : N_regions
w_@{co} = xsteady_multi(@{co});
r_@{co} = xsteady_multi(@{N_regions} + @{co});
l_@{co} = xsteady_multi(2 * @{N_regions} + @{co});
z_@{co} = xsteady_multi(3 * @{N_regions} + @{co});
k_@{co} = xsteady_multi(4 * @{N_regions} + @{co});
invest_@{co} = xsteady_multi(5 * @{N_regions} + @{co});
transf_@{co} = xsteady_multi(6 * @{N_regions} + @{co});
price_fin_@{co} = xsteady_multi(7 * @{N_regions} + @{co});
output_@{co} = xsteady_multi(8 * @{N_regions} + @{co});
h_@{co} = xsteady_multi(9 * @{N_regions} + @{co});
avg_cons_@{co} = xsteady_multi(10 * @{N_regions} + @{co});
v_@{co} = xsteady_multi(11 * @{N_regions} + @{co});
fed_transf_@{co} = xsteady_multi(12 * @{N_regions} + @{co});
fed_spending_@{co} = xsteady_multi(12 * @{N_regions} + @{co});
@#for co2 in 1 : N_regions
mu_@{co}_@{co2} = xsteady_multi(13 * @{N_regions} + @{co} + (@{co2}-1)*(@{N_regions}));
invest_local_@{co}_@{co2} = xsteady_multi(13 * @{N_regions}+ @{N_regions}* @{N_regions} + @{co} + (@{co2}-1)*(@{N_regions}));
price_@{co}_@{co2} = xsteady_multi(13 * @{N_regions}+ 2 * @{N_regions}* @{N_regions} + @{co} + (@{co2}-1)*(@{N_regions}));
inc_@{co}_@{co2} = xsteady_multi(13 * @{N_regions}+ 3 * @{N_regions}* @{N_regions} + @{co} + (@{co2}-1)*(@{N_regions}));
tau_@{co}_@{co2} = xsteady_multi(13 * @{N_regions}+ 4 * @{N_regions}* @{N_regions} + @{co} + (@{co2}-1)*(@{N_regions}));
@#for co3 in 1 : N_regions
cons_@{co}_@{co2}_@{co3} = xsteady_multi(13 * @{N_regions}+ 5 * @{N_regions}* @{N_regions} + @{co} + (@{co2}-1)*(@{N_regions})  + (@{co3}-1)*(@{N_regions}*@{N_regions}));
@#endfor
@#endfor
@#endfor
dividend = xsteady_multi(13 * @{N_regions}+ 5 * @{N_regions}* @{N_regions} + @{N_regions}*@{N_regions}*@{N_regions} + 1);
SDF = xsteady_multi(13 * @{N_regions}+ 5 * @{N_regions}* @{N_regions} + @{N_regions}*@{N_regions}*@{N_regions} + 2);
utilitarian_welfare = xsteady_multi(13 * @{N_regions}+ 5 * @{N_regions}* @{N_regions} + @{N_regions}*@{N_regions}*@{N_regions} + 3);
utilitarian_welfare_pres = xsteady_multi(13 * @{N_regions}+ 5 * @{N_regions}* @{N_regions} + @{N_regions}*@{N_regions}*@{N_regions} + 4);
end;
shocks;
@#for co in 1 : N_regions
var eps_@{co}; stderr sigma_eps;
@#endfor
@#for co in 1 : N_regions
@#for co2 in co : N_regions
@# if co != co2
corr eps_@{co},eps_@{co2} = rho;
@# endif
@#endfor
@#endfor
end;
steady(maxit = 10, solve_algo = 1);

check(qz_zero_threshold=1e-10);

stoch_simul(order = 1, irf = 0, periods = 10003,nograph,nocorr,nofunctions,nomoments,noprint); /*hp_filter = 1600,*/