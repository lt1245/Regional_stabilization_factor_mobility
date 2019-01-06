%cd 'D:\NYU\3rd year\project\Model2\Europe_US'
%addpath D:\Programs\Dynare4.5.1\matlab
addpath 'C:\Program Files (x86)\dynare4.5.3\matlab'
clear; 
close all;
clc;
%% Parameters
N_regions = 13;
betta = 0.99;
betta_diff = 0.00;
nu = 2;
gamma = 2.0; %/*of substitution across time*/
theta= 0.36;
l_bar = 1;
elasticity = 3.0;% /*of substitution across goods*/
% gamma_1 = 1;
% gamma_2 = 1;
gamma_mat = ones(N_regions,1);
gamma_H = 1;
%d = 2.5; %/*iceberg cost*/
d = 1.7* ones(N_regions,N_regions);
sigma_eps=  0.2^2;
rho= 0.265;
auto_corr = 0.906;
cross_autocorr = 0.00;%/*0.088;*/
pphi = ones(N_regions,1)/N_regions; %/*ownership in the capitalist firm*/
ppsi = 2.05;%10;EU % US:2.05;
rich = 	1;
delta = 0.025;
mean_z = 0;
H = (1 - theta) * ((1/betta - 1 + delta)/theta )^(theta/(theta - 1));
% tau_real = 1.6;
tau_real = 1.6* ones(N_regions ,N_regions);
for ii = 1:N_regions
    tau_real(ii,ii) = 0;
    d(ii,ii) = 1;
end
% tau_util_1 = 10.1;
% tau_util_2 = 10.1;%8.99;
tau_util_mat = 14 * ones(N_regions,1);
tau_util_mat = repmat(tau_util_mat',N_regions,1);
for ii = 1:N_regions
    tau_util_mat(ii,ii) = 0;
end
%elasticity_w = -0.9;
elasticity_w = 0;
%elasticity_w = 0.9;
%elasticity_w = 1.15;
%elasticity_h = -1.6;
elasticity_h = 0;
%elasticity_h = 1.6;
S_bar = 0.2;
lower_bound = 10^(-8); % numerical imprecision
params = [betta; betta_diff; nu; gamma; theta; l_bar; elasticity; 
gamma_H;  sigma_eps; rho; auto_corr; cross_autocorr; ppsi; 
rich; delta; mean_z;H;elasticity_w; elasticity_h; S_bar; gamma_mat; pphi;reshape(d,[N_regions*N_regions,1]);reshape(tau_real,[N_regions*N_regions,1]);reshape(tau_util_mat,[N_regions*N_regions,1])];
%% Steady state calculations

x_start = [1*ones(N_regions-1,1);ones(N_regions-1,1)/N_regions];%[1,0.5,1.09413]';
resid1 = @(y)resid2(y,params,N_regions);
options = optimoptions('fsolve','Display','iter','MaxIter', 400,'MaxFunEvals',10000,'FunValCheck','on');
[x_ss,fval,exitflag]  = fsolve(resid1, x_start,options);


res = ones(2 * (N_regions-1),1);
price_mat = zeros(N_regions,N_regions);
l = ones(N_regions,1);
%v = zeros(N_regions,1);
price_mat(1,1) = 1;
for ii = 2:(N_regions)
    price_mat(ii,ii) = x_ss(ii-1);
    l(ii) = x_ss(ii - 1 + N_regions -1);
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

tau = - (inc - repmat(H* h.^(1 + gamma_H)/(1 + gamma_H),1,N_regions)).^(1-gamma)/(1-gamma);
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
avg_cons = l .* diag(inc)/l_bar;
utilitarian_welfare = 0;
for ii = 1:N_regions
    utilitarian_welfare = utilitarian_welfare +  l(ii) *V(ii);
    %for jj = 1:N_regions
        %utilitarian_welfare = utilitarian_welfare - mu(ii,jj) * l(ii) * (tau(ii,jj) - tau_util_mat(ii,jj));
    %end
end
tau_welf = -(tau - tau_util_mat);
utilitarian_welfare_pres = 0;
for ii = 1:N_regions
    for jj = 1:N_regions
        utilitarian_welfare_pres = utilitarian_welfare_pres + mu(ii,jj) * l(ii) * tau_welf(ii,jj);
    end
end
disp(['--------------------------------------------'])
disp(['SS observables            | ', 'Mean '])
disp(['--------------------------------------------'])
disp(['Employment | ', num2str([mean(l) ])])
disp(['Average annual migration rate| ', num2str([sum(mu'*l - diag(mu).*l) ])])
disp(['---------------------------------------------'])
%% Save results
params_multi = params;
xsteady_multi   = [log(w); r; log(l); zeros(N_regions,1); log(k); log(invest);transf;log(price_fin);log(output); log(h); log(avg_cons);V;log(fed_transf);
log(mu(:)) ;log(invest_local(:)); log(price_mat(:)); log(inc(:)); tau(:); 
log(cons(:)); dividend;  SDF_ss;utilitarian_welfare;utilitarian_welfare_pres];
save input params_multi xsteady_multi 
% 
% dynare quadr_dyn_loglin_base
dynare quadr_dyn_loglin_multi1
% d3 = zeros(N_regions,N_regions);
% for co = 1:N_regions
%     for co2 = 1:N_regions
%         d3(co,co2) = params(17 + 2 * N_regions + 2 * N_regions*N_regions + co + (co2-1)*(N_regions));
% %         for co3 = 1:N_regions
% %             d3(co,co2,co3) = xsteady_multi(12 * N_regions+ 5 * N_regions* N_regions + co + (co2-1)*(N_regions)  + (co3-1)*(N_regions*N_regions));
% %         end
%     end
% end
N_regions = 13;
sample_size_total = size(avg_cons_1,1) - 3;
time_length = 26;
burn_period = 100 - time_length;
for i = 1:N_regions
   cons_pc(:,i) = eval(['inc_',num2str(i),'_',num2str(i)]); % log per capita cons but as in the data, assuming migrants consume the same amount
   output(:,i) = eval(['output_',num2str(i)]) ; % log output
   pop(:,i) = eval(['l_',num2str(i)]) ; % log pop
   fed_transf(:,i) = eval(['fed_transf_',num2str(i)]) ; % log fed transf
   fed_spending(:,i) = eval(['fed_spending_',num2str(i)]) ; % log fed spending
   transf(:,i) = eval(['transf_',num2str(i)]) ; % capital income transf
   k(:,i) = eval(['k_',num2str(i)]) ; % log capital stock
   r(:,i) = eval(['r_',num2str(i)]) ; % rental rate
   invest(:,i) = eval(['invest_',num2str(i)]) ; % log investment
   diag_prices(:,i) = eval(['price_',num2str(i),'_',num2str(i)]) ; % log prices of regional goods
   price_fin(:,i) = eval(['price_fin_',num2str(i)]) ; % log final prices 
end

%% Construct variables for the regression
gdp =  output(3:end,:) + diag_prices(3:end,:) - price_fin(3:end,:) - pop(2:(end-1),:);
gdp_bar = output(3:end,:) + diag_prices(3:end,:) - price_fin(3:end,:) - pop(1:(end-2),:);
net_capital_income = (transf(3:end,:) - (exp(k(2:(end-1),:)) .* r(3:end,:)  - exp(price_fin(3:end,:)).*exp(invest(3:end,:)))./ exp(pop(2:(end-1),:))) ./exp(price_fin(3:end,:));
gni = log(exp(gdp) + net_capital_income );
net_fed_transfers =  (exp(fed_transf(3:end,:)) - exp(fed_spending(3:end,:)))./ exp(pop(2:(end-1),:))./exp(price_fin(3:end,:));
di = log(exp(gni) - net_fed_transfers);
cons_pc_real = log(exp(cons_pc(3:end,:)) + exp(fed_spending(3:end,:) - price_fin(3:end,:) - pop(2:(end-1),:))) ;
% gdp =  output(3:end,:) + diag_prices(3:end,:) - price_fin(3:end,:) - pop(3:(end-0),:);
% gdp_bar = output(3:end,:) + diag_prices(3:end,:) - price_fin(3:end,:) - pop(2:(end-1),:);
% net_capital_income = (transf(3:end,:) - (exp(k(2:(end-1),:)) .* r(3:end,:)  - exp(price_fin(3:end,:)).*exp(invest(3:end,:)))./ exp(pop(3:(end-0),:))) ./exp(price_fin(3:end,:));
% gni = log(exp(gdp) + net_capital_income);
% net_fed_transfers =  (exp(fed_transf(3:end,:)) - exp(fed_spending(3:end,:)))./ exp(pop(3:(end),:))./exp(price_fin(3:end,:));
% di = log(exp(gni) + net_fed_transfers);
% cons_pc_real = log(exp(cons_pc(3:end,:)) + exp(fed_spending(3:end,:) - price_fin(3:end,:) - pop(3:(end),:))) ;


d_gdp = gdp(2:end,:) - gdp(1:(end-1),:);
d_gdp_bar = gdp_bar(2:end,:) - gdp_bar(1:(end-1),:);
d_gni = gni(2:end,:) - gni(1:(end-1),:);
d_di = di(2:end,:) - di(1:(end-1),:);
d_di2 = cons_pc(4:end,:) - cons_pc(3:(end-1),:);
d_gdp_vec = reshape(d_gdp',N_regions * sample_size_total,1);
d_gdp_bar_vec = reshape(d_gdp_bar',N_regions * sample_size_total,1);
d_gni_vec = reshape(d_gni',N_regions * sample_size_total,1);
d_di_vec = reshape(d_di',N_regions * sample_size_total,1);
d_di2_vec = reshape(d_di2',N_regions *sample_size_total,1);

d_gdp_vec1 = d_gdp_vec - mean(d_gdp_vec); 
d_gdp_bar_vec1 = d_gdp_bar_vec - mean(d_gdp_bar_vec); 
d_gni_vec1 = d_gni_vec - mean(d_gni_vec); 
d_di_vec1 = d_di_vec - mean(d_di_vec);
d_di2_vec1 = d_di2_vec - mean(d_di2_vec); 

time_fixed_effects = kron(eye(time_length,time_length),ones(N_regions,1));
eta1_sum = zeros(4,1);
for i = 1: sample_size_total/(100)
    start_period = burn_period + (i-1)*100;

    X = d_gdp_bar_vec1((start_period * N_regions + 1):((time_length + start_period)  * N_regions),1);
    X = X -time_fixed_effects * inv(time_fixed_effects'*time_fixed_effects) * time_fixed_effects' * X;
    X = X - mean(X);

    y_1 = d_gni_vec1((start_period * N_regions + 1):((time_length + start_period)  * N_regions),1);
    y_1 = y_1 -time_fixed_effects * inv(time_fixed_effects'*time_fixed_effects) * time_fixed_effects' * y_1;
    y_1 = y_1 - mean(y_1);
    res_1 = X * inv(X'*X) * X' * y_1 - y_1;
    w_1 = res_1.^2;
    w_1_diag = w_1;
    for ii = 1:N_regions
        w_1_diag(1:N_regions:end) = sqrt(1 / mean(w_1(1:N_regions:end)));
    end
    X_1 = diag(w_1_diag) * X;
    y_1_g = diag(w_1_diag) * y_1;
    res_1g = diag(w_1_diag) * res_1;

    y_2 = d_di_vec1((start_period * N_regions + 1):((time_length + start_period)  * N_regions),1);
    y_2 = y_2 -time_fixed_effects * inv(time_fixed_effects'*time_fixed_effects) * time_fixed_effects' * y_2;
    y_2 = y_2 - mean(y_2);
    res_2 = X * inv(X'*X) * X' * y_2 - y_2;
    w_2 = res_2.^2;
    w_2_diag = w_2;
    for ii = 1:N_regions
        w_2_diag(1:N_regions:end) = sqrt(1 / mean(w_2(1:N_regions:end)));
    end
    X_2 = diag(w_2_diag) * X;
    y_2_g = diag(w_2_diag) * y_2;
    res_2g = diag(w_2_diag) * res_2;

    y_3 = d_di2_vec1((start_period * N_regions + 1):((time_length + start_period)  * N_regions),1);
    y_3 = y_3 -time_fixed_effects * inv(time_fixed_effects'*time_fixed_effects) * time_fixed_effects' * y_3;
    y_3 = y_3 - mean(y_3);
    res_3 = X * inv(X'*X) * X' * y_3 - y_3;
    w_3 = res_3.^2;
    w_3_diag = w_3;
    for ii = 1:N_regions
        w_3_diag(1:N_regions:end) = sqrt(1 / mean(w_3(1:N_regions:end)));
    end
    X_3 = diag(w_3_diag) * X;
    y_3_g = diag(w_3_diag) * y_3;
    res_3g = diag(w_3_diag) * res_3;

    y_4 = d_gdp_vec1((start_period * N_regions + 1):((time_length + start_period)  * N_regions),1);
    y_4 = y_4 -time_fixed_effects * inv(time_fixed_effects'*time_fixed_effects) * time_fixed_effects' * y_4;
    y_4 = y_4 - mean(y_4);
    res_4 = X * inv(X'*X) * X' * y_4 - y_4;
    w_4 = res_4.^2;
    w_4_diag = w_4;
    for ii = 1:N_regions
        w_4_diag(1:N_regions:end) = sqrt(1 / mean(w_4(1:N_regions:end)));
    end
    X_4 = diag(w_4_diag) * X;
    y_4_g = diag(w_4_diag) * y_4;
    res_4g = diag(w_4_diag) * res_4;

    %% FGLS - SUR
    X_gls = blkdiag(X_1,X_2,X_3,X_4);
    y_stack = [y_1_g;y_2_g;y_3_g;y_4_g];
    Omega_sur = ones(4,4);
    Omega_sur(1,2) = mean(res_1g .* res_2g); 
    Omega_sur(1,3) = mean(res_1g .* res_3g);
    Omega_sur(1,4) = mean(res_1g .* res_4g);
    Omega_sur(4,1) = Omega_sur(1,4);
    Omega_sur(3,1) = Omega_sur(1,3);
    Omega_sur(2,1) = Omega_sur(1,2);


    Omega_sur(2,3) = mean(res_3g .* res_2g); 
    Omega_sur(2,4) = mean(res_4g .* res_2g); 

    Omega_sur(3,2) = Omega_sur(2,3);
    Omega_sur(4,2) = Omega_sur(2,4);


    Omega_sur(3,4) = mean(res_4g .* res_3g); 

    Omega_sur(4,3) = Omega_sur(3,4);

    Omega_sur_inv = inv(Omega_sur);
    Omega_inv = kron(Omega_sur_inv, eye(time_length * N_regions,time_length  * N_regions));

    beta = inv(X_gls' * Omega_inv *  X_gls) * X_gls' * Omega_inv  * y_stack;

    eta = zeros(4,4);
    eta(1,4) = -1;
    eta(2,4) = 1;
    eta(2,1) = -1;
    eta(3,1) = 1;
    eta(3,2) = -1;
    eta(4,2) = 1;
    eta(4,3) = 0;
    eta1 = eta * beta;
    eta1(1,1) = eta1(1,1) + 1; 
    eta1_sum = eta1 + eta1_sum;
end
eta1_sum/(sample_size_total/(100))
((1 - betta) * (1 - gamma) * xsteady_multi(end-1))^(1 / (1 - gamma)) 
mean(((1 - betta) * (1 - gamma) * utilitarian_welfare).^(1 / (1 - gamma)))
((1 - betta) * (1 - gamma) * mean(utilitarian_welfare)).^(1 / (1 - gamma))
mean(var(cons_pc_real) * gamma / 2)
mean(mean(cons_pc_real))
mean(utilitarian_welfare)
Raws_util = ones(10003,N_regions);
for i=1:N_regions
  eval(sprintf('Raws_util(:,i) = v_%d;', i));
end
mean(min(Raws_util')')
mean(utilitarian_welfare_pres)
%% EU baseline:
% CE: 24.5359 - steady
% CE: 24.5321 - business cycle
% Baseline cost of BSC - 1.87%
% Lucas measure: 0.0122 - volatility * gamma
% Mean real consumption: 3.7277
%mean(utilitarian_welfare) -4.0236
%Rawlsian welfare: -4.0962
%mean(utilitarian_welfare_pres)-0.0625
%% EU federal govt
% CE: 24.5359 - steady
% CE: 22.8102 - business cycle
% Cost of BSC - 0.07%
%   -0.0002
%    0.0747
%    0.2537
%    0.6718
% Lucas measure: 0.0140 - volatility * gamma
% Mean real consumption: 3.7280
%mean(utilitarian_welfare)  -4.0814
%Rawlsian welfare:-4.1371
%% EU better inv
% CE: 24.5359 - steady
% CE: 23.4734 - business cycle
% Cost of BSC - 0.9567 better?
%     0.0006
%     0.3526
%    -0.0000
%     0.6468
% Lucas measure: 0.0168 - volatility * gamma
% Mean real consumption: 3.7350
%mean(utilitarian_welfare) -4.3047
%mean(utilitarian_welfare_pres) -0.0658
%% EU higher migration
% CE: 217.7825 - steady
% CE: 76.3713 - business cycle
% Cost of BSC - 0.9567 better?
%   -0.0004
%    0.0686
%   -0.0000
%    0.9318
% Lucas measure: 0.0120 - volatility * gamma
% Mean real consumption: 3.7278
%mean(utilitarian_welfare) -0.3950