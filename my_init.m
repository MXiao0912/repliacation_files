

v2struct(SettingsTransfor)

%% SET INITIAL VALUES TIME VARYING ELEMENTS 
% f = [ln(sigma_r), ln(sigma_e), ln(sigma_o), atanh(pi_eo)]



%% Initialize 
dp = diff(All_Data(:,1));
[acor, lag] = xcorr(dp);
sigma_e2 = max(0, -acor(lag==1));
sigma_r2_1 = max(0, acor(lag==0)-2*sigma_e2);

dn = diff(All_Data(:,2));
[acor, lag] = xcorr(dn);
sigma_w2 = max(0, -acor(lag==1));
sigma_r2_2 = max(0, acor(lag==0)-2*sigma_w2);

sigma_r2 = 0.5*(sigma_r2_1+sigma_r2_2);

cov_eo = cov(dp,dn);
rho_eo = max(min((cov_eo(1,2)-sigma_r2)/2,sqrt(sigma_e2*sigma_w2)), -sqrt(sigma_e2*sigma_w2))/sqrt(sigma_e2*sigma_w2);

Init_inv_sigma_r2 = log(sigma_r2);
Init_inv_sigma_e2 = log(sigma_e2);
Init_inv_sigma_w2 = log(sigma_w2);

Init_pi = rho_eo;
Init_inv_pi = atanh(Init_pi/MaxValueCorr)/smoothConstantCorr; 

initialF1 = [Init_inv_sigma_r2; Init_inv_sigma_e2; Init_inv_sigma_w2; Init_inv_pi];

