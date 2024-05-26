

v2struct(SettingsTransfor)

%% SET INITIAL VALUES TIME VARYING ELEMENTS 
% f = [pd, g, ln(sigma_d), ln(sigma_g), ln(sigma_mu),atanh(pi_dmu),atanh(pi_gmu)]

try QUANTI_InitObs
    NoInitObs = QUANTI_InitObs;
catch 
    NoInitObs = 40;
end
% NoInitObs = 40; 

Init_pd_bar = mean(All_Data(1:NoInitObs,1)); %% MEAN Returns 
% Init_pd_bar = mean(All_Data(1:NoInitObs,2)); %% MEAN P/D 
Init_g_bar = mean(All_Data(1:NoInitObs,3)); %% MEAN Div Growth 

%% Initialize ExpReturns 
y_init = (All_Data(2:NoInitObs,1));  
x_init = [ones(size(y_init)) (All_Data(1:NoInitObs-1,1:3))];  
% x_init = [ones(size(y_init)) (All_Data(1:NoInitObs-1,2:3))];  
ExpRet = x_init*(x_init\y_init);

y_init = ExpRet(2:end,1); 
x_init = [ones(size(y_init)) (ExpRet(1:end-1,1))];  
ExpRetShock = y_init - x_init*(x_init\y_init); 

sigma_mu = std(ExpRetShock); 

%% Initialize ExpDividendGrowth 
y_init = (All_Data(2:NoInitObs,3));  
x_init = [ones(size(y_init)) (All_Data(1:NoInitObs-1,1:3))];  
% x_init = [ones(size(y_init)) (All_Data(1:NoInitObs-1,2:3))];  
ExpDiv = x_init*(x_init\y_init);
UnexExpDivShock = y_init - ExpDiv; 
UnexExpDivShock = UnexExpDivShock(1:end-1,1);
sigma_d = std(UnexExpDivShock); 

y_init = ExpDiv(2:end,1); 
x_init = [ones(size(y_init)) (ExpDiv(1:end-1,1))];  
ExpDivShock = y_init - x_init*(x_init\y_init); 

sigma_g = std(ExpDivShock); 

corr_dmu = corr(ExpRetShock,UnexExpDivShock) 
corr_gmu = corr(ExpRetShock,ExpDivShock) 

switch CaseTransformationVariance
    case 'Absolute'
        Init_inv_sigma_d = (sigma_d);
        Init_inv_sigma_g = (sigma_g);
        Init_inv_sigma_mu = (sigma_mu);
    case 'Exponential'
        Init_inv_sigma_d = log(sigma_d);
        Init_inv_sigma_g = log(sigma_g);
        Init_inv_sigma_mu = log(sigma_mu);
    case 'Square' 
        Init_inv_sigma_d = sqrt(sigma_d);
        Init_inv_sigma_g = sqrt(sigma_g);
        Init_inv_sigma_mu = sqrt(sigma_mu);
    case 'Logistic'
%         cap_vol = 5; 
        Init_inv_sigma_d = InvMaxZero(cap_vol,sigma_d);
        Init_inv_sigma_g = InvMaxZero(cap_vol,sigma_g);
        Init_inv_sigma_mu = InvMaxZero(cap_vol,sigma_mu);
end



Init_pi_dmu = corr_dmu; 
Init_pi_gmu = corr_gmu/sqrt(1-corr_dmu^2);

switch CaseTransformationCorr 
    case 'Arctan'
        Init_inv_pi_dmu = atanh(Init_pi_dmu/MaxValueCorr)/smoothConstantCorr; 
        Init_inv_pi_gmu = atanh(Init_pi_gmu/MaxValueCorr)/smoothConstantCorr;
    case 'Hypersherical'
        Init_inv_pi_dmu = acos(Init_pi_dmu); 
        Init_inv_pi_gmu = acos(Init_pi_gmu); 
end
initialF1 = [Init_pd_bar; Init_g_bar; Init_inv_sigma_d; Init_inv_sigma_g; Init_inv_sigma_mu; Init_inv_pi_dmu; Init_inv_pi_gmu];

