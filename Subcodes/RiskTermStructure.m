function [RiskTermStructure]=RiskTermStructure(hhoriz,Vpar,yy,initialF1,initialAlfa1,SettingsTransfor,No_Repetitions)

% No_Repetitions = 500; 
TT = max(size(yy));

DecaingPattern = @(ARcoef,n) (1 - ARcoef.^[n-1:-1:1]')./(1 - ARcoef);

for tt=1:TT
    tt/TT
%     rng(1000) % -->> Fix seed for replication 
    for jjrep = 1:No_Repetitions
        [Forecast,Res] = Forecast_TVKF_pd_Pd_Model(hhoriz,Vpar,yy(1:tt,:)',initialF1,initialAlfa1,SettingsTransfor);
        Store_Forecast_Q(:,:,jjrep) = Forecast.QFor; 
    end
%     Forecast_Q = mean(Store_Forecast_Q,3);
    TrimHowMuch = 5; % Percent
    Forecast_Q = trimmean(Store_Forecast_Q,TrimHowMuch,3);
    
    clear Store_Forecast_Q
    
    pdbar = Res.Z(2,1,tt); 
    rhor = exp(pdbar)/(1+exp(pdbar));
    
    ArCoeff_g = Res.T(2,2);
    ArCoeff_mu = Res.T(3,3);
    
    b_mu = abs(Res.Z(2,3,tt)); 
    b_g  = Res.Z(2,2,tt); 
    
    for tau = 1:size(Forecast_Q,2)
        QQ = reshape(Forecast_Q(:,tau),3,3);
        Var_eps_r_t = ((rhor*b_g)^2)*QQ(2,2) + QQ(1,1)+((rhor*b_mu)^2)*QQ(3,3)... 
            -2*(rhor^2)*(b_g*b_mu)*QQ(3,2) -2*(rhor)*(b_mu)*QQ(3,1); 
        Cov_eps_rmu_t = -rhor*b_mu*QQ(3,3) + rhor*b_g*QQ(3,2) + QQ(1,3);
        
        Var_r_tplusj_expret_comp(tau,1) = QQ(3,3); 
        Var_r_tplusj_unexpret_comp(tau,1) = Var_eps_r_t; 
        Var_r_tplusj_cov_comp(tau,1) = Cov_eps_rmu_t;
        
        Variance_CumReturns_tau = (1/tau)*[sum((DecaingPattern(ArCoeff_mu,tau).^2).*Var_r_tplusj_expret_comp(1:tau-1,1),1) ...
            sum(Var_r_tplusj_unexpret_comp(1:tau,1),1) ...
            2*sum((DecaingPattern(ArCoeff_mu,tau)).*Var_r_tplusj_cov_comp(1:tau-1,1),1)];
        Variance_CumReturns(tau,:) = [sum(Variance_CumReturns_tau,2) Variance_CumReturns_tau];
    end
    
    Store_Variance_CumReturns(:,:,tt) = Variance_CumReturns;
    clear Variance_CumReturns Var_r_tplusj_expret_comp Var_r_tplusj_unexpret_comp Var_r_tplusj_cov_comp

    %% DIVIDENDS
    for tau = 1:size(Forecast_Q,2)
        QQ = reshape(Forecast_Q(:,tau),3,3);
        
        Var_Dd_tplusj_expd_comp(tau,1) = (((1-ArCoeff_g^tau)/(1-ArCoeff_g))^2)*QQ(2,2); 
        Var_Dd_tplusj_unexpd_comp(tau,1) = QQ(1,1); 

        Variance_CumDividends_tau = (1/tau)*[sum((DecaingPattern(ArCoeff_mu,tau).^2).*Var_Dd_tplusj_expd_comp(1:tau-1,1),1) ...
            sum(Var_Dd_tplusj_unexpd_comp(1:tau,1),1) ];
        Variance_CumDividends(tau,:) = [sum(Variance_CumDividends_tau,2) Variance_CumDividends_tau];
    
    end
    
    Store_Variance_CumDividends(:,:,tt) = Variance_CumDividends;
    clear Variance_CumDividends Var_Dd_tplusj_expd_comp Var_Dd_tplusj_unexpd_comp

end

RiskTermStructure.CumReturns = Store_Variance_CumReturns;
RiskTermStructure.CumDividends = Store_Variance_CumDividends;

