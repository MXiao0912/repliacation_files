function [Decompositions, Expectations] = Decomposition_TVKF_pd_Pd_Model(Res,maxHoriz)

TT = size(Res.RT_alfa_t,2); 
Pd_CondVarDecomp = NaN*zeros(TT,4);

for tt = 1:TT
        %% TERM STRUCTURE OF EXPECTATIONS 
        pdbar = Res.Z(2,1,tt); 
        rhor = exp(pdbar)/(1+exp(pdbar));
        
        gtilda_t = Res.RT_alfa_t(2,tt);
        mutilda_t = Res.RT_alfa_t(3,tt); 
        
        gbar_t = Res.Z(1,1,tt);
        pdbar_t = Res.Z(2,1,tt); 
        mubar_t = log(exp(gbar_t-pdbar_t) + exp(gbar_t));
        
        ArCoeff_g = Res.T(2,2);
        ArCoeff_mu = Res.T(3,3);
        
        for jj = 1:maxHoriz
%             Annualized_Exp_g_t(jj,1) = (1/jj)*((1-(rhor)^jj)/(1-rhor))*gbar_t ...
%                 +(1/jj)*((1-(rhor*ArCoeff_g)^jj)/(1-rhor*ArCoeff_g))*gtilda_t; 
%             Annualized_Exp_mu_t(jj,1) = (1/jj)*((1-(rhor)^jj)/(1-rhor))*mubar_t ...
%                 +(1/jj)*((1-(rhor*ArCoeff_mu)^jj)/(1-rhor*ArCoeff_mu))*mutilda_t;             

            Annualized_Exp_g_t(jj,1) = gbar_t ...
                +(1/jj)*((1-(ArCoeff_g^jj))/(1-ArCoeff_g))*gtilda_t; 
            Annualized_Exp_mu_t(jj,1) = mubar_t ...
                +(1/jj)*((1-(ArCoeff_mu^jj))/(1-ArCoeff_mu))*mutilda_t;             

%             Annualized_Exp_g_t(jj,1) = (1/jj)*((1-(rhor*ArCoeff_g)^jj)/(1-rhor*ArCoeff_g))*gtilda_t; 
%             Annualized_Exp_mu_t(jj,1) = (1/jj)*((1-(rhor*ArCoeff_mu)^jj)/(1-rhor*ArCoeff_mu))*mutilda_t;             
        end
        
        Annualized_Exp_g(tt,:) = Annualized_Exp_g_t'; 
        Annualized_Exp_mu(tt,:) = Annualized_Exp_mu_t'; 
        
        b_mu = abs(Res.Z(2,3,tt)); 
        b_g  = Res.Z(2,2,tt); 
        
        %% PD COND. VARIANCE DECOMPOSITION 
        Pd_CondVarDecomp_Gcomp = (b_g^2)*Res.Q(2,2,tt); 
        Pd_CondVarDecomp_Rcomp = (b_mu^2)*Res.Q(3,3,tt); 
        Pd_CondVarDecomp_Covcomp = -2*(b_g*b_mu)*Res.Q(3,2,tt); 
        Pd_CondVarDecomp_all = Pd_CondVarDecomp_Gcomp + Pd_CondVarDecomp_Rcomp + Pd_CondVarDecomp_Covcomp;
        Pd_CondVarDecomp(tt,:) = [Pd_CondVarDecomp_all Pd_CondVarDecomp_Gcomp...
            Pd_CondVarDecomp_Rcomp Pd_CondVarDecomp_Covcomp];

        %% PD COND. COVARIANCE DECOMPOSITION With RETURNS 
        Pd_CondCovDecomp_Gcomp = rhor*(b_g^2)*Res.Q(2,2,tt); 
        Pd_CondCovDecomp_Rcomp = rhor*(b_mu^2)*Res.Q(3,3,tt); 
        Pd_CondCovDecomp_Covcomp = -2*rhor*(b_g*b_mu)*Res.Q(3,2,tt) - b_mu*Res.Q(3,1,tt); 
        Pd_CondCovDecomp_all = Pd_CondCovDecomp_Gcomp + Pd_CondCovDecomp_Rcomp + Pd_CondCovDecomp_Covcomp;
        Pd_CondCovDecomp(tt,:) = [Pd_CondCovDecomp_all Pd_CondCovDecomp_Gcomp...
            Pd_CondCovDecomp_Rcomp Pd_CondCovDecomp_Covcomp];
        
        %% PD COND. COVARIANCE DECOMPOSITION With RETURNS 
        Pd_CondCovDecompWithD(tt,:) = -(b_mu)*Res.Q(3,1,tt); 
        
        %% RET COND. VARIANCE DECOMPOSITION 
        Ret_CondVarDecomp_Dcomp = ((rhor*b_g)^2)*Res.Q(2,2,tt) + Res.Q(1,1,tt); 
        Ret_CondVarDecomp_Rcomp = ((rhor*b_mu)^2)*Res.Q(3,3,tt); 
        Ret_CondVarDecomp_Covcomp = -2*(rhor^2)*(b_g*b_mu)*Res.Q(3,2,tt) -2*(rhor)*(b_mu)*Res.Q(3,1,tt); 
        Ret_CondVarDecomp_all = Ret_CondVarDecomp_Dcomp + Ret_CondVarDecomp_Rcomp + Ret_CondVarDecomp_Covcomp;
        Ret_CondVarDecomp(tt,:) = [Ret_CondVarDecomp_all Ret_CondVarDecomp_Dcomp...
            Ret_CondVarDecomp_Rcomp Ret_CondVarDecomp_Covcomp];

        %% RET COND. COVARIANCE DECOMPOSITION 
        Ret_CondCovDecomp_Rcomp = -((rhor*b_mu))*Res.Q(3,3,tt); 
        Ret_CondCovDecomp_CovMugcomp = rhor*b_g*Res.Q(3,2,tt); 
        Ret_CondCovDecomp_CovMudcomp = Res.Q(3,1,tt); 
        Ret_CondCovDecomp_all = Ret_CondCovDecomp_Rcomp + Ret_CondCovDecomp_CovMugcomp + Ret_CondCovDecomp_CovMudcomp;
        Ret_CondCovDecomp(tt,:) = [Ret_CondCovDecomp_all Ret_CondCovDecomp_Rcomp...
            Ret_CondCovDecomp_CovMugcomp Ret_CondCovDecomp_CovMudcomp];
        
        %% RET COND. CORRELATION DECOMPOSITION 
        Ret_CondCorrDecomp(tt,:) = Ret_CondCovDecomp(tt,:)./sqrt(Ret_CondVarDecomp_all*Res.Q(3,3,tt)); 
        Pd_CondCorrDecomp(tt,:) = Pd_CondCovDecomp(tt,:)./sqrt(Pd_CondVarDecomp_all*Ret_CondVarDecomp_all); 
        Pd_CondCorrDecompWithD(tt,:) = Pd_CondCovDecompWithD(tt,:)./sqrt(Pd_CondVarDecomp_all*Res.Q(1,1,tt));
end

Decompositions.Pd_CondVarDecomp = Pd_CondVarDecomp; 
Decompositions.Pd_CondCovDecomp = Pd_CondCovDecomp; 
Decompositions.Pd_CondCorrDecomp = Pd_CondCorrDecomp; 
Decompositions.Pd_CondCovDecompWithD = Pd_CondCovDecompWithD; 
Decompositions.Pd_CondCorrDecompWithD = Pd_CondCorrDecompWithD;

Decompositions.Ret_CondVarDecomp = Ret_CondVarDecomp; 
Decompositions.Ret_CondCovDecomp = Ret_CondCovDecomp; 
Decompositions.Ret_CondCorrDecomp = Ret_CondCorrDecomp; 

Expectations.Annualized_Exp_g = Annualized_Exp_g; 
Expectations.Annualized_Exp_mu = Annualized_Exp_mu; 
