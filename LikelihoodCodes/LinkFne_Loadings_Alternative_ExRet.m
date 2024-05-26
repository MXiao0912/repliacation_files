function [TransformmedZ,JacobianZ] = LinkFne_Loadings_Alternative_ExRet(UntransformmedZ,ARcoeffs,SettingsTransfor)

v2struct(SettingsTransfor);

switch CaseTransformationMu 
    case 'NoTransformation'
        PositiveTrans = @(x) x;
        DerPositiveTrans = @(x) 1;                
    case 'Exponential'
        PositiveTrans = @(x) exp(x);
        DerPositiveTrans = @(x) exp(x);        
    case 'Square' 
        PositiveTrans = @(x) (x.^2);
        DerPositiveTrans = @(x) 2*x;        
    case 'Logistic'
        MaxZero = @(maxVal,x) maxVal*(1./(1+exp(-x)));
        ChosenCap = MaxMu;
        PositiveTrans = @(x) MaxZero(ChosenCap,x);
        DerMaxZero = @(maxVal,x) maxVal*MaxZero(1,x).*(1-MaxZero(1,x));
        DerPositiveTrans = @(x) DerMaxZero(ChosenCap,x);
end

mu_bar = PositiveTrans(UntransformmedZ(1)); 
g_bar = UntransformmedZ(2); 
pd_bar = g_bar - log(exp(mu_bar)-exp(g_bar));

gamma1 = ARcoeffs(1); 
fi1 = ARcoeffs(2); 

Dpd_bar_Dmu_bar = - exp(mu_bar)/(exp(mu_bar)-exp(g_bar));
Dpd_bar_Dg_bar = 1 + (exp(g_bar)/(exp(mu_bar)-exp(g_bar)));

rho_bar = exp(pd_bar)/(1+exp(pd_bar));
%Drho_bar_Dpd_bar = exp(pd_bar)/(1+exp(pd_bar))^2;
Drho_bar_Dpd_bar =rho_bar*(1-rho_bar);


TransformmedZ = [g_bar; pd_bar; 1/(1-rho_bar*gamma1); -1/(1-rho_bar*fi1)];

JacobianZ1 = diag([DerPositiveTrans(UntransformmedZ(1)),1]);
JacobianZ = [0 1; ...
             Dpd_bar_Dmu_bar Dpd_bar_Dg_bar;...
             (1/(1-rho_bar*gamma1)^2)*Drho_bar_Dpd_bar*gamma1*Dpd_bar_Dg_bar 0;...
             -(1/(1-rho_bar*fi1)^2)*Drho_bar_Dpd_bar*fi1*Dpd_bar_Dmu_bar 0];
JacobianZ = JacobianZ*JacobianZ1;


