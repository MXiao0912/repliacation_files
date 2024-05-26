function [TransformmedZ,JacobianZ] = LinkFne_Loadings_Alternative(UntransformmedZ,ARcoeffs)

mu_bar = UntransformmedZ(1); 
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

JacobianZ = [0 1; ...
             Dpd_bar_Dmu_bar Dpd_bar_Dg_bar;...
             (1/(1-rho_bar*gamma1)^2)*Drho_bar_Dpd_bar*gamma1*Dpd_bar_Dg_bar 0;...
             -(1/(1-rho_bar*fi1)^2)*Drho_bar_Dpd_bar*fi1*Dpd_bar_Dmu_bar 0];

