function [TransformmedZ,JacobianZ] = LinkFne_Loadings(UntransformmedZ,ARcoeffs)

pd_bar = UntransformmedZ(1); 
g_bar = UntransformmedZ(2); 

gamma1 = ARcoeffs(1); 
fi1 = ARcoeffs(2); 

rho_bar = exp(pd_bar)/(1+exp(pd_bar));
%Drho_bar_Dpd_bar = exp(pd_bar)/(1+exp(pd_bar))^2;
Drho_bar_Dpd_bar =rho_bar*(1-rho_bar);


TransformmedZ = [g_bar; pd_bar; 1/(1-rho_bar*gamma1); -1/(1-rho_bar*fi1)];

JacobianZ = [0 1; ...
             1 0;...
             (1/(1-rho_bar*gamma1)^2)*Drho_bar_Dpd_bar*gamma1 0;...
             -(1/(1-rho_bar*fi1)^2)*Drho_bar_Dpd_bar*fi1 0];

