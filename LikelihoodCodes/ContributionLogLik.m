function [Contribution] = ContributionLogLik(Vpar,yy,initialF1,initialAlfa1,SettingsTransfor); 

[LL,Res] = TVKF_pd_Pd_Model(Vpar,yy,initialF1,initialAlfa1,SettingsTransfor); 
Contribution = Res.ContributionLogLik;
