function [Qt,Qdot] = LinkFne_Q(TVparam,SettingsTransfor)

v2struct(SettingsTransfor);

% CaseTransformationVariance = 'Exponential'; %'Exponential'; %'Square'; 'Logistic';
% MinVol = 10^-5; 

switch CaseTransformationVariance
    case 'Absolute'
        PositiveTrans = @(x) MinVol + abs(x);
        DerPositiveTrans = @(x) sign(x);        
    case 'Exponential'
        PositiveTrans = @(x) MinVol + exp(x);
        DerPositiveTrans = @(x) exp(x);        
    case 'Square' 
        PositiveTrans = @(x) MinVol + (x.^2);
        DerPositiveTrans = @(x) 2*x;        
    case 'Logistic'
        MaxZero = @(maxVal,x) MinVol + maxVal*(1./(1+exp(-x)));
%         cap_vol = 5; 
        PositiveTrans = @(x) MaxZero(cap_vol,x);
        DerMaxZero = @(maxVal,x) maxVal*MaxZero(1,x).*(1-MaxZero(1,x));
        DerPositiveTrans = @(x) DerMaxZero(cap_vol,x);
end
        


% 

Ssigma=[1 0 0;...
        zeros(3,3);...
        0 1 0;...
        zeros(3,3);...
        0 0 1];
Sdelta =[zeros(3,2) eye(3) zeros(3,2)];   
Srho = [zeros(2,2); ...
        1 0;...
        zeros(2,2);...
        0 1;...
        1 0;...
        0 1;...
        0 0];
Sgamma = [zeros(2,5) eye(2)];
deltas = Sdelta*TVparam; 
D = diag(PositiveTrans(deltas)); 

gammas = Sgamma*TVparam;
[R,Pi,Psi31,Psi32] = MappingGammaToCorrelations(gammas,SettingsTransfor);

Qt = D*R*D;

II = eye(size(Qt,1));

dvecQ_dvecD = kron(D*R,II) + kron(II,D*R);
dvecD_dSigma = Ssigma; 
dSigma_ddelta = diag(DerPositiveTrans(deltas)); 
ddelta_df = Sdelta; 

dvecQ_dvecR = kron(D,D);
dvecR_drho = Srho; 
dgamma_df = Sgamma; 

drho_dpi = [1 0; -Pi(2,3)*Pi(1,3)/sqrt(1-Pi(1,3)^2) sqrt(1-Pi(1,3)^2)];
dpi_dgamma = diag([Psi31,Psi32]);

Qdot = dvecQ_dvecD*dvecD_dSigma*dSigma_ddelta*ddelta_df +...
        dvecQ_dvecR*dvecR_drho*drho_dpi*dpi_dgamma*dgamma_df;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_t,Pi_t,Psi31,Psi32] = MappingGammaToCorrelations(gammas,SettingsTransfor)

v2struct(SettingsTransfor);

% smoothConstantCorr = 10^-1;%10^-1;%1;%10^-2; 
% MaxValueCorr = .95; 

switch CaseTransformationCorr
    case 'Arctan'
        CorrTrans = @(x) MaxValueCorr*tanh(smoothConstantCorr*x);
        DerCorrTrans = @(x) MaxValueCorr*smoothConstantCorr*(1-x.^2);   
    case 'Hypersherical'
        CorrTrans = @(x) MaxValueCorr*cos(x);
        DerCorrTrans = @(x) -MaxValueCorr*sin(x);
end

Pi_vect = CorrTrans(gammas)';
Pi_t = eye(3) + [zeros(2,3); Pi_vect 0] +[zeros(2,3); Pi_vect 0]';

Psi31 =DerCorrTrans(Pi_t(3,1)); 
Psi32 =DerCorrTrans(Pi_t(3,2));

rho_31 = Pi_t(3,1); 
rho_32 = Pi_t(3,2)*sqrt(1-Pi_t(3,1)^2);

R_t = eye(3) + [zeros(2,3); [rho_31 rho_32] 0] + [zeros(2,3); [rho_31 rho_32] 0]';





