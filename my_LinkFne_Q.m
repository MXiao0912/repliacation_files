function [Qt,Qdot] = my_LinkFne_Q(TVparam,SettingsTransfor)

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
% this S select the param from vec(D)
Sdelta =[eye(3) zeros(3,1)];   
% this S select from ft
Srho = [zeros(5,1); ...
        1 ;...
        0 ;...
        1 ;...
        0];
% select from vec(R)
Sgamma = [zeros(1,3) 1];
% select from ft
deltas = Sdelta*TVparam; 
D = diag(PositiveTrans(deltas)); 

gammas = Sgamma*TVparam;
[R,Pi,Psi31] = MappingGammaToCorrelations(gammas,SettingsTransfor);

Qt = D*R*D;

II = eye(size(Qt,1));

dvecQ_dvecD = kron(D*R,II) + kron(II,D*R);
dvecD_dSigma = Ssigma; 
dSigma_ddelta = diag(DerPositiveTrans(deltas)); 
ddelta_df = Sdelta; 

dvecQ_dvecR = kron(D,D);
dvecR_drho = Srho; 
dgamma_df = Sgamma; 

Qdot = dvecQ_dvecD*dvecD_dSigma*dSigma_ddelta*ddelta_df +...
        dvecQ_dvecR*dvecR_drho*Psi31*dgamma_df;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_t,Pi_t,Psi31] = MappingGammaToCorrelations(gammas,SettingsTransfor)

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

Pi_vect = CorrTrans(gammas);
Pi_t = eye(3) + [zeros(2,3); 0 Pi_vect 0] +[zeros(2,3);0 Pi_vect 0]';

Psi31 =DerCorrTrans(Pi_t(3,2)); 

rho_32 = Pi_t(3,2);

R_t = eye(3) + [zeros(2,3); 0 rho_32 0] + [zeros(2,3); 0 rho_32 0]';





