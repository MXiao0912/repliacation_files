function [LL,Res] = my_TVKF_pd_Pd_Model(Vpar,yy,initialF1,initialAlfa1,SettingsTransfor)

% global CNm

% get information on the problem
N           = size(yy,1);               % Total no. of variables
TT          = size(yy,2);               % Observations to be used

m           = 5;                        % dimension of the state space (1 unobserved components)
k           = 4;                        % No of score driven TV parameters (2 SV) 

%%
run Useful_Transformations

%% 
% f = [pd, g, ln(sigma_d), ln(sigma_g), ln(sigma_mu),atanh(pi_dmu),atanh(pi_gmu)]
%

%%
% the time varying parameters are collected in the vector f_t that follows
% the following law of motion
% f_t = A_coeff f_(t-1) + B_coeff s_t where s_t is scaled score
% B_coeff has a block diagonal structure and A_coeff is Identity (random walks)
% parameters, enter f_t  in the following order: H, Q
 
%% 
if size(Vpar,1) ==1 
    Vpar = Vpar'; 
end

%% SET UP THE DYNAMICS OF THE SCORE DRIVEN PARAMTERS 

% a) Parameters of the score driven process 
omega = Vpar(1:4); 
%TransfParam(1:5) = Vpar(1:5); 

A_coeff = diag(MaxZero(1,Vpar(5:8)));  
%TransfParam(6:10) = MaxZero(1,Vpar(6:10)); 

B_coeff     = Vpar(8:11);
%TransfParam(11:12) = .015+MaxZero(.3,Vpar(11:12)); 
%TransfParam(13:17) = MaxZero(.1,Vpar(13:17)); 
kap_hes     = MaxZero(.5,Vpar(9));
%TransfParam(18) = MaxZero(.5,Vpar(18)); 
% kap_hes     = .015;
% kap_hes     = .015;

% b) Paramters of the SS 
%TransfParam(19) = PositiveTrans(Vpar(19)); 
ARcoeffs_Trans = MaxZero(.99,Vpar([10,11]));
%TransfParam([20,21]) = MaxZero(.99,Vpar([20,21]));
% ARcoeffs_Trans = UnoMenoUno(Vpar([20,21]));

% c) Initialize the score driven process
f_1 = initialF1;
f_t = [f_1 NaN*zeros(k,TT)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now setup the selection matrices to be used below
% Given a general matrix Mt in the system, considering that some coefficients are constant
% you can write the vec(Mt) in general form as:
% vec(Mt) = Sm + SM1*(SM2*f_t) where
% 1-	Sm collects the time constant coefficients
% 2-	SM1 and SM2 are two selection matrices
% 3-	SM2 selects within f_t the elements that belong to Mt
% 4-	SM1 selects within Mt the elements that are time varying and places them in the appropriate position in vec(Mt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. Selection for Z (fixed Z =1)
%Loading = 1;
Z = [1 1 0 0 0; 0 0 1 1 0];

H = diag([0,0]);

%% 3. T loadings of the transitions equation (fixed T=1)
psi = ARcoeffs_Trans(1);
phi = ARcoeffs_Trans(2);
T = [1 zeros(1,4); zeros(2,1) diag([psi, phi]) zeros(2,2);...
    zeros(2,3) [1+phi -phi; 1 0]];

%% 4. Selection for Q
SS = [eye(2) zeros(2,1);zeros(1,2) 1;...
    1-phi 0 0;zeros(1,3)];

%% INITIALIZE THE RECURSIONS FOR THE SCORE DRIVEN MODEL 

Init_Q = (eye(4)-A_coeff)\omega;
f_t(:,1) = Init_Q;

[Q,Qdot] = my_LinkFne_Q(f_t(:,1),SettingsTransfor);
Q = SS*Q*SS';
kronSS = kron(SS,SS); 
Qdot = kronSS*Qdot;

% Nm = (.5*(eye(m^2)+CommMatrixK(m)));

  
  
%% INITIALIZE THE KF

% starting values for state vectors and state variances
a_1 = zeros(m,1);
if isempty(initialAlfa1)==0; a_1 =initialAlfa1; end
P_1 = reshape((eye((m)^2)-kron(T,T))\vec(Q),m,m);


alfa_t = [a_1 NaN*zeros(m,TT)];
P_t    = [P_1(:) NaN*zeros(m^2,TT)];
Lik_t  = NaN*zeros(1,TT);
  
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORGANIZE OUTPUT

Res.Z        = NaN(2,5,TT+1);
Res.Q        = NaN(3,3,TT+1);
Res.H        = zeros(2,2);
Res.T        = T;
% Res.FI       = NaN(p_lag,TT+1);
% Res.FI0       = NaN(1,TT+1);

Res.RT_alfa_t= NaN(m,TT);
Res.RT_P_t   = NaN(m,m,TT);
% Res.Lik_t    = NaN(1,TT);
% Res.f_t      = NaN(k,TT);
Res.alfa_t   = NaN(m,TT+1);
Res.P_t      = NaN(m,m,TT+1);
% Res.F_t = zeros(N,N,TT);
        
% additional output for Kalman smoother
% Res.v_t         = NaN(N,TT);
% Res.Pt          = NaN(m,m,TT+1);

Res.alfa_t(:,1) = a_1; 
Res.P_t(:,:,1) = P_1;

for tt=1:TT+1


%% STEP 0: Fill the matrices of the Meas. Eq. of the State Space before filtering
% 1. Matrix Z    
Zdot = zeros(10,k);
% 2. Matrix H: defined outside
    
%% Kalman Filter 

    P = reshape(P_t(:,tt),m,m);
    a_t = alfa_t(:,tt);
    
    if tt<=TT
        
        %% STEP 1: KF STEP
        yt = yy(:,tt);
        [Wt,nt]=SelectMatW(yt);
        Wt = sparse(Wt);
        yt(isnan(yt)) = 0;
        v = Wt*(yt - Z*a_t);
%         Res.v_t(:,tt) = v;
        
        F =  Wt*(Z*P*Z' + H)* Wt';
%         Res.F_t(~isnan(yt),~isnan(yt),tt) = F;
%         invF = inv(F);
        invF = F\eye(size(F));
        
        G = sparse(P*Z'*Wt'*invF); %G = P*Z'\F;
        
        %% STEP 1b: Real Time Filtering
        a_tt = a_t + G*v;
%         RT_alfa_t(:,tt) = a_tt;
        Res.RT_alfa_t(:,tt)=a_tt;
%         FitY(:,tt) = Z*RT_alfa_t(:,tt);
%         Res.FitY(:,tt)=FitY(:,tt);
                
        RT_Ptemp = P - G*F*G';
        P_tt = RT_Ptemp;
%         RT_P_t(:,tt) = RT_Ptemp(:);
        Res.RT_P_t(:,:,tt)=RT_Ptemp;
        
        %% STEP 2: UPDATE Likelihood
        Lik_t(1,tt) = -.5*(log(2*pi) + log(det(F)) + v'*invF*v);
%         Res.Lik_t(1,tt)=Lik_t(1,tt);
        
        %% STEP 3: COMPUTE SCORE
        

        Vdot = -(kron(a_t',eye(N))*Zdot);
        Nnt = sparse(.5*(eye(nt^2)+CommMatrixK(nt)));
        Fdot =  2*Nnt*(kron(Wt*Z*P,Wt))*Zdot + kron(Wt*Z,Wt*Z)*Qdot;

        invFkroninvF=kron(invF,invF);

        GradD = (Fdot'*(invFkroninvF)*(kron(v,v)-F(:)) - 2*Vdot'*invF*v);
        if tt==1
            InfMat_temp = Fdot'*(invFkroninvF)*Fdot + 2*Vdot'*invF*Vdot;
%             InfMat = (1-kap_hes) * (eye(size(InfMat_temp))*1) + kap_hes * InfMat_temp;
            InfMat = (1-kap_hes) * diag(diag(InfMat_temp)) + kap_hes * InfMat_temp;
        else
            InfMat_temp = Fdot'*(invFkroninvF)*Fdot + 2*Vdot'*invF*Vdot;
            InfMat =  (1-kap_hes) * InfMat + kap_hes * InfMat_temp;   %smoothing the Hessian of coeffs
        end
%     InfMat = Fdot'*(invFkroninvF)*Fdot + 2*Vdot'*invF*Vdot;
%         InfMat = diag(diag(InfMat));
        
        %% CHOOSE SCALING SCORE
        Score = InfMat\GradD;
%        Score = pinv(InfMat)*GradD;
%        Score = (InfMat^(-1))*GradD;
%        Score = (InfMat^(.5))\GradD;
         
%        Score = GradD;
        
        
        %% STEP 4: UPDATE PARAMETERS
        f_t(:,tt+1) = omega + A_coeff*f_t(:,tt) + B_coeff*Score;
%         f_t(:,tt+1) = f_t(:,tt) + B_coeff*Score;
%         Res.f_t(:,tt+1)=f_t(:,tt+1);
        
        
        %% STEP 5: UPDATE MATRICES T AND Q AND THEIR JACOBIANS GIVEN THE NEW PARAMETERS IN f_t(:,tt+1)
        % 1. Matrix T: defined outside
        % 2. Matrix Q
        [Q,Qdot] = LinkFne_Q(f_t(:,tt+1),SettingsTransfor);
        Res.Q(:,:,tt) = Q;
        Q = SS*Q*SS';
        Qdot = kronSS*Qdot;
        
        %% UPDATE KF
        K = T*G; L = (T-K*Z);  

        alfa_t(:,tt+1)          = T*a_t + K*v;
        Res.alfa_t(:,tt+1)      = alfa_t(:,tt+1);

        P_temp                  = T*P*L'+Q;
        Res.P_t(:,:,tt+1)       = P_temp;
        P_t(:,tt+1)             = P_temp(:);

        
    end
end

Res.TransfParam = TransfParam;
Res.omega = omega;
Res.A_coeff = A_coeff;  
Res.B_coeff = B_coeff;  
Res.kap_hes = kap_hes; 
Res.Sigma_Meas_pd = Sigma_Meas_pd;

Res.ContributionLogLik = Lik_t;
cutFirstObservations = 0; 
% LL = -sum(Lik_t(1,cutFirstObservations+1:end))/(TT);
LL = sum(Lik_t(1,cutFirstObservations+1:end));


if isreal(LL)==0 
    LL = -10^10; 
end

if isfinite(LL)==0 
    LL = -10^10; 
end



