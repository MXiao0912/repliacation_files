%% common kf with arbitrarily large init variance 1e9
% specify likelihood function under a particular Vpar
function [LL,Res] = base_kf_missing(Vpar,yy)

% get information on the problem
N           = size(yy,1);               % Total no. of variables
TT          = size(yy,2);               % Observations to be used
m           = 5;                        % dimension of the state space (1 unobserved components)

run Useful_Transformations

%% 
if size(Vpar,1) ==1 
    Vpar = Vpar'; 
end

% map vpar to structural variables (4 variables)
psi = UnoMenoUno(Vpar(1));
phi = UnoMenoUno(Vpar(2));
H = zeros(N);
HL = diag(Vpar(3:5));
HL(2,1) = Vpar(6);
HL(3,1) = Vpar(7);
HL(3,2) = Vpar(8);
R = [1 0 0;0 1 0;0 0 1;1-phi 0 0;0 0 0];
Q = R*HL*(HL')*(R');

% H = HL*(HL');
% Q = zeros(2);
% Q(1,1) = Vpar(6)^2;


% non-time-varying system matrices
% Z = [1, -psi; 1-phi 0];
% T = [1 0;1 0];
Z = [1 1 0 0 0;0 0 1 1 0];
T = [1 0 0 0 0;0 psi 0 0 0;0 0 phi 0 0;0 0 0 1+phi -phi;0 0 0 1 0];  
  
%% INITIALIZE THE KF

% starting values for state vectors and state variances
a_1 = zeros(m,1);
kap = 1e9;

A = [1 0 0;0 0 0;0 0 0;0 1 0;0 0 1];
Pinf = A*A';
R0 = [0 0;1 0;0 1;0 0;0 0];
Q0 = reshape((eye((2^2))-kron(T(2:3,2:3),T(2:3,2:3)))\vec(Q(2:3,2:3)),2,2);
Pstar = R0*Q0*R0';
P_1 = kap*Pinf+Pstar;


alfa_t = [a_1 NaN*zeros(m,TT)];
P_t    = [P_1(:) NaN*zeros(m^2,TT)];
Lik_t  = NaN*zeros(1,TT);
  
% yy = yy(:,2:end) - [psi 0;0 phi]*yy(:,1:(end-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORGANIZE OUTPUT

Res.Z = cell(TT);
Res.v_t   = cell(TT);
Res.invF_t   = cell(TT);
Res.L_t = NaN(m,m,TT);
Res.alfa_t   = NaN(m,TT+1);
Res.P_t      = NaN(m,m,TT+1);
        
Res.alfa_t(:,1) = a_1; 
Res.P_t(:,:,1) = P_1;

for tt=1:TT+1

    
%% Kalman Filter 

    P = reshape(P_t(:,tt),m,m);
    a_t = alfa_t(:,tt);
    
    if tt<=TT
        
        %% STEP 1: KF STEP
        yt = yy(:,tt);
        [Wt,nt]=SelectMatW(yt);
        Wt = sparse(Wt);
        yt(isnan(yt)) = 0;
        Res.Z{tt} = full(Wt)*Z;
        
        v = full(Wt)*(yt - Z*a_t);
        Res.v_t{tt} = v;
        
        F =  full(Wt)*(Z*P*Z' + H)* full(Wt)';

        invF = F\eye(size(F));
        Res.invF_t{tt} = invF;

        G = sparse(P*Z'*full(Wt)'*invF);
        
        %% STEP 2: UPDATE Likelihood
        Lik_t(1,tt) = -.5*(N*log(2*pi) + log(det(F)) + v'*invF*v);
            
        %% UPDATE KF
        K = T*full(G); L = (T-K*full(Wt)*Z);  
        Res.L_t(:,:,tt) = L;

        alfa_t(:,tt+1)          = T*a_t + K*v;
        Res.alfa_t(:,tt+1)      = alfa_t(:,tt+1);

        P_temp                  = T*P*L'+Q;
        Res.P_t(:,:,tt+1)       = P_temp;
        P_t(:,tt+1)             = P_temp(:);

        
    end
end

Res.ContributionLogLik = Lik_t;
cutFirstObservations = 0; 
LL = sum(Lik_t(1,cutFirstObservations+1:end));


if isreal(LL)==0 
    LL = -10^10; 
end

if isfinite(LL)==0 
    LL = -10^10; 
end
