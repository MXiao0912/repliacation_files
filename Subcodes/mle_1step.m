function [loglik,param_out,StdErr,SKF,PHIQ,PHI,SIGMA,Atau,Btau,Ctau,Dtau,Gtau,Htau] = mle_1step( param_in, data_yields, data_futures, info )

% [loglik,param_out,StdErr,SKF,PHIQ,PHI,SIGMA,Atau,Btau,Ctau,Dtau,Gtau,Htau] = mle_1step( param_in, data_yields, data_futures, info )
% This version of the code is for the model with stochastic seasonality


tau_yie     = info.tau_yie;     % Vector with the yield maturities
tau_fut     = info.tau_fut;     % Vector with the futures maturities
nmat_yie    = length(tau_yie);  % Number of maturities yields
nmat_fut    = length(tau_fut);  % Number of maturities futures
nobs        = info.nobs;        % Number of time-series observations
qqq         = info.qqq;         % Number of pairs of seasonal factors

Rescale_yields = 20;            % Rescaling parameter for yield curve block (20 works fine)

% Put all consecutive maturities 
tau_yie_all = (1:tau_yie(end))';

% Put all data together
data = [Rescale_yields*data_yields data_futures];

% To construct diffuse priors for the states:
small_eps = 10^-9;
large_eps = 10^4;

if info.estimate>0      % If do estimation:
      
    % If perform simulated annealing:
    if info.do_anneal > 0;
        fprintf('\n Optimizing using Simulated Annealing\n') 
        [ param_out, loglik ] = anneal( @log_likelihood, param_in, info.options_annealing );
        param_in = param_out;
    end
    
    % Check if do Nelder-Mead simplex method
    if info.do_simplex > 0
        fprintf('\n Optimizing using Nelder-Mead Simplex Method\n')
        [ param_out, loglik ] = fminsearch( @log_likelihood, param_in, info.options_fminsearch );
        param_in = param_out;
    end

    % Check if do quasi-newton optimization
    if info.do_qnewton >0
        if info.qn_type == 1    % use matlab "fminunc"
            fprintf('\n Optimizing using Quasi-Newton Method: fminunc \n')
            [ param_out, loglik ] = fminunc( @log_likelihood, param_in, info.options_fminunc );
        elseif info.qn_type==2   % use McGrattan's "uncmin"
            fprintf('\n Optimizing using Quasi-Newton Method: uncmin \n')
            [ param_out, loglik] = uncmin( param_in, @log_likelihood );
        else
            fprintf('\n Wrong input in info.qn_type: you must choose 1 or 2. \n')
            error(' ')
        end
    end
       
else  % Do not perform estimation (for graphs, computing standard errors, etc)
   loglik    = log_likelihood( param_in ); 
   param_out = param_in;
end

% Now compute filtered states and things like that. Everything is in the
% structure SKF
[ ~, ~, SKF, PHI_star ] = log_likelihood( param_out );

% Now compute smoothed estimates and add them to the structure SKF
SKF = FIS(PHI_star,SKF);

% COMPUTATION OF STANDARD ERRORS
if info.do_stderr >0
    fprintf('\n Computing asymptotic standard errors. This may take some time... ')
    switch info.sderr_type
        case 1   % use fdhess
            HHess = fdhess( @log_likelihood, param_out, 1, loglik );
        case 2 
            HHess = fdhess( @log_likelihood, param_out, 0, loglik );
        case 3   % use fdhess2
            HHess = fdhess2( @log_likelihood, param_out );
        case 4   % use outer product of the score
            jac = fdjacob( @contrib_likelihood, param_out, 1 ); 
            HHess = jac'*jac;
        otherwise
            fprintf('\n Wrong input in info.sderr_type \');
            error(' ')
    end
    StdErr = sqrt( diag( inv(HHess) ) );
    fprintf('Done. \n ')
else
    StdErr = NaN(length(param_out),1);
end


% NESTED FUNCTIONS: Likelihood function.
    function [loglik, contrib, SS, PHI_Trans ] = log_likelihood( param_in )

        % PREPARE PARAMETERS
        
        %   YIELD CURVE BLOCK PARAMETERS (except measurement errors)
        
        theta1 = exp(param_in(1));      % Nelson siegel parameter
        mu1Q   = param_in(2)/1000;      % Parameter of the risk neutral intercept
        
        % Parameters of the physical measure
        MU_delta  = param_in(3:5)/1000;             % Intercept physical measure
        
        PHI_delta = reshape(param_in(6:14),3,3);    % Matrix on lagged values physical measure
        
        GAMMA_delta = zeros(3,3);                   % Cholesky covariance matrix
        GAMMA_delta(1:3,1) = param_in(15:17)/1000;
        GAMMA_delta(2:3,2) = param_in(18:19)/1000;
        GAMMA_delta(3,3)   = param_in(20)/1000;
        
       
        %   FUTURES BLOCK
        theta2 = exp(param_in(21));
        omega  = exp(param_in(22));
        
        % Parameters of the physical measure
        MU_beta  = param_in(23:26)/1000;            % Intercept physical measure
        
        PHI_beta = reshape(param_in(27:42),4,4);    % Matrix on lagged values physical measure
        
        GAMMA_beta = zeros(4);                      % Cholesky futures factors (except seasonal)
        GAMMA_beta(1:4,1) = param_in(43:46)/1000;
        GAMMA_beta(2:4,2) = param_in(47:49)/1000;
        GAMMA_beta(3:4,3) = param_in(50:51)/1000;
        GAMMA_beta(4,4)   = param_in(52)/1000;

        % PARAMETERS OF THE PHYSICAL MEASURE
        MU = [MU_delta; MU_beta; 0; 0];     % The last two 0 are associated with the random walk
        
        PHI = zeros(9);
        PHI(1:3,1:3) = PHI_delta;
        PHI(4:7,4:7) = PHI_beta;
        PHI(8,8) = 1;
        PHI(9,9) = 1;
        
        % PARAMETERS OF THE RISK NEUTRAL MEASURE
        MUQ_yields  = [mu1Q; 0; 0];
        MUQ_futures = [0;0;0;0;0;0];
        
        MUQ = [MUQ_yields;MUQ_futures];
        
        PHIQ = zeros(9);
        
        PHIQ(1,1) = 1;
        PHIQ(2,2) = exp(-theta1);
        PHIQ(2,3) = theta1*exp(-theta1);
        PHIQ(3,3) = exp(-theta1);
        
        PHIQ(4,1) = 1;
        PHIQ(4,2) = (1-exp(-theta1))/theta1;
        PHIQ(4,3) = (1-exp(-theta1))/theta1-exp(-theta1);
        PHIQ(4,4) = 1;
        PHIQ(4,5) = 1;
        PHIQ(4,6) = (1-exp(-theta2))/theta2;
        PHIQ(4,7) = (1-exp(-theta2))/theta2-exp(-theta2);
        
        PHIQ(5,5) = 1;
        
        PHIQ(6,6) = exp(-theta2);
        PHIQ(6,7) = theta2*exp(-theta2);
        PHIQ(7,7) = exp(-theta2);
        
        PHIQ(8,8) = exp(-omega);
        PHIQ(9,9) = exp(-omega);
        
        % Variance of the seasonality components in the transition equation
        var_seas1 = (param_in(53)^2)/100;  % Constrained to be > 0.
        var_seas2 = (param_in(54)^2)/100;  % Constrained to be > 0.
        
        % Big OMEGA matrix
        SIGMA = zeros(9);
        SIGMA(1:3,1:3) = GAMMA_delta*GAMMA_delta';
        SIGMA(4:7,4:7) = GAMMA_beta*GAMMA_beta';
        SIGMA(8,8)     = var_seas1;
        SIGMA(9,9)     = var_seas2;
        
        % Parameters of the shortest yield
        rho0 = 0;
        rho1 = [1;(1-exp(-theta1))/theta1;(1-exp(-theta1))/theta1-exp(-theta1);0;0;0;0;0;0];
        
        % COMPUTE RECURSION FOR BONDS
        Btau = zeros(9,tau_yie(end));
        Atau = zeros(1,tau_yie(end));
        
        Atau(1)   = -rho0;
        Btau(:,1) = -rho1;
        
        for itau=2:tau_yie(end)
            Atau(itau)   = Atau(itau-1)-rho0 + MUQ'*Btau(:,itau-1) + 0.5*Btau(:,itau-1)'*SIGMA*Btau(:,itau-1);
            Btau(:,itau) = PHIQ'*Btau(:,itau-1) - rho1;
        end
        
        % Constant and loadings of yield on factors
        atau = -Atau./tau_yie_all';
        btau = bsxfun(@rdivide,-Btau,tau_yie_all');
        
        % Parameters of the log-spot price
        gamma0 = 0;
        gamma1 = cell(1,12);
        gamma1(:) = {[0;0;0;1;0;0;0;0;0]};   % Initialize to zeros and a 1 in the location of the deseasonalized spot price
        for im=1:12;
            gamma1{im}(8) = cos((pi/6)*im);
            gamma1{im}(9) = sin((pi/6)*im);
        end

        % COMPUTE RECURSION FOR FUTURES
       
        % Initialize Gtau^m and Htau^m, one cell per month and set to zero
        Gtau = cell(1,12);
        Htau = cell(1,12);
        
        % Initialize each cell to zeros
        Gtau(:) = {zeros(1,tau_fut(end))};
        Htau(:) = {zeros(9,tau_fut(end))};

        % Iterate over futures maturities
        for itau = 1:tau_fut(end)
    
            for imonth=1:12;
            
                imonthp1 = imonth+1 - 12*floor(imonth/12);
                
                if itau==1;     % First maturity is different
                    
                    H0 = gamma1{imonthp1};
                    
                    Htau{imonth}(:,1) = PHIQ'*H0-rho1;
                    
                    Gtau{imonth}(1)   = gamma0-rho0 + MUQ'*H0+0.5*H0'*SIGMA*H0;
                    
                else
                    
                    Htau{imonth}(:,itau) = PHIQ'*Htau{imonthp1}(:,itau-1)-rho1;

                    Gtau{imonth}(itau)   = Gtau{imonthp1}(itau-1) - rho0 + MUQ'*Htau{imonthp1}(:,itau-1) + 0.5*Htau{imonthp1}(:,itau-1)'*SIGMA*Htau{imonthp1}(:,itau-1);

                end
                
            end
        end
            
        % Construct Ctau^m and Dtau^m (parameters log-futures equation)
        Ctau = cell(1,12);
        Dtau = cell(1,12);
        DDtau = cell(1,12);
        for imonth = 1:12;
            Ctau{imonth} = Gtau{imonth} - Atau(1:tau_fut(end));
            Dtau{imonth} = Htau{imonth} - Btau(:,1:tau_fut(end));
            
            % The following is a loading vector that includes the constant
            % as a state variable
            DDtau{imonth} = [ Ctau{imonth} ; Dtau{imonth}];
        end
        
        % Here construct covariance matrix of measurement errors
        meas_errors_yields  = (param_in(55:61).^2)/10;
        meas_errors_futures = (param_in(62:85)/100).^2;
        
        SIGMA_Meas = diag([meas_errors_yields; meas_errors_futures]);

        % NOW WRITE TRANSITION EQUATION AS A VAR(1) ADDING THE CONSTANT AS 
        % ANOTHER STATE VARIABLE. First element of the state equation is the constant.
        
        PHI_Trans   = [ 1     zeros(1,9) ; 
                        MU    PHI       ];
                    
        SIGMA_Trans = zeros(10);
        SIGMA_Trans(2:10,2:10) = SIGMA;
                    
        % New yield curve loadings after expanding the state vector with the constant
        %       Yield curve for the relevant maturities in tau_yie vector): 
        bbtau = Rescale_yields*[atau(tau_yie) ; btau(:,tau_yie)];   %Rescale_yields is the rescaling parameter for yields. Here we rescale the loadings and constant
        
        % CONSTRUCT LOADING MATRIX OF STACKED MEASUREMENT EQUATION:
        %   Yields come first, futures come second
        
        % Loading matrix is time varying
        % columns are for [ constant, delta1, delta2, delta3, beta0, beta1, beta2, beta3, xi1t, xi2t ]
        LoadMatrix = zeros(nmat_yie + nmat_fut, 8+min(2*qqq,11), nobs );
        
        for it=1:nobs
            find_month=it-12*floor((it-1)/12);
            
            %   Loadings yields block
            LoadMatrix(1:nmat_yie,:,it) = bbtau';   % Loadings for yield curve equation
            
            %   Loadings futures block
            LoadMatrix(nmat_yie+1:nmat_yie+nmat_fut,:,it) = DDtau{find_month}';
        end
        
        % Diffuse prior for the states
        x_0   = [1 0 0 0 0 0 0 0 0 0]';  
        Sig_0 = diag([small_eps; large_eps*ones(9,1)]); 
        
        % Evaluate (negative) of log-likelihood function using Kalman filter
        % with time-varying matrices
        SS = KfilterMD(data',LoadMatrix,SIGMA_Meas,PHI_Trans,SIGMA_Trans,x_0,Sig_0);
        
        % Compute negative of log-likelihood (for minimization) and add the
        % potential penalty
        loglik = -SS.loglik;
        contrib = SS.contrib_loglik;

    end

    function contrib = contrib_likelihood( theta_in )
        % This function returns a vector with the contribution of each
        % observation to the log-likelihood. This is useful to compute
        % standard errors using the outer product of the gradient
        [~, contrib] = log_likelihood( theta_in );
    end % end contrib-likelihood nested function


    function S = KfilterMD(Y,Z,R,T,Q,A_0,P_0)
    %               --------------------------------
    % Computes Kalman filter for systems with time-varying matrices and
    % missing observations
    %
    % The model is   A[t+1] = T[t] * A[t] +  u[t]   => State equation
    %                Y[t]   = Z[t] * X[t] +  eps[t] => Observ. equation
    %
    % where u[t]~N(0,Q) and eps[t]~N(0,R).
    %
    % Data must be arranged so that each column is a time period and rows
    % refers to different maturities
    %
    %  INPUTS:
    %  Y: Data.
    %  Z: Loading matrix on unobserved states in observation equation
    %  R: Covariance matrix of observation errors
    %  T: Matrix on lagged values, state equation
    %  Q: Covariance matrix shocks in state equation
    %
    %  OUTPUT: Structure S with
    %  S.Am    :   Predicted state vector E[A[t]|t-1]
    %  S.AmU   :   Filtered state vector E[A[t]|t]
    %  S.Pm    :   Predicted covariance matrix of A[t]|t-1
    %  S.Pmu   :   Filtered covariance matrix of A[t]|t
    %  S.loglik:   Negative of value of log-likelihood (for minimization)
    %  S.contrib_loglik: Stores contribution to log-likelihood of each
    %                    observation
    %              -----------------------------------
    %
        % Grab dimensions
        nstates = size(Z,2);  

        % Initialize arrays
        S.Am  = nan(nstates, nobs);
        S.AmU = nan(nstates, nobs+1);
        S.Pm  = nan(nstates,nstates,nobs);
        S.PmU = nan(nstates,nstates,nobs+1);
        S.loglik = 0;
        S.contrib_loglik = zeros(nobs,1);    % To store contribution to log-likelihood function

        % Initial conditions state variables:
        Au = A_0;   % E[A[0]|t=0]
        Pu = P_0;   % Covariance matrix of A[0]|t=0

        S.AmU(:,1)   = Au;
        S.PmU(:,:,1) = Pu;

        IdP = eye(size(P_0));

        % Begin Kalman iteration
        for it=1:nobs

            % PREDICTION STEP: compute A = A[t]|t-1 and P = P[t]|t-1
            A = T*Au;
            P = T*Pu*T' + Q;
            %P = 0.5*(P+P');     % Make sure P is positive definite

            % HANDLING OF MISSING DATA
            ix  = ~isnan(Y(:,it));  % Index of variables with data
            y_t = Y(ix,it);         % Vector with data
            Z_t = Z(ix,:,it);       % Elements of Z with data
            R_t = R(ix,ix);         % Covariance of measurement errors with data

            % UPDATING STEP
            if isempty(y_t);        % If there is no data in this period
                Au = A;         
                Pu = P; 
            else
                PZ = P*Z_t';
                Ome = Z_t*PZ + R_t;
                PZF = PZ/Ome;
                V   = y_t - Z_t*A;
                Au  = A  + PZF*V;
                Pu  = (IdP- PZF*Z_t)*P*(IdP- PZF*Z_t)' + PZF*R_t*PZF';
                
                % Contribution of observation "it" to the log-likelihood
                % (ignoring the constant term)
                S.contrib_loglik(it) = -0.5*( log(det(Ome)) + V'*(Ome\V) );
            end
            
            S.Am(:,it)   = A;
            S.Pm(:,:,it) = P;

            % Au = A_t|t   & Pu = P_t|t
            S.AmU(:,it+1)    = Au;
            S.PmU(:,:,it+1)  = Pu;
            
        end % it
        S.loglik = sum( S.contrib_loglik );   % Log-likelihood function without the constant
        
        % The following is used to intialize the Kalman smoother
        if isempty(y_t)
            S.KZ = zeros(nstates);
        else
            S.KZ = PZF*Z_t;
        end
        S.T = T;

    end

    %______________________________________________________________________
    function S = FIS(T,S)
    %______________________________________________________________________
    % Fixed intervall smoother (see Harvey, 1989, p. 154)
    % FIS returns the smoothed state vector AmT and its covar matrix PmT             
    % Use this in conjnuction with function SKF
    %______________________________________________________________________
    % INPUT  
    %        Y         Data                                 (nobs x n)  
    %        S Estimates from Kalman filter SKF                                                          
    %          S.Am   : Estimates     a_t|t-1                  (nobs x m) 
    %          S.Pm   : P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
    %          S.AmU  : Estimates     a_t|t                    (nobs x m) 
    %          S.PmU  : P_t|t   = Cov(a_t|t)               (nobs x m x m)       
    % OUTPUT 
    %        S Smoothed estimates added to above
    %          S.AmT  : Estimates     a_t|T                    (nobs x m) 
    %          S.PmT :  P_t|T   = Cov(a_t|T)               (nobs x m x m)
    %          S.PmT_1 : Cov(a_ta_t-1|T)
    %        where m is the dim of state vector and t = 1 ...T is time

      nstates         = size(S.Am,1);
      S.AmT           = zeros(nstates,nobs+1);
      S.PmT           = zeros(nstates,nstates,nobs+1);
      S.AmT(:,nobs+1)   = squeeze(S.AmU(:,nobs+1))  ;
      S.PmT(:,:,nobs+1) = squeeze(S.PmU(:,:,nobs+1));
      S.PmT_1(:,:,nobs) = (eye(nstates)-S.KZ) *T*squeeze(S.PmU(:,:,nobs));

      J_2 = squeeze(S.PmU(:,:,nobs)) * T' * pinv(squeeze(S.Pm(:,:,nobs)));

      for t = nobs:-1:1 
          PmU = squeeze(S.PmU(:,:,t));
          Pm1 = squeeze(S.Pm(:,:,t));
          P_T = squeeze(S.PmT(:,:,t+1));
          P_T1 = squeeze(S.PmT_1(:,:,t));

          J_1 = J_2;

          S.AmT(:,t)   = S.AmU(:,t) + J_1 * (S.AmT(:,t+1) - T * S.AmU(:,t)) ; 
          S.PmT(:,:,t) = PmU        +  J_1 * (P_T - Pm1) * J_1'; 

         if t>1
              J_2 = squeeze(S.PmU(:,:,t-1)) * T' * pinv(squeeze(S.Pm(:,:,t-1)));
              S.PmT_1(:,:,t-1) = PmU*J_2'+J_1*(P_T1-T*PmU)*J_2';
          end

      end

    end % end FIS

end