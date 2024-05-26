
clear; clc; 

path(path,'BreakTestsSubcodes');

WhereAmI = cd; 

%% LOAD DATA 

PickDataFreq = 'Annual'; %'Monthly'; %'Annual'; %'Quarterly'
% load([WhereAmI '\Data\SaveData_' PickDataFreq '.mat']); 
load([WhereAmI '\Data\SaveDataReal_' PickDataFreq '.mat']); 

%% PICK DATA 
%  y = [d, pd]

SelectVariables = [1,3,2];
StartDate =find(Dates == 1873);% 1;% 
EndDate = find(Dates == 2018); %size(All_Data,1); %
y = All_Data(StartDate:EndDate,SelectVariables); 
y_Names = All_Names(1,SelectVariables)'


% Options for Bai and Perron

% The followings are options if p > 0 -------------------------------    

fixb=0;                   % set to 1 if use fixed initial values for beta
betaini=0;                % if fixb=1, load the initial value of beta
maxi=20;                  % maximum number of iterations for the nonlinear procedure to 
                          % obtain global minimizers
printd=1;                 % set to 1 if want the output from the iterations to be printed
eps=0.0001;               % criterion for the convergence

robust=1;                % set to 1 if want to allow for heterogeneity and autocorrelation 
                         % in the residuals, 0 otherwise. The method used is Andrews(1991) 
                         % automatic bandwidth with AR(1) approximation and the quadratic 
                         % kernel. Note: Do not set to 1 if lagged dependent variables are 
                         % included as regressors.
prewhit=1;               % set to 1 if want to apply AR(1) prewhitening prior to estimating 
                         % the long run covariance matrix.
hetdat=1;                % option for the construction of the F tests. Set to 1 if want to
                         % allow different moment matrices of the regressors across segments. 
                         % If hetdat=0, the same moment matrices are assumed for each segment 
                         % and estimated from the ful sample. It is recommended to set 
                         % hetdat=1 if p>0.
hetvar=1;                % option for the construction of the F tests.Set to 1 if want to allow 
                         % for the variance of the residuals to be different across segments. 
                         % If hetvar=0, the variance of the residuals is assumed constant 
                         % across segments and constructed from the ful sample. This option 
                         % is not available when robust=1.  
hetomega=1;              % used in the construction of the confidence intervals for the break 
                         % dates. If hetomega=0, the long run covariance matrix of zu is 
                         % assumed identical across segments (the variance of the errors u 
                         % if robust=0)
hetq=1;                  % used in the construction of the confidence intervals for the break 
                         % dates. If hetq=0, the moment matrix of the data is assumed identical 
                         % across segments.
doglobal=1;              % set to 1 if want to cal the procedure to obtain global minimizers
dotest=1;                % set to 1 if want to construct the supF, UDmax and WDmax tests 
                         % doglobal must be set to 1 to run this procedure.
dospflp1=1;              % set to 1 if want to construct the supF(l+1|l) tests where under
                         % the null the l breaks are obtained using global minimizers. 
                         % doglobal must be set to 1 to run this procedure.
doorder=1;               % set to 1 if want to cal the procedure that selects the number of
                         % breaks using information criteria. doglobal must be set to 1 to 
                         % run this procedure.
dosequa=1;               % set to 1 if want to estimate the breaks sequentialy and estimate 
                         % the number of breaks using supF(l+1|l) test   
dorepart=1;              % set to 1 if want to modify the break dates obtained from the 
                         % sequential method using the repartition method of Bai (1995),
                         % Estimating breaks one at a time. This is needed for the confidence 
                         % intervals obtained with estim below to be valid.
estimbic=1;              % set to 1 if want to estimate the model with the number of breaks 
                         % selected by BIC.
estimlwz=0;              % set to 1 if want to estimate the model with the number of breaks  
                         % selected by LWZ
estimseq=1;              % set to 1 if want to estimate the model with the number of breaks
                         % selected using the sequential procedure
estimrep=1;              % set to 1 if want to estimate the model with the breaks selected
                         % using the repartition method
estimfix=0;              % set to 1 if want to estimate the model with a prespecified number
                         % of breaks equal to fixn set below
fixn=1;



%% Run Tests

lag_Nyblom = 1; 
No_variables = size(y,2); 
Nyblom_Test_Individual = nan(No_variables,lag_Nyblom+2);
Nyblom_Test_All = nan(No_variables,1);

Bai_Perron_Test_Results  = nan(No_variables,4);

for VariablePick = 1:No_variables
    
Variable_i = y(:,VariablePick);

            
            yy = Variable_i(isfinite(Variable_i(:,1)),1);
            lags=0;
            bigt=size(yy,1);
            z=ones(bigt,1);
            x=[];

            
            q=1;                      % number of regressors z
            p=lags;                   % number of regressors x 
            m=3;                      % maximum number of structural changes allowed
            eps1=.15;                 % value of the trimming (in percentage) for the construction 
                                      % and critical values of the supF type tests (used in the 
                                      % supF test, the Dmax, the supF(l+1|l) and the sequential 
                                      % procedure). If these tests are used, h below should be set 
                                      % at int(eps1*bigt). But if the tests are not required, 
                                      % estimation can be done with an arbitrary h. There are five 
                                      % options: eps1= .05, .10, .15, .20, or .25. For each option, 
                                      % the maximal value of m above is: 10 for eps=.05, 8 for 
                                      % eps1=.10, 5 for eps1=.15, 3 for eps1=.20, and 2 for eps1=.25.

            h=round(eps1*bigt);       % minimal length of a segment (h >= q). Note: if robust=1, h 
                                      % should be set at a larger value. 


           [ResultsTest,EstimatedBreakDates,DateMinSSbreak] = BPbreaktest(bigt,yy,z,q,m,h,eps1,robust,prewhit,hetomega,hetq,doglobal,dotest,dospflp1,...
                doorder,dosequa,dorepart,estimbic,estimlwz,estimseq,estimrep,estimfix,fixb,x,p,eps,maxi,betaini,printd,hetdat,hetvar,fixn,Dates);

           Bai_Perron_Test_Results(VariablePick,:) = [VariablePick ResultsTest.ftest(1) ResultsTest.cvftest(1:2,1)'];
            
           Table_Bai_Perror_Test(:,:,VariablePick) = [ResultsTest.ftest ResultsTest.cvftest';...
                                                      ResultsTest.supfl(:,1) ResultsTest.cvsupfl(:,[2,3])';...
                                                      ResultsTest.UDmax ResultsTest.cvUDmax'];
           Table_Bai_Perror_Dates(:,:,VariablePick) = [DateMinSSbreak];                                          
            
           WhereCell = {'b3','c3','d3'};
           xlswrite('EstimatesModelFinal',Table_Bai_Perror_Test(:,:,VariablePick),'Break_Tests',WhereCell{VariablePick}); 

           WhereCell = {'n3','n7','n11'};
           xlswrite('EstimatesModelFinal',Table_Bai_Perror_Dates(:,:,VariablePick),'Break_Tests',WhereCell{VariablePick}); 
           
        % Nyblom Test
        
            yy = Variable_i(isfinite(Variable_i(:,1)),1);
            yy = yy;
            lags=lag_Nyblom;
            
            if lags==0
                datnames=['______DY';'Constant';];
                Y = [yy(1:end,1) ones(size(yy(1:end,1)))]; 
            elseif lags==1
                datnames=['______DY';'Constant';'__DY(-1)'];
                Y = [yy(2:end,1) ones(size(yy(2:end,1))) yy(2-1:end-1,1)];
            elseif lags==2
                datnames=['______DY';'Constant';'__DY(-1)';'__DY(-2)'];
                Y = [yy(3:end,1) ones(size(yy(3:end,1))) yy(3-1:end-1,1) yy(3-2:end-2,1)];
            end        
            [li,lj]=LC(Y,datnames);
            
            Nyblom_Test_Individual(VariablePick,:) = li';
            Nyblom_Test_All(VariablePick) = lj;

            
            if VariablePick ==3; 
            Table_Nyblom_Test = [Nyblom_Test_Individual Nyblom_Test_All]';                                          
            xlswrite('EstimatesModelFinal',Table_Nyblom_Test,'Break_Tests','b13'); 
            end
            
end

CV_lines = ones(size(Nyblom_Test_Individual(:,1)))*[.47 .35];
testResults = [Bai_Perron_Test_Results Nyblom_Test_Individual(:,1) CV_lines];

