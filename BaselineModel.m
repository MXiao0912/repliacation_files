
clear; clc; 

FileName_ResultsSave = 'TVP_BAselineModel_Estimates_10Sept_Prova';
FileName_RiskTermStructSave = 'RiskTermStructure_Final_10Sept_Prova';

%% DECIDE IF ESTIMATE OR USE SAVED DATA 
EstimateModel = 'N'; 
clear EstimateRiskTermStructure
EstimateRiskTermStructure = 'N'; 

%%
CaseToEstimate = 'AllCoeffs'; %'AllCoeffs'; % 'MinimumCoeffs'; 'MediumCoeffs'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CODES FOR TVP MODEL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SettingsTransfor.CaseTransformationVariance = 'Exponential'; % 'Absolute';%'Exponential'; %'Square'; 'Logistic';
SettingsTransfor.MinVol = 10^-5;
SettingsTransfor.cap_vol = 5; 
SettingsTransfor.CaseTransformationCorr = 'Arctan'; %'Arctan'; % 'Hypersherical';
SettingsTransfor.smoothConstantCorr = 1*10^-1;%10^-1;%1;%10^-2; 
SettingsTransfor.MaxValueCorr = .95; 

%% 
path(path,'Subcodes'); 
path(path,'LikelihoodCodes'); 
path(path,'Sims Solver');

WhereAmI = cd; 


if EstimateModel == 'Y'
%% LOAD DATA 

PickDataFreq = 'Annual'; %'Monthly'; %'Annual'; %'Quarterly'
% load([WhereAmI '\Data\SaveData_' PickDataFreq '.mat']); 
load([WhereAmI '\Data\SaveDataReal_' PickDataFreq '.mat']); 

%% PICK DATA 
%  y = [d, pd]

SelectVariables = [3,2];
StartDate =find(Dates == 1873);% 1;% 
EndDate = find(Dates == 2018); %size(All_Data,1); %
y = All_Data(StartDate:EndDate,SelectVariables); 
y_Names = All_Names(1,SelectVariables)';


%% USEFUL TRANSFORMATIONS 
run Useful_Transformations

%% SET INITIAL VALUES TIME VARYING ELEMENTS 
% f = [pd, g, ln(sigma_d), ln(sigma_g), ln(sigma_mu),atanh(pi_dmu),atanh(pi_gmu)]
Initialize_TVP

%% INITIAL PARAMETERS 

Init_A_coeff = .9*ones(5,1); 
Init_A_coeff = InvMaxZero(1,Init_A_coeff); 

Init_Omega = (eye(5)-diag(MaxZero(1,Init_A_coeff)))*initialF1(3:end);

Init_B_coeff_ss = [0.05; 0.05;]; % 0.2/0.08; SQRT: .05; .01
Init_B_coeff_QQ = 0.015*ones(5,1); % 0.015; SQRT: .07
Init_B_coeff = [InvMaxZero(.3,[Init_B_coeff_ss]); InvMaxZero(.1,[Init_B_coeff_QQ])]; 

Init_kap_hes = InvMaxZero(.5,.015);
Init_Sigma_Meas_pd = InvPositiveTrans(.001);

Init_gamma1 = .4; %.3
Init_gamma1 = InvMaxZero(.99,Init_gamma1);  

Init_fi1 = .7; %.9
Init_fi1 = InvMaxZero(.99,Init_fi1);  

switch CaseToEstimate 
    case 'AllCoeffs'
        InitialParams = [Init_Omega; Init_A_coeff; Init_B_coeff; Init_kap_hes;...
                         Init_Sigma_Meas_pd; Init_gamma1; Init_fi1]; 
    case 'MinimumCoeffs'
        InitialParams = [Init_Omega; Init_A_coeff(1); Init_A_coeff(1); Init_B_coeff(1); Init_B_coeff(3);...
                         Init_B_coeff(6); Init_kap_hes;...
                         Init_Sigma_Meas_pd; Init_gamma1; Init_fi1]; 
    case 'MediumCoeffs'
        InitialParams = [Init_Omega; Init_A_coeff(1); Init_A_coeff(1); Init_B_coeff(1); Init_B_coeff(2); Init_B_coeff(3);...
                         Init_B_coeff(6); Init_kap_hes;...
                         Init_Sigma_Meas_pd; Init_gamma1; Init_fi1]; 
end

%% DEFINE LIKELIHOOD FUNCTION

switch CaseToEstimate 
    case 'AllCoeffs'
        LossToMinmize = @(vparam) -TVKF_pd_Pd_Model([vparam],y',initialF1,[],SettingsTransfor)/size(y,1);
    case 'MinimumCoeffs'
        LossToMinmize = @(vparam) -TVKF_pd_Pd_Model([vparam(1:5,1); ...
                        vparam(6,1)*ones(3,1);vparam(7,1)*ones(2,1);vparam(8,1)*ones(2,1);...
                        vparam(9,1)*ones(3,1);vparam(10,1)*ones(2,1);vparam(11:14,1)],y',initialF1,[],SettingsTransfor)/size(y,1);
    case 'MediumCoeffs'
        LossToMinmize = @(vparam) -TVKF_pd_Pd_Model([vparam(1:5,1); ...
                        vparam(6,1)*ones(3,1);vparam(7,1)*ones(2,1);vparam(8,1);vparam(9,1);...
                        vparam(10,1)*ones(3,1);vparam(11,1)*ones(2,1);vparam(12:15,1)],y',initialF1,[],SettingsTransfor)/size(y,1);
end

             
%% MAXIMIZE LIKELIHOOD

% optionsIVAN = optimset('Display', 'iter-detailed','LargeScale', 'off','MaxFunEvals',5000);
% [EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);
% [EstimParams] = fminsearch(LossToMinmize,EstimParams,optionsIVAN);
% [EstimParams] = fminunc(LossToMinmize,EstimParams,optionsIVAN);

% [EstimParams] = fminsearch(LossToMinmize,InitialParams,optionsIVAN);
% [EstimParams] = fminunc(LossToMinmize,EstimParams,optionsIVAN);

[fhat,EstimParams] = csminwel(LossToMinmize,InitialParams,eye(size(InitialParams,1))*.5,[] ,1e-14,100);

%% PLOT RESULTS
vparam = EstimParams; 

switch CaseToEstimate 
    case 'AllCoeffs'
        EstimParamsNew = [vparam];
    case 'MinimumCoeffs'
        EstimParamsNew = [vparam(1:5,1); ...
                        vparam(6,1)*ones(3,1);vparam(7,1)*ones(2,1);vparam(8,1)*ones(2,1);...
                        vparam(9,1)*ones(3,1);vparam(10,1)*ones(2,1);vparam(11:14,1)];
    case 'MediumCoeffs'
        EstimParamsNew = [vparam(1:5,1); ...
                        vparam(6,1)*ones(3,1);vparam(7,1)*ones(2,1);vparam(8,1); vparam(9,1);...
                        vparam(10,1)*ones(3,1);vparam(11,1)*ones(2,1);vparam(12:15,1)];
end

[LL,Res] = TVKF_pd_Pd_Model(EstimParamsNew,y',initialF1,[],SettingsTransfor); 

Collect(:,1) = squeeze(Res.Q(3,2,:))./sqrt(squeeze(Res.Q(2,2,:)).*squeeze(Res.Q(3,3,:)));
Collect(:,2) = squeeze(Res.Q(3,1,:))./sqrt(squeeze(Res.Q(1,1,:)).*squeeze(Res.Q(3,3,:)));
Collect(:,3) = squeeze(sqrt(Res.Q(1,1,:)));
Collect(:,4) = squeeze(sqrt(Res.Q(2,2,:)));
Collect(:,5) = squeeze(sqrt(Res.Q(3,3,:)));
Collect(:,6) = squeeze(Res.Z(1,1,:)); 
Collect(:,7) = squeeze(Res.Z(2,1,:)); 

% CALCULATE mu_bar 
g_bar = Collect(:,6);  
pd_bar = Collect(:,7);
mu_bar = log(exp(g_bar-pd_bar) + exp(g_bar));

paramTransform = Res.TransfParam'; open paramTransform

%% DECOMPOSE PD 

TT = size(Res.RT_alfa_t,2); 
for tt=1:TT; 
    % Trend 
    Fit_Pd(tt,1) = Res.Z(2,1,tt)*Res.RT_alfa_t(1,tt); 
    % Mu Component
    Fit_Pd(tt,3) = Res.Z(2,3,tt)*Res.RT_alfa_t(3,tt); 
    % G Component
    Fit_Pd(tt,4) = Res.Z(2,2,tt)*Res.RT_alfa_t(2,tt); 
end

Fit_Pd(:,2) = sum([Fit_Pd(:,[3:4])],2); 

%% CALCULATE SE ESTIMATES WITH SIMULATIONS (coded up only 'AllCoeffs' specification)
% Fix seed 
rng(1000);
LossToMinmize = @(vparam) TVKF_pd_Pd_Model([vparam],y',initialF1,[],SettingsTransfor);
ContrLossToMinmize = @(vparam) ContributionLogLik([vparam],y',initialF1,[],SettingsTransfor);
TT = size(y,1); 

HHessian = fdhess2(LossToMinmize, EstimParamsNew );
jac = fdjacob(ContrLossToMinmize, EstimParamsNew, 1 );
HHessianSand = jac'*jac; 

[TransParam,Jacob] = TransformParametersLik(EstimParamsNew);
TransParam = TransParam';

VARIANZA_Coeffs =(1/TT)*inv(HHessianSand);
% VARIANZARob_Coeffs =(1/TT)*inv(HHessian)*HHessianSand*inv(HHessian);

VARIANZA_TransfCoeffs =(1/TT)*diag(Jacob)*inv(HHessianSand)*diag(Jacob);
% VARIANZARob_TransfCoeffs =(1/TT)*diag(Jacob)*inv(HHessian)*HHessianSand*inv(HHessian)*diag(Jacob);

PickNoSimulations = 1000; 

Save_Collect.Corr32 = NaN*zeros(TT+1,PickNoSimulations);
Save_Collect.Corr31 = NaN*zeros(TT+1,PickNoSimulations);
Save_Collect.Std11 = NaN*zeros(TT+1,PickNoSimulations);
Save_Collect.Std22 = NaN*zeros(TT+1,PickNoSimulations);
Save_Collect.Std33 = NaN*zeros(TT+1,PickNoSimulations);
Save_Collect.gbar = NaN*zeros(TT+1,PickNoSimulations);
Save_Collect.mubar = NaN*zeros(TT+1,PickNoSimulations);
Save_Collect.pdbar = NaN*zeros(TT+1,PickNoSimulations);
Save_Collect.TransfParam = NaN*zeros(size(TransParam,1),PickNoSimulations);

Save_Collect.Corr32(:,1) = Collect(:,1);
Save_Collect.Corr31(:,1) = Collect(:,2);
Save_Collect.Std11(:,1) = Collect(:,3);
Save_Collect.Std22(:,1) = Collect(:,4);
Save_Collect.Std33(:,1) = Collect(:,5);
Save_Collect.gbar(:,1) = g_bar;
Save_Collect.mubar(:,1) = mu_bar;
Save_Collect.pdbar(:,1) = pd_bar;
Save_Collect.TransfParam(:,1) = TransParam;

for Simu = 2:PickNoSimulations
    
    Simu 
    CheckFinite = 0; 
    
    while CheckFinite==0
%     Paramter_i = mvnrnd(EstimParamsNew,VARIANZARob_Coeffs);
    Paramter_i = mvnrnd(EstimParamsNew,VARIANZA_Coeffs);
    
    [LL_i,Res_i] = TVKF_pd_Pd_Model(Paramter_i,y',initialF1,[],SettingsTransfor); 
    
        if LL_i~=-1.0000e+10
        Save_Collect.Corr32(:,Simu) = squeeze(Res_i.Q(3,2,:))./sqrt(squeeze(Res_i.Q(2,2,:)).*squeeze(Res_i.Q(3,3,:)));
        Save_Collect.Corr31(:,Simu) = squeeze(Res_i.Q(3,1,:))./sqrt(squeeze(Res_i.Q(1,1,:)).*squeeze(Res_i.Q(3,3,:)));
        Save_Collect.Std11(:,Simu) = squeeze(sqrt(Res_i.Q(1,1,:)));
        Save_Collect.Std22(:,Simu) = squeeze(sqrt(Res_i.Q(2,2,:)));
        Save_Collect.Std33(:,Simu) = squeeze(sqrt(Res_i.Q(3,3,:)));

        g_bar_i = squeeze(Res_i.Z(1,1,:)); 
        pd_bar_i = squeeze(Res_i.Z(2,1,:)); 
        mu_bar_i = log(exp(g_bar_i-pd_bar_i) + exp(g_bar_i));
        paramTransform_i = Res_i.TransfParam'; 

        Save_Collect.gbar(:,Simu) = g_bar_i;
        Save_Collect.mubar(:,Simu) = mu_bar_i;
        Save_Collect.pdbar(:,Simu) = pd_bar_i;
        Save_Collect.TransfParam(:,Simu) = paramTransform_i;
        
        if all(Save_Collect.Std11(1:end-1,Simu)<100) && all(Save_Collect.Std22(1:end-1,Simu)<100)...
                && all(Save_Collect.Std33(1:end-1,Simu)<100);  
            CheckFinite = 1;  
        end
        else 
            'Salta Draw'
        end
    
    end
end

%% SAVE ALL RESULTS 
save(FileName_ResultsSave)

elseif EstimateModel == 'N'  
%% LOAD RESULTS 
load(FileName_ResultsSave)

end

%% SAVE TABLE WITH ESTIMATES 

% A) TABLE MAIN ESTIMATES = [mubar; gbar; rhomu; rhog; sig11; sig22; sig33; Corr31; Corr32; SigME; LL]
TableMainEstimates(:,1) = [NaN; NaN; TransParam(21); TransParam(20); ...
    nanmean(Save_Collect.Std11(:,1)); nanmean(Save_Collect.Std22(:,1)); nanmean(Save_Collect.Std33(:,1));...
    nanmean(Save_Collect.Corr31(:,1)); nanmean(Save_Collect.Corr32(:,1)); TransParam(19); LL]; 

SE_TransfCoeffs = diag(VARIANZA_TransfCoeffs).^.5; 
TableMainEstimates(:,2) = [NaN; NaN; SE_TransfCoeffs(21); SE_TransfCoeffs(20); ...
    quantile(nanmean(Save_Collect.Std11),[.16]); quantile(nanmean(Save_Collect.Std22),[.16]); quantile(nanmean(Save_Collect.Std33),[.16]);...
    quantile(nanmean(Save_Collect.Corr31),[.16]); quantile(nanmean(Save_Collect.Corr32),[.16]); SE_TransfCoeffs(19); NaN]; 
TableMainEstimates(:,3) = [NaN; NaN; NaN; NaN; ...
    quantile(nanmean(Save_Collect.Std11),[.84]); quantile(nanmean(Save_Collect.Std22),[.84]); quantile(nanmean(Save_Collect.Std33),[.84]);...
    quantile(nanmean(Save_Collect.Corr31),[.84]); quantile(nanmean(Save_Collect.Corr32),[.84]); NaN; NaN]; 

% B) TABLE GAS COEFFICIENTS = [diag(A); diag(B); kappa_h]
TableGASEstimates(:,1) = [TransParam(6:18);]; 
TableGASEstimates(:,2) = [SE_TransfCoeffs(6:18);]; 

xlswrite('EstimatesModelFinal',TableMainEstimates,'TVP_Model_Prova','c4'); 
xlswrite('EstimatesModelFinal',TableGASEstimates,'TVP_Model_Prova','h4'); 

%% PLOTS 
WhereAmI = cd; 
WhereToSave = [WhereAmI '\Store_Figures\'];

SPESSORE = 1; 
FONTsizePick = 18;


%% PLOT LONG RUN STEADY STATES 
Quantili = [.16 .5 .84];
SelectTiming = [1:TT]; 
% SelectTiming = [2:TT+1]; 

ColorLine{1} = 'b';  
ColorLine{2} = 'r'; 
ColorLine{3} = 'g';  

ColorShade{1} = [.7 .9 1];  
ColorShade{2} = [1 .7 .8]; 
ColorShade{3} = [.8 1 .8];  

            figure
            
            plot1 = plot(DatesNum,NaN*zeros(TT,2),'linewidth',SPESSORE); 
            set(plot1(1),'color',ColorLine{1},'linestyle','-');
            set(plot1(2),'color',ColorLine{2},'linestyle','--');
            hold on
            bandplotColor(DatesNum,quantile(Save_Collect.mubar(SelectTiming,:),[Quantili(1) Quantili(3)],2)',0,ColorShade{1}); 
            hold on
            plot1 = plot(DatesNum,[quantile(Save_Collect.mubar(SelectTiming,:),Quantili(2),2)],'linewidth',SPESSORE+.5,'color',ColorLine{1}); 
            ylim([0 .1]); 
            xlim([min(DatesNum) max(DatesNum)]); 
            hold on
            bandplotColor(DatesNum,quantile(Save_Collect.gbar(SelectTiming,:),[Quantili(1) Quantili(3)],2)',0,ColorShade{2}); 
            hold on
            plot1 = plot(DatesNum,[quantile(Save_Collect.gbar(SelectTiming,:),Quantili(2),2)],'linewidth',SPESSORE+.5,'color',ColorLine{2},'linestyle','--'); 
            ylim([0 .1]); 
            datetick
            recessionplot
            xlim([min(DatesNum) max(DatesNum)]); 
            grid on
            legend1 = legend('$\overline{\mu}_{t|t-1}$','$\overline{g}_{t|t-1}$'); 
            set(legend1,'Edgecolor',[1 1 1],'location','northeast','Interpreter','Latex','fontsize',FONTsizePick);
            legend('boxoff');
%             title([NamesPlot{jjii}]);
            hold off
            tightfig
            fig = gcf;
            fig.PaperPositionMode = 'auto'; 
%             fig.PaperSize = [5.85 4.5];
            print(fig,[WhereToSave 'Plot_Figure_2_a'],'-dpdf');



            figure
            plot1 = plot(DatesNum,NaN*zeros(TT,2),'linewidth',SPESSORE); 
            set(plot1(1),'color','k','linestyle','-');
            set(plot1(2),'color',ColorLine{3},'linestyle','--');
            hold on
            bandplotColor(DatesNum,quantile(Save_Collect.pdbar(SelectTiming,:),[Quantili(1) Quantili(3)],2)',0,ColorShade{3}); 
            hold on
            plot1 = plot(DatesNum,[quantile(Save_Collect.pdbar(SelectTiming,:),Quantili(2),2)],'linewidth',SPESSORE+.5,'color',ColorLine{3},'linestyle','--'); 
            ylim([2 4.5]); 
            xlim([min(DatesNum) max(DatesNum)]); 
            hold on
            plot1 = plot(DatesNum,y(:,2),'linewidth',SPESSORE,'color','k','linestyle','-'); 
            ylim([2 4.5]); 
            datetick
            recessionplot
            xlim([min(DatesNum) max(DatesNum)]); 
            grid on
            legend1 = legend('$pd_t$','$\overline{pd}_{t|t-1}$'); 
            set(legend1,'Edgecolor',[1 1 1],'location','southeast','Interpreter','Latex','fontsize',FONTsizePick);
            legend('boxoff');
%             title([NamesPlot{jjii}]);
            hold off
            tightfig
            fig = gcf;
            fig.PaperPositionMode = 'auto'; 
%             fig.PaperSize = [5.85 4.5];
            print(fig,[WhereToSave 'Plot_Figure_2_b'],'-dpdf');

%% PLOT EXPECTED DIVIDENDS AND RETURNS 

SelectTiming = [1:TT]; 
% SelectTiming = [2:TT+1]; 

Dividends = All_Data(:,3); 
ExpectedDividends = g_bar(SelectTiming) + [NaN;Res.RT_alfa_t(2,1:end-1)']; 
SD_ExpectedDividends = [NaN; sqrt(squeeze(Res.RT_P_t(2,2,1:end-1)))]; 
CI_ExpectedDividends = ExpectedDividends*[1 1]+SD_ExpectedDividends*[-1 1]; 

Returns = All_Data(:,1); 
ExpectedReturns = mu_bar(SelectTiming) + [NaN;Res.RT_alfa_t(3,1:end-1)']; 
SD_ExpectedReturns = [NaN; sqrt(squeeze(Res.RT_P_t(3,3,1:end-1)))]; 
CI_ExpectedReturns = ExpectedReturns*[1 1]+SD_ExpectedReturns*[-1 1]; 


            figure
            plot1 = plot(DatesNum,NaN*zeros(TT,3),'linewidth',SPESSORE); 
            set(plot1(1),'color','k','linestyle','-');
            set(plot1(2),'color',ColorLine{1},'linestyle',':');
            set(plot1(3),'color',ColorLine{1},'linestyle','--');
            hold on
            bandplotColor(DatesNum,CI_ExpectedReturns',-1,ColorShade{1}); 
            hold on
            plot1 = plot(DatesNum,[mu_bar(SelectTiming) ExpectedReturns],'linewidth',SPESSORE+.5,'color',ColorLine{1},'linestyle','--'); 
            ylim([-.4 .7]); 
            set(plot1(1),'color',ColorLine{1},'linestyle',':','linewidth',SPESSORE);
            set(plot1(2),'color',ColorLine{1},'linestyle','--');
            xlim([min(DatesNum) max(DatesNum)]); 
            hold on
            plot1 = plot(DatesNum,Returns,'linewidth',SPESSORE,'color','k','linestyle','-'); 
            ylim([-.4 .7]); 
            datetick
            recessionplot
            xlim([min(DatesNum) max(DatesNum)]); 
%             grid on
            legend1 = legend('$r_t$','$\overline{\mu}_{t|t-1}$','$\mathtt{E}_{t-1}(r_t)$'); 
            set(legend1,'Edgecolor',[1 1 1],'location','northeast','Interpreter','Latex','fontsize',FONTsizePick);
            legend('boxoff');
            title('Expected Return');
            hold off
            tightfig
            fig = gcf;
            fig.PaperPositionMode = 'auto'; 
%             fig.PaperSize = [5.85 4.5];
            print(fig,[WhereToSave 'Plot_Figure_1_a'],'-dpdf');

            
            figure
            plot1 = plot(DatesNum,NaN*zeros(TT,3),'linewidth',SPESSORE); 
            set(plot1(1),'color','k','linestyle','-');
            set(plot1(2),'color',ColorLine{2},'linestyle',':');
            set(plot1(3),'color',ColorLine{2},'linestyle','--');
            hold on
            bandplotColor(DatesNum,CI_ExpectedDividends',-1,ColorShade{2}); 
            hold on
            plot1 = plot(DatesNum,[g_bar(SelectTiming) ExpectedDividends],'linewidth',SPESSORE+.5,'color',ColorLine{2},'linestyle','--'); 
            ylim([-.5 .5]); 
            set(plot1(1),'color',ColorLine{2},'linestyle',':','linewidth',SPESSORE);
            set(plot1(2),'color',ColorLine{2},'linestyle','--');
            xlim([min(DatesNum) max(DatesNum)]); 
            hold on
            plot1 = plot(DatesNum,Dividends,'linewidth',SPESSORE,'color','k','linestyle','-'); 
            ylim([-.5 .5]); 
            datetick
            recessionplot
            xlim([min(DatesNum) max(DatesNum)]); 
%             grid on
            legend1 = legend('$\Delta d_t$','$\overline{g}_{t|t-1}$','$\mathtt{E}_{t-1}(\Delta d_t)$'); 
            set(legend1,'Edgecolor',[1 1 1],'location','northeast','Interpreter','Latex','fontsize',FONTsizePick);
            legend('boxoff');
            title('Expected Dividend Growth');
            hold off
            tightfig
            fig = gcf;
            fig.PaperPositionMode = 'auto'; 
%             fig.PaperSize = [5.85 4.5];
            print(fig,[WhereToSave 'Plot_Figure_1_b'],'-dpdf');
            
            
%% CONDITIONAL VARIANCE AND CORRELATION OF RETURNS 

[Decompositions, Expectations] = Decomposition_TVKF_pd_Pd_Model(Res,100); 


%% PLOT EXPECTATIONS             
SelectYears = [2,10];
QualeMaturityAll = [1:10];

LevelAndSlope_G = [Expectations.Annualized_Exp_g(:,SelectYears(2))...
    Expectations.Annualized_Exp_g(:,SelectYears(2))-Expectations.Annualized_Exp_g(:,SelectYears(1))]; 
LevelAndSlope_MU = [Expectations.Annualized_Exp_mu(:,SelectYears(2))...
    Expectations.Annualized_Exp_mu(:,SelectYears(2))-Expectations.Annualized_Exp_mu(:,SelectYears(1))]; 


            figure
            plot1 = plot(DatesNum,[LevelAndSlope_G],'linewidth',SPESSORE+.5);
            set(plot1(1),'color',ColorLine{1},'linestyle','--');
            set(plot1(2),'color',ColorLine{2},'linestyle','-');
            datetick
            recessionplot
            xlim([min(DatesNum) max(DatesNum)]); 
            grid on
            legend1 = legend('Level: $g^{(10)}_t$','Slope: $g^{(10)}_t-g^{(2)}_t$'); 
            set(legend1,'Edgecolor',[1 1 1],'location','northeast','Interpreter','Latex','fontsize',FONTsizePick-5);
            legend('boxoff');
            title(['Term Structure of Expected Dividend Growth']);
            hold off
            tightfig
            fig = gcf;
            fig.PaperPositionMode = 'auto'; 
%             fig.PaperSize = [5.85 4.5];
            print(fig,[WhereToSave 'Plot_Figure_4_b'],'-dpdf');


            figure
            plot1 = plot(DatesNum,[LevelAndSlope_MU],'linewidth',SPESSORE+.5);
            set(plot1(1),'color',ColorLine{1},'linestyle','--');
            set(plot1(2),'color',ColorLine{2},'linestyle','-');
            datetick
            recessionplot
            xlim([min(DatesNum) max(DatesNum)]); 
            grid on
            legend1 = legend('Level: $\mu^{(10)}_t$','Slope: $\mu^{(10)}_t-\mu^{(2)}_t$'); 
            set(legend1,'Edgecolor',[1 1 1],'location','southeast','Interpreter','Latex','fontsize',FONTsizePick-5);
            legend('boxoff');
            title(['Term Structure of Expected Return']);
            hold off
            tightfig
            fig = gcf;
            fig.PaperPositionMode = 'auto'; 
%             fig.PaperSize = [5.85 4.5];
            print(fig,[WhereToSave 'Plot_Figure_4_a'],'-dpdf');

            
            
%% PD Decomposition 
% A) MEAN 


            figure
            plot1 = plot(DatesNum,[[Fit_Pd(:,3:4)]],'linewidth',SPESSORE+.5);
            set(plot1(1),'color',ColorLine{1},'linestyle','--');
            set(plot1(2),'color',ColorLine{2},'linestyle','-');
            datetick
            recessionplot
            xlim([min(DatesNum) max(DatesNum)]); 
            grid on
            legend1 = legend('Return Component (i.e. $-b_{1,t|t-1}\tilde{\mu}_t)$','Dividend Growth Component (i.e. $b_{2,t|t-1}\tilde{g}_t)$'); 
            set(legend1,'Edgecolor',[1 1 1],'location','northwest','Interpreter','Latex','fontsize',FONTsizePick-5);
            legend('boxoff');
            title1 = title('$pd_t - \overline{pd}_{t|t-1}$');
            set(title1,'Interpreter','Latex');
            hold off
            tightfig
            fig = gcf;
            fig.PaperPositionMode = 'auto'; 
%             fig.PaperSize = [5.85 4.5];
            print(fig,[WhereToSave 'Plot_Figure_3_a'],'-dpdf');

            
% B: Conditional VARIANCE 

% Share equally the covariance component 
PD_CondVarDec = Decompositions.Pd_CondVarDecomp(:,[2,3])+Decompositions.Pd_CondVarDecomp(:,4)*ones(1,2)*.5; 
PD_CondVarDec = [Decompositions.Pd_CondVarDecomp(:,[1]) PD_CondVarDec]; 

            figure
            bar1 = bar(DatesNum,[PD_CondVarDec(:,[2,3])],'BarWidth',1.2,'BarLayout','stacked');
            set(bar1(2),'FaceColor',ColorShade{1},'EdgeColor',ColorShade{1});
            set(bar1(1),'FaceColor',[1 0 0],'EdgeColor',[1 0 0]);
            hold on 
            plot1 = plot(DatesNum,[PD_CondVarDec(:,1)],'linewidth',SPESSORE,'color','k'); 
            datetick
            recessionplot
            xlim([min(DatesNum) max(DatesNum)]); 
            grid on
            legend1 = legend('Dividend Growth Component','Return Component'); 
            set(legend1,'Edgecolor',[1 1 1],'location','northwest','Interpreter','Latex','fontsize',FONTsizePick-5);
            legend('boxoff');
            title1 = title('$Var_t(pd_{t+1} - \overline{pd}_{t+1|t})$');
            set(title1,'Interpreter','Latex');
            hold off
            tightfig
            fig = gcf;
            fig.PaperPositionMode = 'auto'; 
%             fig.PaperSize = [5.85 4.5];
            print(fig,[WhereToSave 'Plot_Figure_3_b'],'-dpdf');

            