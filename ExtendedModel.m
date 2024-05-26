
clear; clc; 


% FileName_ResultsSave = 'TVP_ExRetModel_Estimates_10Sept';
FileName_ResultsSave = 'TVP_ExRetModel_Estimates_10Sept_70obs_init';
% FileName_ResultsSave = 'TVP_ExRetModel_Estimates_10Sept_85obs_init';

%% DECIDE IF ESTIMATE OR USE SAVED DATA 
EstimateModel = 'N'; 


%% LOAD THE ADDITIONAL DATA YOU NEED 
%% FROM BASELINE MODEL AND r* 

load('TVP_BAselineModel_Estimates_10Sept_Prova.mat','Save_Collect');  
Quantili = [.16 .5 .84];

mubar_Baseline= [Save_Collect.mubar(:,1) quantile(Save_Collect.mubar,[Quantili([2,1,3])],2)];
gbar_Baseline= [Save_Collect.gbar(:,1) quantile(Save_Collect.gbar,[Quantili([2,1,3])],2)];
pdbar_Baseline= [Save_Collect.pdbar(:,1) quantile(Save_Collect.pdbar,[Quantili([2,1,3])],2)];
clear Save_Collect; 

WhereAmI = cd; 
load([WhereAmI '\Data\Data_Rstar.mat']); 


%%

CaseToEstimate = 'AllCoeffs'; %'AllCoeffs'; % 'MinimumCoeffs'; 'MediumCoeffs'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CODES FOR TVP MODEL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SettingsTransfor.CaseTransformationMu = 'Exponential'; % 'NoTransformation';%'Exponential'; %'Square'; 'Logistic';
SettingsTransfor.MaxMu = .1;
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
% load([WhereAmI '\Data\SaveDataReal_' PickDataFreq '.mat']); 
load([WhereAmI '\Data\SaveDataPremium_' PickDataFreq '.mat']); 

%% PICK DATA 
%  y = [d, pd]

SelectVariables = [3,2];
StartDate =find(Dates == 1873);% 1;% 
EndDate = find(Dates == 2018); %size(All_Data,1); %
y = All_Data(StartDate:EndDate,SelectVariables); 
y_Names = All_Names(1,SelectVariables)'


%% USEFUL TRANSFORMATIONS 
run Useful_Transformations

%% SET INITIAL VALUES TIME VARYING ELEMENTS 
% f = [pd, g, ln(sigma_d), ln(sigma_g), ln(sigma_mu),atanh(pi_dmu),atanh(pi_gmu)]
QUANTI_InitObs = 75; 
Initialize_TVP

% initialF1(1) = .03; 
switch CaseTransformationMu 
    case 'Exponential'
        initialF1(1) = InvPositiveTrans(initialF1(1)); 
    case 'Square' 
        initialF1(1) = sqrt(initialF1(1)); 
    case 'Logistic'
        ChosenCap = SettingsTransfor.MaxMu;
        initialF1(1) = InvMaxZero(ChosenCap,initialF1(1)); 
end

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

Init_fi1 = .4; %.9
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
        LossToMinmize = @(vparam) -TVKF_pd_Pd_ExRetModel([vparam],y',initialF1,[],SettingsTransfor)/size(y,1);
    case 'MinimumCoeffs'
        LossToMinmize = @(vparam) -TVKF_pd_Pd_ExRetModel([vparam(1:5,1); ...
                        vparam(6,1)*ones(3,1);vparam(7,1)*ones(2,1);vparam(8,1)*ones(2,1);...
                        vparam(9,1)*ones(3,1);vparam(10,1)*ones(2,1);vparam(11:14,1)],y',initialF1,[],SettingsTransfor)/size(y,1);
    case 'MediumCoeffs'
        LossToMinmize = @(vparam) -TVKF_pd_Pd_ExRetModel([vparam(1:5,1); ...
                        vparam(6,1)*ones(3,1);vparam(7,1)*ones(2,1);vparam(8,1);vparam(9,1);...
                        vparam(10,1)*ones(3,1);vparam(11,1)*ones(2,1);vparam(12:15,1)],y',initialF1,[],SettingsTransfor)/size(y,1);
end

             
%% MAXIMIZE LIKELIHOOD

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

[LL,Res] = TVKF_pd_Pd_ExRetModel(EstimParamsNew,y',initialF1,[],SettingsTransfor); 

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

% %% SOME PRELIMINARY PLOTS 
% Dates_prova = [1872; Dates]; 
% figure
% plot(Dates_prova,[pd_bar pdbar_Baseline(:,1)]);
% figure
% plot(Dates_prova,[mu_bar mubar_Baseline(:,1)]);
% figure
% plot(Dates_prova,[mu_bar]);
% Implied_Real_Rates = [mubar_Baseline(:,1)-mu_bar gbar_Baseline(:,1)-g_bar];
% figure
% plot(Dates_prova,Implied_Real_Rates);
% figure
% plot(Dates_prova,[Implied_Real_Rates [NaN*RealRfRate(1,:); RealRfRate]]);
% figure
% plot(Dates_prova,[Implied_Real_Rates [NaN*Data_Rstar_A(1,:); Data_Rstar_A]]);
% 
% paramTransform(end-1:end)

% return
%% CALCULATE SE ESTIMATES WITH SIMULATIONS (coded up only 'AllCoeffs' specification)

LossToMinmize = @(vparam) TVKF_pd_Pd_ExRetModel([vparam],y',initialF1,[],SettingsTransfor);
ContrLossToMinmize = @(vparam) ContributionLogLik_ExReturns([vparam],y',initialF1,[],SettingsTransfor);
TT = size(y,1); 

% HHessian = fdhess2(LossToMinmize, EstimParamsNew );
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
    
    [LL_i,Res_i] = TVKF_pd_Pd_ExRetModel(Paramter_i,y',initialF1,[],SettingsTransfor); 
    
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

xlswrite('EstimatesModelFinal',TableMainEstimates,'TVP_Model_ExRet','c4'); 
xlswrite('EstimatesModelFinal',TableGASEstimates,'TVP_Model_ExRet','h4'); 

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
ColorLine{4} = 'y';  

ColorShade{1} = [.7 .9 1];  
ColorShade{2} = [1 .7 .8]; 
ColorShade{3} = [.8 1 .8];  
ColorShade{4} = [1,1,.8];




%% Calculate the implied r*            
HowManyPeriods = 10; 
MA_RealRate = tsmovavg(RealRfRate','s',HowManyPeriods)'; 

r_star_from_mu = mubar_Baseline(:,1)*ones(1,size(Save_Collect.mubar,2)) - Save_Collect.mubar; 
r_star_from_g = gbar_Baseline(:,1)*ones(1,size(Save_Collect.gbar,2)) - Save_Collect.gbar; 
r_star_ave = (r_star_from_mu + r_star_from_g)*.5; 

% Which_Rstar_toPlot = r_star_from_mu; 
% Which_Rstar_toPlot = r_star_from_g; 
Which_Rstar_toPlot = r_star_ave; 

            

            figure
            plot1 = plot(DatesNum,NaN*zeros(TT,2),'linewidth',SPESSORE); 
            set(plot1(1),'color',ColorLine{1},'linestyle','-');
            set(plot1(2),'color',ColorLine{2},'linestyle','--');
            hold on
            plotConfidenceBandsBlue(DatesNum,[mubar_Baseline(2:end,[3,2,4])],ColorLine{1}); 
            hold on
            plotConfidenceBandsBlueBroken(DatesNum,[quantile(Save_Collect.mubar(SelectTiming,:),Quantili,2)],ColorLine{2}); 
%             ylim([-0.04 0.08]); 
            ylim([0 .1]); 
            xlim([min(DatesNum) max(DatesNum)]); 
            hold on
            datetick
            recessionplot
            xlim([min(DatesNum) max(DatesNum)]); 
            grid on
            legend1 = legend('$\overline{\mu}_{t|t-1}$','$\overline{\mu}^{ex}_{t|t-1}$'); 
            set(legend1,'Edgecolor',[1 1 1],'location','northeast','Interpreter','Latex','fontsize',FONTsizePick);
            legend('boxoff');
            title('LONG-RUN RETURN AND EQUITY PREMIUM');
            hold off
            tightfig
            fig = gcf;
            fig.PaperPositionMode = 'auto'; 
%             fig.PaperSize = [5.85 4.5];
            print(fig,[WhereToSave 'Plot_Figure_5_a'],'-dpdf');
            
            figure
            plot1 = plot(DatesNum,NaN*zeros(TT,3),'linewidth',SPESSORE+1); 
            set(plot1(1),'color',ColorLine{1},'linestyle','-');
            set(plot1(2),'color',ColorLine{2},'linestyle','--');
            set(plot1(3),'color',ColorLine{3},'linestyle',':');
            hold on
            bandplotColor(DatesNum,quantile(Which_Rstar_toPlot(SelectTiming,:),[Quantili(1) Quantili(3)],2)',0,ColorShade{1}); 
            hold on
            plot1 = plot(DatesNum,[quantile(Which_Rstar_toPlot(SelectTiming,:),Quantili(2),2)],'linewidth',SPESSORE+.5,'color',ColorLine{1},'linestyle','-'); 
            hold on
            plot1 = plot(DatesNum,Data_Rstar_A(SelectTiming,1),'linewidth',SPESSORE+1,'color',ColorLine{2},'linestyle','--'); 
%             plot1 = plot(DatesNum,Data_Rstar_A(SelectTiming,3),'linewidth',SPESSORE+1,'color',ColorLine{2},'linestyle','--'); 
            hold on
            plot1 = plot(DatesNum,Data_Rstar_A(SelectTiming,3),'linewidth',SPESSORE+1,'color',ColorLine{3},'linestyle',':'); 
%             plot1 = plot(DatesNum,Rbar_DGGT(SelectTiming,3),'linewidth',SPESSORE+1,'color',ColorLine{3},'linestyle',':'); 
            hold on
            ylim([0 0.06]); 
            datetick
            recessionplot
            xlim([min(DatesNum) max(DatesNum)]); 
            grid on
            legend1 = legend('$\overline{r}_{t|t-1}$','$r^{*}_t$ (LW)','$r^{*}_t$ (HLW)'); 
%             legend1 = legend('$\bar{r}_t$','$r^{*}_t$ (LW)','$r^{*}_t$ (DDGT)'); 
            set(legend1,'Edgecolor',[1 1 1],'location','southwest','Interpreter','Latex','fontsize',FONTsizePick);
            legend('boxoff');
            title('LONG-RUN RISKLESS REAL RATE');
            hold off
            tightfig
            fig = gcf;
            fig.PaperPositionMode = 'auto'; 
%             fig.PaperSize = [5.85 4.5];
            print(fig,[WhereToSave 'Plot_Figure_5_b'],'-dpdf');
            
            
