function [TransParam,Jacob] = TransformParametersLik(Vpar);

run Useful_Transformations

%% 
if size(Vpar,1) ==1 
    Vpar = Vpar'; 
end

%% SET UP THE DYNAMICS OF THE SCORE DRIVEN PARAMTERS 

% a) Parameters of the score driven process 
TransParam(1:5) = [Vpar(1:5);]; %% OMEGA 
Jacob(1:5) = ones(5,1); 

TransParam(6:10) = MaxZero(1,Vpar(6:10)); %% A  
Jacob(6:10) = DerMaxZero(1,Vpar(6:10)); 

TransParam(11:17) = [[.015+MaxZero(.3,Vpar(11:12)); MaxZero(.1,Vpar(13:17))]]; %% B
Jacob(11:17) = [[DerMaxZero(.3,Vpar(11:12)); DerMaxZero(.1,Vpar(13:17))]]; 

TransParam(18)    = MaxZero(.5,Vpar(18)); %% kap_hes
Jacob(18) = DerMaxZero(.5,Vpar(18)); 

TransParam(19) = PositiveTrans(Vpar(19)); %% Sigma_Meas_pd
Jacob(19) = DerPositiveTrans(Vpar(19));

TransParam(20:21) = MaxZero(.99,Vpar([20,21])); %% ARcoeffs_Trans
Jacob(20:21) = DerMaxZero(.99,Vpar([20,21]));



