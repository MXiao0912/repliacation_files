function [SigmaMat,StdMat,CorrMat,PCorrMat] = Variance_DR_Partial_Mapping3by3(SigVect,PiVect)

%% Define size 
nn = 3; 

%% Rearrange the partial correlations 
A = tril(ones(nn),-1);
A(A~=0) = PiVect; 
PCorrMat = A + A' +eye(nn);


%% Map the partial corr to corr 
RhoVect = PiVect; 
RhoVect(end) = PiVect(end)*sqrt((1-PiVect(1)^2)*(1-PiVect(2)^2))+PiVect(1)*PiVect(2);

%% Rearrange the correlations 
A = tril(ones(nn),-1);
A(A~=0) = RhoVect; 
CorrMat = A + A' +eye(nn);

%% Rearrange the Std 
StdMat = diag(SigVect); 

%% Compute Variance 
SigmaMat = StdMat*CorrMat*StdMat';
