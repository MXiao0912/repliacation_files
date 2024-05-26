

%% USEFUL TRANSFORMATIONS 
PositiveTrans = @(x) exp(x); 
InvPositiveTrans = @(x) log(x); 
DerPositiveTrans = @(x) exp(x);

UnoMenoUno = @(x) tanh(x); 
InvUnoMenoUno = @(x) atanh(x);
DerUnoMenoUno = @(x) 1-(tanh(x)^2);

% MaxZero = @(maxVal,x) maxVal*(x^2/(x^2+1)); InvMaxZero = @(maxVal,x) sqrt(x/(maxVal-x)); 
% MaxZero = @(maxVal,x) maxVal*0.5*(tanh(x)+1); InvMaxZero = @(maxVal,x) atanh((2*x-maxVal)./maxVal);
LogisticFun = @(x) (1./(1+exp(-x))); 
MaxZero = @(maxVal,x) maxVal*(1./(1+exp(-x))); 
InvMaxZero = @(maxVal,x) -log((maxVal./x)-1);
DerMaxZero = @(maxVal,x) maxVal*LogisticFun(x).*(1-LogisticFun(x));

%% SOME MATRIX OPERATIONS 
vec = @(x) x(:);




