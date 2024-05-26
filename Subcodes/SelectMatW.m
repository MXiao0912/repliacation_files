function [Wt,nt]=SelectMatW(yt) 
% yt - is an Nx1 Vector
%
% Wt - is NtxN


Select = isfinite(yt);
nt = sum(Select);
II = eye(size(yt,1));
Wt = II(Select,:);

