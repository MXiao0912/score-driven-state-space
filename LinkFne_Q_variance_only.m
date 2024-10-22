function [Qt,Qdot] = LinkFne_Q_variance_only(TVparam)

% CaseTransformationVariance = 'Exponential'; %'Exponential'; %'Square'; 'Logistic';
MinVol = 10^-5; 
PositiveTrans = @(x) MinVol + exp(x);
DerPositiveTrans = @(x) exp(x);        

% MaxZero = @(maxVal,x) MinVol + maxVal*(1./(1+exp(-x)));
% cap_vol = 5; 
% PositiveTrans = @(x) MaxZero(cap_vol,x);
% DerMaxZero = @(maxVal,x) maxVal*MaxZero(1,x).*(1-MaxZero(1,x));
% DerPositiveTrans = @(x) DerMaxZero(cap_vol,x);

Ssigma=[1 0 0;...
        zeros(3,3);...
        0 1 0;...
        zeros(3,3);...
        0 0 1];                 %S1d
Sdelta =eye(3);   %S2d
deltas = Sdelta*TVparam; 
D = diag(PositiveTrans(deltas)); 

R = eye(size(D,1));
Qt = D*R*D;

II = eye(size(Qt,1));

dvecQ_dvecD = kron(D*R,II) + kron(II,D*R);
dvecD_dSigma = Ssigma; 
dSigma_ddelta = diag(DerPositiveTrans(deltas)); 
ddelta_df = Sdelta; 


Qdot = dvecQ_dvecD*dvecD_dSigma*dSigma_ddelta*ddelta_df;






