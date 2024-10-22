function [Qt,Qdot] = LinkFne_Q(TVparam)

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
Sdelta =[eye(3) zeros(3,2)];   %S2d
Srho = [zeros(2,2); ...
        1 0;...
        zeros(2,2);...
        0 1;...
        1 0;...
        0 1;...
        0 0];                   %   S1r
Sgamma = [zeros(2,3) eye(2)]; %S2r
deltas = Sdelta*TVparam; 
D = diag(PositiveTrans(deltas)); 

gammas = Sgamma*TVparam;
[R,Pi,Psi31,Psi32] = MappingGammaToCorrelations(gammas);

Qt = D*R*D;

II = eye(size(Qt,1));

dvecQ_dvecD = kron(D*R,II) + kron(II,D*R);
dvecD_dSigma = Ssigma; 
dSigma_ddelta = diag(DerPositiveTrans(deltas)); 
ddelta_df = Sdelta; 

dvecQ_dvecR = kron(D,D);
dvecR_drho = Srho; 
dgamma_df = Sgamma; 

drho_dpi = sqrt(Psi31);
dpi_dgamma = [sqrt(Psi31), 0; -Pi(3,1)*Pi(3,2), Psi32];

Qdot = dvecQ_dvecD*dvecD_dSigma*dSigma_ddelta*ddelta_df +...
        dvecQ_dvecR*dvecR_drho*drho_dpi*dpi_dgamma*dgamma_df;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_t,Pi_t,Psi31,Psi32] = MappingGammaToCorrelations(gammas)

smoothConstantCorr = 10^-1;%10^-1;%1;%10^-2; 
MaxValueCorr = .95; 

CorrTrans = @(x) MaxValueCorr*tanh(smoothConstantCorr*x);
DerCorrTrans = @(x) MaxValueCorr*smoothConstantCorr*(1-x.^2);   


Pi_vect = CorrTrans(gammas)';
Pi_t = eye(3) + [zeros(2,3); Pi_vect 0] +[zeros(2,3); Pi_vect 0]';

Psi31 =DerCorrTrans(Pi_t(3,1)); 
Psi32 =DerCorrTrans(Pi_t(3,2));

rho_31 = Pi_t(3,1); 
rho_32 = Pi_t(3,2)*sqrt(1-Pi_t(3,1)^2);

R_t = eye(3) + [zeros(2,3); [rho_31 rho_32] 0] + [zeros(2,3); [rho_31 rho_32] 0]';





