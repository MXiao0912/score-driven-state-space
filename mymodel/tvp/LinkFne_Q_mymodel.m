function [Qt,Qdot] = LinkFne_Q(TVparam)

% CaseTransformationVariance = 'Exponential'; %'Exponential'; %'Square'; 'Logistic';
MinVol = 10^-5; 
PositiveTrans = @(x) MinVol + exp(x);
DerPositiveTrans = @(x) exp(x);        

S1d = [1 0 0;...
        zeros(3,3);...
        0 1 0;...
        zeros(3,3);...
        0 0 1];  
S2d = [eye(3) zeros(3,3)];
S1r = [0 0 0; ...
        1 0 0;...
        0 1 0;...
        1 0 0;...
        0 0 0;...
        0 0 1;...
        0 1 0;...
        0 0 1;...
        0 0 0];
S2r = [zeros(3,3) eye(3)];
deltas = S2d*TVparam; 
D = diag(PositiveTrans(deltas)); 

gammas = S2r*TVparam;
[R,Pi,Psi12,Psi31,Psi32] = MappingGammaToCorrelations(gammas);

Qt = D*R*D;

II = eye(size(Qt,1));

DRIIDR = kron(D*R,II) + kron(II,D*R);
Ddot = S1d*diag(DerPositiveTrans(deltas))*S2d;
DD = kron(D,D);

drho_dpi = [1 0 0;0 1 0;...
Pi(1,3)-Pi(1,2)*Pi(2,3)*sqrt((1-Pi(1,3)^2)/(1-Pi(1,2)^2)),...
Pi(1,2)-Pi(2,3)*Pi(1,3)*sqrt((1-Pi(1,2)^2)/(1-Pi(1,3)^2)),...
sqrt((1-Pi(1,3)^2)*(1-Pi(1,2)^2))];

dpi_dgamma = diag([Psi12,Psi31,Psi32]);
Rdot = S1r*drho_dpi*dpi_dgamma*S2r;

Qdot = DRIIDR*Ddot+DD*Rdot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_t,Pi_t,Psi12,Psi31,Psi32] = MappingGammaToCorrelations(gammas)

smoothConstantCorr = 10^-1;%10^-1;%1;%10^-2; 
MaxValueCorr = .95; 

CorrTrans = @(x) MaxValueCorr*tanh(smoothConstantCorr*x);
DerCorrTrans = @(x) MaxValueCorr*smoothConstantCorr*(1-x.^2);   

Pi_vect = CorrTrans(gammas)';
Pi_t = eye(3) + [zeros(1,3);Pi_vect(1) 0 0; Pi_vect(2:3) 0] +[zeros(1,3);Pi_vect(1) 0 0; Pi_vect(2:3) 0]';

Psi12 =DerCorrTrans(Pi_t(1,2));
Psi31 =DerCorrTrans(Pi_t(3,1)); 
Psi32 =DerCorrTrans(Pi_t(3,2));

rho_12 = Pi_t(1,2); 
rho_31 = Pi_t(3,1); 
rho_32 = Pi_t(3,2)*sqrt((1-Pi_t(1,2)^2)*(1-Pi_t(3,1)^2))+Pi_t(1,2)*Pi_t(3,1);

R_t = eye(3) + [zeros(1,3);rho_12 0 0;[rho_31 rho_32] 0] + [zeros(1,3);rho_12 0 0;[rho_31 rho_32] 0]';


