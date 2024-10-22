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

