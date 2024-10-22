% specify likelihood function under a particular Vpar
function [LL,Res] = kf_simple_mymodel(Vpar,yy)

% get information on the problem
N           = size(yy,1);               % Total no. of variables
TT          = size(yy,2);               % Observations to be used
m           = 10;                        % dimension of the state space (1 unobserved components)

run Useful_Transformations

%% 
if size(Vpar,1) ==1 
    Vpar = Vpar'; 
end

% ar coefficients
psi_u = MaxZero(1,Vpar(1));
psi_p = MaxZero(1,Vpar(2));
psi_n = MaxZero(1,Vpar(3));
psi_t = MaxZero(1,Vpar(4));
theta_u = MaxZero(-1,Vpar(5));
theta_p = MaxZero(-1,Vpar(6));

% error mean
% mean_err = Vpar(7:9); 
mean_err = [0,0,0]';

% error covariance
% H = diag(Vpar(6:7).^2);
% Q = [Vpar(8)^2, 0; 0, 0];
% cov_half = [Vpar(10),0,0;Vpar(11),Vpar(13),0;Vpar(14),Vpar(15),Vpar(12)];
cov_half = [Vpar(7),0,0;Vpar(10),Vpar(8),0;Vpar(11),Vpar(12),Vpar(9)];
Q = cov_half*(cov_half');
H = zeros(2,2);

% system matrices
Z = [1,zeros(1,4),1, zeros(1,4);zeros(1,2),1,zeros(1,3),1, zeros(1,3)];
T = blkdiag([1+psi_u, -psi_u;1, 0],[psi_n+psi_u+1, -(psi_u*psi_n+psi_u+psi_n), psi_u*psi_n; 1,0,0;0,1,0], [psi_p],[psi_n+psi_t,-psi_n*psi_t;1,0],zeros(2,2));
T(1,8)=1;
T(3,9)=1-psi_n;
T(6,10)=1;
SS = [1,0,0;0,0,0;1-psi_n,0,0;zeros(2,3);0,1,0;0,0,1;0,0,0;theta_u,0,0;0,theta_p,0];
  
%% INITIALIZE THE KF

% starting values for state vectors and state variances
a_1 = zeros(m,1);
kap = 1e9;

% selection for nonstationary part of the state vector, dimension m*(number
% of nonstationary variables)
Ainit = [eye(5);zeros(5,5)];
Rinit = [zeros(5,5);eye(5)];
Pinf = Ainit*Ainit';

% selection for stationary part of the state vector, dimension m*(number of
% stationary variables) Pstar 0 in our case
initOmega = SS*Q*(SS');

initP = ((eye(m^2)-kron(T+(1e-6)*eye(m),T+(1e-6)*eye(m)))\eye(m^2))*initOmega(:);
initP = reshape(initP,m,m);
P_1 = kap*Pinf+Rinit*initP(6:10,6:10)*Rinit';

alfa_t = [a_1 NaN*zeros(m,TT)];
P_t    = [P_1(:) NaN*zeros(m^2,TT)];
Lik_t  = NaN*zeros(1,TT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORGANIZE OUTPUT

Res.Z = cell(TT);
Res.v_t   = cell(TT);
Res.invF_t   = cell(TT);
Res.L_t = NaN(m,m,TT);
Res.alfa_t   = NaN(m,TT+1);
Res.P_t      = NaN(m,m,TT+1);
        
Res.alfa_t(:,1) = a_1; 
Res.P_t(:,:,1) = P_1;

% alpha_t here is at|t-1 in Harvey's book

for tt=1:TT+1

    
%% Kalman Filter 

    P = reshape(P_t(:,tt),m,m);
    a_t = alfa_t(:,tt);
    
    if tt<=TT
        
        %% STEP 1: KF STEP
        yt = yy(:,tt);
        [Wt,nt]=SelectMatW(yt);
        Wt = sparse(Wt);
        yt(isnan(yt)) = 0;
        Res.Z{tt} = full(Wt)*Z;
        
        v = full(Wt)*(yt - Z*a_t);
        Res.v_t{tt} = v;
        
        F =  full(Wt)*(Z*P*Z' + H)* full(Wt)';

        invF = F\eye(size(F));
        Res.invF_t{tt} = invF;

        G = sparse(P*Z'*full(Wt)'*invF);

        %% STEP 2: UPDATE Likelihood
        Lik_t(1,tt) = -.5*(N*log(2*pi) + log(det(F)) + v'*invF*v);
            
        %% UPDATE KF
        K = T*full(G); L = (T-K*full(Wt)*Z);  
        Res.L_t(:,:,tt) = L;

        alfa_t(:,tt+1)          = T*a_t + K*v + SS*mean_err;
        Res.alfa_t(:,tt+1)      = alfa_t(:,tt+1);

        % P_temp                  = T*P*L'+Q;
        Omega = SS*Q*SS';
        P_temp                  = L*P*L'+K*full(Wt)*H*full(Wt)'*K'+ Omega;
        Res.P_t(:,:,tt+1)       = P_temp;
        P_t(:,tt+1)             = P_temp(:);

        
    end
end

Res.ContributionLogLik = Lik_t;
cutFirstObservations = 0; 
LL = sum(Lik_t(1,cutFirstObservations+1:end));


if isreal(LL)==0 
    LL = -10^10; 
end

if isfinite(LL)==0 
    LL = -10^10; 
end
