% specify likelihood function under a particular Vpar
disp("Fix Q bug: initialize P")
function [LL,Res] = kf_tvp_corr1(Vpar,yy,initialF1, pred_ft)
% as the state space takes as obs the differenced y, the obs size for the
% state space is T-1, which would correspond to true ft from 3:end as ft start 
% with a sample size of T+1 and the first 2 periods are dropped. For Q and
% alfa_t, and est ft, count from the second one as the first is the initial
% values.

clearFileContents('log.txt');

% get information on the problem
N           = size(yy,1);               % Total no. of variables
TT          = size(yy,2);               % Observations to be used
m           = 4;                        % dimension of the state space (1 unobserved components)
k           = 4;

run Useful_Transformations

%% 
if size(Vpar,1) ==1 
    Vpar = Vpar'; 
end

% ar coefficients
psi = UnoMenoUno(Vpar(1));
phi = UnoMenoUno(Vpar(2));

% error mean
mean_err = Vpar(3:5); 

% score transition parameters
kap_hes = MaxZero(1,Vpar(6));
c = Vpar(7:10);
A = diag(MaxZero(0.99,Vpar(11:14)));
B = diag(MaxZero(0.1,Vpar(15:18)));
% c = (eye(5)-A)*initialF1;

% system matrices
Z = [1, -psi, 1, 0; 1-phi, 0, 0, 1];
T = [1, 0, 0, 0; 1, 0, 0, 0;0 0 0 0;0 0 0 0]; 

%% Selection matrix for Q in the state equations
SS = [0 0 1;0 0 0;1 0 0; 0 1 0];

%% INITIALIZE THE RECURSIONS FOR THE SCORE DRIVEN MODEL 
[Q,Qdot] = LinkFne_Q_corr1(initialF1);
% prep the Q for step1 of ssm
Omega = SS*Q*SS';
kronSS = kron(SS,SS); 
Omegadot = kronSS*Qdot;

%% INITIALIZE THE KF

% starting values for state vectors and state variances
a_1 = zeros(m,1);
kap = 1e9;

% selection for nonstationary part of the state vector, dimension m*(number
% of nonstationary variables)
Ainit = [1 0; 0 1; zeros(2,2)];
Rinit = [zeros(2,2); 1 0; 0 1];
Pinf = Ainit*Ainit';

% selection for stationary part of the state vector, dimension m*(number of
% stationary variables) Pstar 0 in our case
P_1 = kap*Pinf+Rinit*Q(1:2,1:2)*Rinit';
  
yy = yy(:,2:end) - [psi;phi].*yy(:,1:(end-1)); % MODIFY FOR NAN VALUES
TT          = size(yy,2);               % Observations to be used

alfa_t = [a_1 NaN*zeros(m,TT)];
P_t    = [P_1(:) NaN*zeros(m^2,TT)];
Lik_t  = NaN*zeros(1,TT);

% initialize score
if ~isempty(pred_ft)
    f_t = pred_ft;
else
    f_1 = initialF1;
    f_t = [f_1 NaN*zeros(k,TT)];
end
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

Res.Q        = NaN(3,3,TT+1);
Res.Qdot        = NaN(9,4,TT+1);

Res.Q(:,:,1) = Q; 
Res.Qdot(:,:,1) = Qdot;

H = zeros(N);

Res.f_t      = f_t;

Res.fdot = cell(TT);

Res.F = cell(TT);
Res.score = NaN(k, TT);
Res.GradD = cell(TT);
Res.InfMat = cell(TT);
Res.InfMat_temp = cell(TT);

% alpha_t here is at|t-1 in Harvey's book

Score = zeros(4,1);
% v_score = zeros(3,1);
% mom = 0.9;
mom = MaxZero(1, Vpar(19));


% try
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
        v = Wt*(yt - Z*a_t);
        Res.Z{tt} = Wt*Z;
        Res.v_t{tt} = v;
        
        F =  Wt*(Z*P*Z' + H)* Wt';

        epsilon_F = 1e-6;
        F = F + epsilon_F * speye(size(F)); 

        % try
        % % Cholesky decomposition for stable inversion
            % R = chol(F);
            % invF = R \ (R' \ speye(size(F)));
        % catch
        %     % If Cholesky fails, fall back to pseudo-inverse
        %     invF = pinv(full(F));
        % end
        Res.F{tt} = F;
        invF = pinv(F);
        % invF = inv(F);

        Res.invF_t{tt} = invF;

        G = sparse(P*Z'*Wt'*invF);
        
        %% STEP 2: UPDATE Likelihood
        Lik_t(1,tt) = -.5*(N*log(2*pi) + log(det(F)) + v'*invF*v);

        %% STEP 3: COMPUTE SCORE
        
        % Nnt = sparse(.5*(eye(nt^2)+CommMatrixK(nt)));
        Fdot =  kron(Wt*Z,Wt*Z)*Omegadot;

        Res.fdot{tt} = Fdot;

        invFkroninvF=kron(invF, invF);
        % invFkroninvF = invFkroninvF + epsilon_F * speye(size(invFkroninvF));

        GradD = (Fdot'*(invFkroninvF)*(kron(v,v)-F(:)));
        Res.GradD{tt} = GradD;

        Res.F{tt} = F;
        if tt==1
            InfMat_temp = Fdot'*(invFkroninvF)*Fdot;
            InfMat = (1-kap_hes) * (eye(size(InfMat_temp))*1) + kap_hes * InfMat_temp;
            % InfMat = (1-kap_hes) * diag(diag(InfMat_temp)) + kap_hes * InfMat_temp;
        else
            InfMat_temp = Fdot'*(invFkroninvF)*Fdot;
            InfMat =  (1-kap_hes) * InfMat + kap_hes * InfMat_temp;   %smoothing the Hessian of coeffs
        end

        % InfMat_temp = Fdot'*(invFkroninvF)*Fdot;
        % InfMat = (1-kap_hes) * (eye(size(InfMat_temp))*1) + kap_hes * InfMat_temp;
        % InfMat_temp = Fdot'*(invFkroninvF)*Fdot;
        % InfMat = (1-kap_hes) * diag(diag(InfMat_temp)) + kap_hes * InfMat_temp;
        Res.InfMat{tt} = InfMat;
        Res.InfMat_temp{tt} = InfMat_temp;

        % Score = InfMat\GradD;
        % try
        %     Score = chol(InfMat)\eye(size(InfMat))*GradD;
        % catch ME
        %     fprintf('eig(Q): %s\n',eig(Q));
        %     fprintf('eig(P): %s\n',eig(P));
        %     fprintf('eig(F): %s\n',eig(F));
        %     fprintf('eig(invF): %s\n',eig(invF));
        %     fprintf('eig(invFkroninvF): %s\n',eig(invFkroninvF));
        %     fprintf('eig(InfMat_temp): %s\n',eig(InfMat_temp));
        %     fprintf('eig(InfMat): %s\n',eig(InfMat));
        %     rethrow(ME);
        % end

        compute_score = svd_pseudoinverse(InfMat)*GradD;
        Score = (1-mom)*compute_score+mom*Score;
        
        Res.score(:,tt) = Score;
        %        Score = pinv(InfMat)*GradD;
        %        Score = (InfMat^(-1))*GradD;
        %        Score = (InfMat^(.5))\GradD;
                
        %        Score = GradD;
        
        
        %% STEP 4: UPDATE PARAMETERS
        if isempty(pred_ft)
            f_t(:,tt+1) = c + A*f_t(:,tt) + B*Score;
    %         f_t(:,tt+1) = f_t(:,tt) + B_coeff*Score;
            Res.f_t(:,tt+1)=f_t(:,tt+1);
        end
        
        %% STEP 5: UPDATE MATRICES T AND Q AND THEIR JACOBIANS GIVEN THE NEW PARAMETERS IN f_t(:,tt+1)
        % 1. Matrix T: defined outside
        % 2. Matrix Q
        [Q,Qdot] = LinkFne_Q_corr1(f_t(:,tt+1));
        Res.Q(:,:,tt+1) =  Q;
    

        Omega = SS*Q*SS';
        Omegadot = kronSS*Qdot;   
            
        %% UPDATE KF
        K = T*G; L = (T-K*Wt*Z);  
        Res.L_t(:,:,tt) = L;

        alfa_t(:,tt+1)          = T*a_t + K*v + SS*mean_err;
        Res.alfa_t(:,tt+1)      = alfa_t(:,tt+1);

        % P_temp                  = T*P*L'+Q;
        % P_temp                  = L*P*L'+K*Wt*H*Wt'*K'+ Omega;
        P_temp                  = L*P*L'+ Omega;
        Res.P_t(:,:,tt+1)       = P_temp;
        P_t(:,tt+1)             = P_temp(:);
    
        % if any(eig(Q)<0)
        %     % Open the log file for appending
        %     logFile = fopen('log.txt', 'a');  % 'a' mode opens for appending
        % 
        %     if logFile == -1
        %         error('Failed to open the log file.');
        %     end
        % 
        %     % Write some entries to the log file
        %     fprintf(logFile, 'Log entry: %s\n', datestr(now));  % Log the current date and time
        %     fprintf(logFile, 'psi: %d\n', UnoMenoUno(Vpar(1)));
        %     fprintf(logFile, 'phi: %d\n', UnoMenoUno(Vpar(2)));
        %     fprintf(logFile, 'means: %d\n', UnoMenoUno(Vpar(3:5)));
        %     fprintf(logFile, 'kap_hes: %.2f\n', MaxZero(.5,Vpar(6)));
        %     fprintf(logFile, 'A: %s\n', UnoMenoUno(Vpar(12:16)));
        %     fprintf(logFile, 'B: %s\n', MaxZero(0.1, Vpar(17:21)));
        %     fprintf(logFile, 'c: %s\n', (eye(5)-diag(A))*initialF1);
        %     fprintf(logFile, '----------\n');
        % 
        %     % Close the log file
        %     fclose(logFile);
        %     LL = -10^10; 
        %     break;
        % end

        if isnan(Lik_t(1,tt))
            LL = -10^10; 
            break;
        end
        
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

% catch ME
%     if strcmp(ME.identifier, 'MATLAB:warn')
%         logWarning(ME.message);
%     else
%         rethrow(ME);
%     end
% end


% function logWarning(message)
%     % Open the log file in append mode
%     logFile = fopen('log.txt', 'a');
    
%     % Write the warning message to the log file
%     fprintf(logFile, '%s: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), message);
    
%     % Close the log file
%     fclose(logFile);
% end

end

