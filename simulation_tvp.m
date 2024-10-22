% Time-Varying State-Space Model Parameters
T = 10000; % Number of time steps
psi = 0.6;
phi = 0.8;


% Initialize ft parameters
mean_ft = [-2; -2; -2; 0.1; 0.1];

f_ar = [0.9,0.9,0.9,0.3,0.3]; % AR(1) coefficient for all parameters
% f_ar = zeros(1, 5);
f_b = 0.3;

f_intercept = (eye(5)-diag(f_ar))*mean_ft;

% Simulate ft process
ft = zeros(5, T);
ft(:,1) = f_intercept;
for t = 2:T
    ft(:, t) = f_intercept+ diag(f_ar) * (ft(:, t-1)) + f_b*randn(5, 1);
end
% ft(1, 500:510) = ft(1, 500:510)+2;
% ft(2, 500:510) = ft(2, 500:510)+2;
% ft(3, 500:510) = ft(3, 500:510)+2;
% ft(4, 500:510) = ft(4, 500:510)+1;
% ft(5, 500:510) = ft(5, 500:510)+2;

fig = figure();
subplot(5, 1, 1);
plot(1:T, exp(ft(1,1:T)), 'LineWidth', 1);
subplot(5, 1, 2);
plot(1:T, exp(ft(2,1:T)), 'LineWidth', 1);
subplot(5, 1, 3);
plot(1:T, exp(ft(3,1:T)), 'LineWidth', 1);
subplot(5, 1, 4);
plot(1:T, tanh(ft(4,1:T)), 'LineWidth', 1);
subplot(5, 1, 5);
plot(1:T, tanh(ft(5,1:T)).*sqrt(1-tanh(ft(4,1:T)).^2), 'LineWidth', 1);


% Allocate memory for states and observations
obs = table(zeros(T, 1), zeros(T, 1), zeros(T, 1), zeros(T, 1), ...
    'VariableNames', {'price', 'nav', 'dprice', 'dnav'});
state = zeros(4, T);

% Simulate the system
for t = 2:T
    sigma_e = exp(ft(1, t));
    sigma_w = exp(ft(2, t));
    sigma_r = exp(ft(3, t));
    pi_re = tanh(ft(4, t));
    pi_rw = tanh(ft(5, t));
    corr_re = pi_re;
    corr_rw = pi_rw*sqrt(1-pi_re^2);

    corr_ew = 0; % Assume correlation between epsilon and omega is zero

    % Construct the covariance matrix
    covMatrix = [sigma_e^2, corr_ew * sigma_e * sigma_w, corr_re * sigma_e * sigma_r; 
                 corr_ew * sigma_e * sigma_w, sigma_w^2, corr_rw * sigma_w * sigma_r; 
                 corr_re * sigma_e * sigma_r, corr_rw * sigma_w * sigma_r, sigma_r^2];

    L = chol(covMatrix, 'lower');  % Cholesky decomposition
    noise = L * randn(3, 1);
    epsilon = noise(1);
    omega = noise(2);
    r = noise(3);

    % state transition
    state(:, t) = [1, 0, 0, 0; 1, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0] * state(:, t-1) + ...
                   [0, 0, 1; 0, 0, 0; 1, 0, 0; 0, 1, 0] * [epsilon; omega; r];

    % obs transition
    price_nav = [psi, 0; 0, phi] * [obs.price(t-1); obs.nav(t-1)] + ...
                [1, -psi, 1, 0; 1-phi, 0, 0, 1] * [state(1,t); state(2, t); epsilon; omega];
    
    % Calculate first differences
    dprice = price_nav(1) - obs.price(t-1);
    dnav = price_nav(2) - obs.nav(t-1);
    
    % Update states
    obs.price(t) = price_nav(1);
    obs.nav(t) = price_nav(2);
    obs.dprice(t) = dprice;
    obs.dnav(t) = dnav;
    
end


fig = figure();
plot(1:(T-1), obs.price(2:end), 'LineWidth', 1);
hold on
plot(1:(T-1), obs.nav(2:end), 'LineWidth', 1);
hold on
plot(1:(T-1), state(1,2:end), 'LineWidth', 1);
legend('px', 'nav','state');
hold off


[psi, phi, initF1] = nonparam_tvp(obs);

run Useful_Transformations.m

% user-specify A and B
A = InvMaxZero(1,.9*ones(5,1)); 
B = InvMaxZero(0.1, 0.015*ones(5,1));
kap_hes = InvMaxZero(.5,.015);
c = (eye(5)-diag(MaxZero(1,A)))*initF1;
InitialParams = [InvUnoMenoUno(psi);InvUnoMenoUno(phi);0;0;0;kap_hes;c;A;B];
y = table2array(obs(2:end,{'price','nav'}))';
 
% TVP
LossToMinmize = @(vparam) -kf_tvp(vparam,y, initF1)/(size(y,2)-1);
optionsIVAN = optimoptions('fminunc', 'Display', 'iter-detailed','MaxFunEvals',5000, 'StepTolerance', 1e-10);
[EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);
% [fhat,EstimParams] = csminwel(LossToMinmize,InitialParams,eye(size(InitialParams,1))*.5,[] ,1e-14,100);

psi = UnoMenoUno(EstimParams(1))
phi = UnoMenoUno(EstimParams(2))
kap_hes = MaxZero(.5,EstimParams(6));
A = MaxZero(1,EstimParams(11:15))
B = MaxZero(0.1, EstimParams(16:20))
c = (eye(5)-diag(A))*initF1



% Plot Results
figure;
subplot(2,1,1);
plot(1:T, states(1, :), 'r', 1:T, x_est(1, :), 'b--');
title('State p_t');
legend('True State', 'Estimated State');

subplot(2,1,2);
plot(1:T, states(2, :), 'r', 1:T, x_est(2, :), 'b--');
title('State n_t');
legend('True State', 'Estimated State');