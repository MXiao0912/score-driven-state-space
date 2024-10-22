function [obs, mean_ft, ft] = sim_dgp_corr1(T, psi, phi, mean_ft)

% The observation sample size, state and ft sizes are T+1 as they all need
% initial values. Then we extract the sample from the second period and the
% sample size we use is of length T. 


% % clear all;
% T = 1000; % Number of time steps
% psi = 0.6;
% phi = 0.8;
% % Initialize ft parameters
% mean_ft = [-1; -1; -1; 0.5];

f_ar = 0.5*ones(1,4); % AR(1) coefficient for all parameters
f_b = 0.01*ones(1,4);
f_intercept = (eye(4)-diag(f_ar))*mean_ft;

% Simulate ft process
ft = zeros(4, T+1);
ft(:,1) = mean_ft;
% ft(4,:) = 0.5;
for t = 2:(T+1)
    % ft(:, t) = f_intercept+ diag(f_ar) * (ft(:, t-1)) + diag(f_b)*randn(4, 1);

    ft(1, t) = sin(2*pi*t/T);
    ft(2, t) = sin(2*pi*t/(T/2));
    ft(3, t) = sin(2*pi*(t-200)/T);
    ft(4, t) = sin(4*pi*(t-800)/(T));

    % ft(4,t) = (mean_ft(1:3));
end
%  ft(1, 500:end) = ft(1, 500:end)+1;
%  ft(2, 510:end) = ft(2, 510:end)+1;
%  ft(3, 520:end) = ft(3, 520:end)+1;
%  ft(4, 400:end) = ft(4, 400:end)+0.2;

% ft(4, :) = 0.5;
% ft(5, :) = 0;


fig = figure("Visible","off");
subplot(4, 1, 1);
plot(1:(T+1), exp(ft(1,1:(T+1))), 'LineWidth', 1);
subplot(4, 1, 2);
plot(1:(T+1), exp(ft(2,1:(T+1))), 'LineWidth', 1);
subplot(4, 1, 3);
plot(1:(T+1), exp(ft(3,1:(T+1))), 'LineWidth', 1);
subplot(4, 1, 4);
plot(1:(T+1), tanh(ft(4,1:(T+1))), 'LineWidth', 1);
saveas(fig, 'recovered_true_ft.png');


% Allocate memory for states and observations
obs = table(zeros(T+1, 1), zeros(T+1, 1), zeros(T+1, 1), zeros(T+1, 1), ...
    'VariableNames', {'price', 'nav', 'dprice', 'dnav'});
state = zeros(4, T+1);

% Simulate the state and obs
for t = 2:(T+1)
    sigma_e = exp(ft(1, t));
    sigma_w = exp(ft(2, t));
    sigma_r = exp(ft(3, t));
    pi_ew = tanh(ft(4, t));
    corr_ew = pi_ew;
    corr_re = 0;
    corr_rw = 0;

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




fig = figure("Visible","off");
plot(2:(T+1), obs.price(2:end), 'LineWidth', 1);
hold on
plot(2:(T+1), obs.nav(2:end), 'LineWidth', 1);
hold on
plot(2:(T+1), state(1,2:end), 'LineWidth', 1);
legend('px', 'nav','state');
hold off
saveas(fig, 'true_price_states.png');


obs = obs(2:end, :);


end
