clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate dgp

T = 1000; % Number of time steps
psi = 0.6;
phi = 0.8;
% Initialize ft parameters
mean_ft = [-1; -1; -1; 0.5];

% generate dgp
[obs, mean_ft, ft] = sim_dgp_corr1(T, psi, phi, mean_ft);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate parameters
run Useful_Transformations

% user-specify A and B
% [psi, phi, initF1] = nonparam_tvp(obs);
A = InvMaxZero(0.99, 0.5*ones(4,1)); 
B = InvMaxZero(0.1, [0.015*ones(4,1)]);
c = (eye(4)-diag(MaxZero(0.99, A)))*mean_ft;
kap_hes = InvMaxZero(1,.015);
mom_ = InvMaxZero(1, 0.5);

% InitialParams = [kap_hes;c;A;B;mom_];
% fixedParams = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1); NaN(1); NaN(4,1); NaN(4,1); NaN(4,1); NaN(1)];


InitialParams = [InvUnoMenoUno(0.5);InvUnoMenoUno(0.5); zeros(3,1); kap_hes;c;A;B;mom_];
% fixedParams = [NaN(5,1); [kap_hes;c;A;B;mom_]];
fixedParams  = NaN(19,1);


y = table2array(obs(:,{'price','nav'}))';

[EstimParams, FullParams] = sim_est_corr1(y, mean_ft, "patternsearch", false, InitialParams, fixedParams, []);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examine estimated results

psi = UnoMenoUno(FullParams(1))
phi = UnoMenoUno(FullParams(2))
err_mean = FullParams(3:5)

kap_hes = MaxZero(1,FullParams(6))
c = FullParams(7:10)
A = MaxZero(0.99, FullParams(11:14))
B = MaxZero(0.01, FullParams(15:18))
mom = MaxZero(1,FullParams(19))

[LL,Res]=kf_tvp_corr1(FullParams,y, mean_ft,[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% true parameter likelihood
[LL_T,Res_T]=kf_tvp_corr1(insertParams([InvUnoMenoUno(0.6);InvUnoMenoUno(0.8); zeros(3,1); kap_hes;c;A;B;mom_], fixedParams),y, mean_ft, ft(:,2:end));

% compare est ft with true param ft
fig = figure();
subplot(4,1,1);
plot(1:(T-1), squeeze(Res.f_t(1,2:end)));
hold on;
plot(1:(T-1), squeeze(Res_T.f_t(1,2:end)));
subplot(4,1,2);
plot(1:(T-1), squeeze(Res.f_t(2,2:end)));
hold on;
plot(1:(T-1), squeeze(Res_T.f_t(2,2:end)));
subplot(4,1,3);
plot(1:(T-1), squeeze(Res.f_t(3,2:end)));
hold on;
plot(1:(T-1), squeeze(Res_T.f_t(3,2:end)));
subplot(4,1,4);
plot(1:(T-1), squeeze(Res.f_t(4,2:end)));
hold on;
plot(1:(T-1), squeeze(Res_T.f_t(4,2:end)));
saveas(fig, strcat('Diagnostics/trail_ft','.png'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check se of ft

% compare est ft with true param ft
color = "g";
alpha=0.3;
fig = figure();
subplot(4,1,1);
plot(1:(T-1), squeeze(Res.Q(1,1,2:end)));
hold on;
confplot(1:(T-1), squeeze(Res.Q(1,1,2:end)),squeeze(mad(Save_Collect.vare,0,2) * 1.96)', color, alpha);
hold on;
plot(1:(T-1), exp(ft(1,3:end)).^2);
subplot(4,1,2);
plot(1:(T-1), squeeze(Res.Q(2,2,2:end)));
hold on;
confplot(1:(T-1), squeeze(Res.Q(2,2,2:end)),squeeze(mad(Save_Collect.varw,0,2) * 1.96)', color, alpha);
hold on;
plot(1:(T-1), exp(ft(2,3:end)).^2);
subplot(4,1,3);
plot(1:(T-1), squeeze(Res.Q(3,3,2:end)));
hold on;
confplot(1:(T-1), squeeze(Res.Q(3,3,2:end)),squeeze(mad(Save_Collect.varr,0,2) * 1.96)', color, alpha);
hold on;
plot(1:(T-1), exp(ft(3,3:end)).^2);
subplot(4,1,4);
plot(1:(T-1), squeeze(Res.Q(2,1,2:end)./(Res.Q(1,1,2:end).*Res.Q(2,2,2:end)).^(0.5)));
hold on;
confplot(1:(T-1), squeeze(Res.Q(2,1,2:end)./(Res.Q(1,1,2:end).*Res.Q(2,2,2:end)).^(0.5)), squeeze(mad(Save_Collect.correw,0,2) * 1.96)', color, alpha);
hold on;
plot(1:(T-1), tanh(ft(4,3:end)));
saveas(fig, strcat('Diagnostics/trail_ft','.png'));

% patternsearch:    -3.4159e+03
% csminwel:  -4.0251e+03
% fminunc: 4.07385 

% LL_T =

%   -3.3346e+03
%


fig = figure();
histogram(Save_Collect.vare(200,:));
saveas(fig, strcat('Diagnostics/simu_var','.png'));