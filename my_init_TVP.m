v2struct(SettingsTransfor)

px_nav = subTables{j}(1:30,:);
px_nav.lprice = log(px_nav.price);
px_nav.lnav = log(px_nav.nav);
px_nav.dprice = [NaN; diff(px_nav.lprice)];
px_nav.dnav = [NaN; diff(px_nav.lnav)];
px_nav.dprice = fillmissing(px_nav.dprice,'previous');
px_nav.dnav = fillmissing(px_nav.dnav,'previous');

model_arma_11 = arima('ARLags', 1, 'D', 0, 'MALags', 1);
[estModelpx, estParamspx] = estimate(model_arma_11, px_nav.dprice);
[estModelnav, estParamsnav] = estimate(model_arma_11, px_nav.dnav);
psi = estModelpx.AR{1};
phi = estModelnav.AR{1};

px_nav.ldp = [NaN; px_nav.dprice(1:end-1)];
px_nav.ldn = [NaN; px_nav.dnav(1:end-1)];
px_nav.ptilde = px_nav.dprice-psi*px_nav.ldp;
px_nav.ntilde = px_nav.dnav-phi*px_nav.ldn;

knn = nanmean(px_nav.ntilde.^2);
knp = nanmean(px_nav.ntilde.*px_nav.ptilde);
kpp = nanmean(px_nav.ptilde.^2);
knln = nanmean(px_nav.ntilde.*[NaN;px_nav.ntilde(1:end-1)]);
knlp = nanmean(px_nav.ntilde.*[NaN;px_nav.ptilde(1:end-1)]);
kpln = nanmean(px_nav.ptilde.*[NaN;px_nav.ntilde(1:end-1)]);
kplp = nanmean(px_nav.ptilde.*[NaN;px_nav.ptilde(1:end-1)]);


objfunc = @(vparam) moment_rec(vparam,[knn, kpp, knp, knln,...
kplp, knlp, kpln], [psi, phi]);
x0 = [0,0,0,0,0,0];
lb = [1e-10, 1e-10, 1e-10, -Inf, -Inf, -Inf];
ub = [Inf, Inf, Inf, Inf, Inf, Inf];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = @nonlcon;
problem = createOptimProblem('fmincon', ...
'objective', objfunc, ...
'x0', x0, ...
'lb', lb, ...
'ub', ub, ...
'nonlcon', @nonlcon,...
'options', optimoptions('fmincon', 'Algorithm', 'sqp','TolFun',1e-12, 'TolX', 1e-7, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 5000, 'TolCon', 1e-20));
gs = GlobalSearch;
[t] = evalc('[x, fval, flag, out, allmins] = run(gs, problem);');

cov_init = [x(1) x(4) x(5); x(4) x(2) x(6); x(5) x(6) x(3)];
HL_init = chol(cov_init)';

% Init TVP
y = table2array(px_nav(:,{'lprice','lnav'}))';
InitialParams = [InvUnoMenoUno(psi),InvUnoMenoUno(phi),HL_init(1,1), HL_init(2,2), HL_init(3,3), HL_init(2,1), HL_init(3,1), HL_init(3,2)];
LossToMinmize = @(vparam) -base_kf_missing(vparam,y)/(size(y,2)-1);
optionsIVAN = optimset('Display', 'iter-detailed','LargeScale', 'off','MaxFunEvals',5000);
[EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);

% Initial covariance matrix
init_cov_half = [EstimParams(3), 0, 0; EstimParams(6), EstimParams(4), 0; EstimParams(7), EstimParams(8), EstimParams(5)]
init_cov = init_cov_half*(init_cov_half')

sigma_r = sqrt(init_cov(1,1));
sigma_e = sqrt(init_cov(2,2));
sigma_w = sqrt(init_cov(3,3));
corr_re = init_cov(2,1)/(sigma_r*sigma_e);
corr_rw = init_cov(3,1)/(sigma_r*sigma_w);

% case transformation
Init_inv_sigma_r = log(sigma_r);
Init_inv_sigma_e = log(sigma_e);
Init_inv_sigma_w = log(sigma_w);

pi_re = corr_re;
pi_rw = corr_rw;
pi_ew = -pi_re*pi_rw/sqrt((1-pi_re^2)*(1-pi_rw^2));

Init_inv_pi_re = atanh(pi_re); 
Init_inv_pi_rw = atanh(pi_rw);
Init_inv_pi_ew = atanh(pi_ew); 

initialF1 = [Init_inv_sigma_r; Init_inv_sigma_e; Init_inv_sigma_w; Init_inv_pi_re; Init_inv_pi_rw; Init_inv_pi_ew ];

