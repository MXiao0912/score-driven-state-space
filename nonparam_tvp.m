function[psi, phi, initialF1] = nonparam_tvp(px_nav)

[PositiveTrans, InvPositiveTrans, DerPositiveTrans, ...
    UnoMenoUno, InvUnoMenoUno, DerUnoMenoUno, ...
    LogisticFun, MaxZero, InvMaxZero, DerMaxZero, vec] = Useful_Transformations();

model_arma_11 = arima('ARLags', 1, 'D', 0, 'MALags', 1);
[estModelpx, estParamspx] = estimate(model_arma_11, removeOutliersFillMedian(px_nav.dprice));
[estModelnav, estParamsnav] = estimate(model_arma_11, removeOutliersFillMedian(px_nav.dnav));
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
kplp = nanmean(px_nav.ptilde.*[NaN;px_nav.ptilde(1:end-1)]);


objfunc = @(vparam) moment_rec_tvp(vparam,[knn, kpp, knp, knln,...
kplp], [psi, phi]);
x0 = [0.1,0.1,0.1,0,0];
lb = [1e-10, 1e-10, 1e-10, -Inf, -Inf];
ub = [Inf, Inf, Inf, Inf, Inf];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = @nonlcon_tvp;

problem = createOptimProblem('fmincon', ...
'objective', objfunc, ...
'x0', x0, ...
'lb', lb, ...
'ub', ub, ...
'nonlcon', @nonlcon_tvp,...
'options', optimoptions('fmincon', 'Algorithm', 'sqp','TolFun',1e-12, 'TolX', 1e-7, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 5000, 'TolCon', 1e-20));
gs = GlobalSearch;
[t] = evalc('[x, fval, flag, out, allmins] = run(gs, problem);');

% Initial covariance matrix
sigma_e = sqrt(x(1));
sigma_w = sqrt(x(2));
sigma_r = sqrt(x(3));


corr_re = x(4)/(sigma_r*sigma_e);
corr_rw = x(5)/(sigma_r*sigma_w);

% case transformation
Init_inv_sigma_r = log(sigma_r);
Init_inv_sigma_e = log(sigma_e);
Init_inv_sigma_w = log(sigma_w);

% cap_vol = 5; 
% Init_inv_sigma_e = InvMaxZero(cap_vol,sigma_e);
% Init_inv_sigma_w = InvMaxZero(cap_vol,sigma_w);
% Init_inv_sigma_r = InvMaxZero(cap_vol,sigma_r);


pi_re = corr_re;
pi_rw = corr_rw/sqrt(1-pi_re^2);

Init_inv_pi_re = atanh(pi_re); 
Init_inv_pi_rw = atanh(pi_rw);

initialF1 = [Init_inv_sigma_e; Init_inv_sigma_w; Init_inv_sigma_r; Init_inv_pi_re; Init_inv_pi_rw];


end