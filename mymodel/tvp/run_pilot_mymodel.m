function run_pilot_mymodel(j)
    j = str2double(j);
    % j=12;
    cd("/rds/user/mx235/hpc-work/px_nav_ssm/mytvp/mymodel/tvp/")
    addpath("../")
    run Useful_Transformations.m

    px_nav_tot = readtable("../../../mydata/px_nav.csv");
    isin_list = unique(px_nav_tot.isin);

    subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

    for i = 1:length(isin_list)
        subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
    end

    px_nav = subTables{j};
    px_nav.lprice = log(px_nav.price);
    px_nav.lnav = log(px_nav.nav);
    y = table2array(px_nav(:,{'lprice','lnav'}))';
    save(strcat('yall/',px_nav.isin{j}), 'y');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial parameter values
    params = load(strcat('../param_est/',px_nav.isin{j})).param_tab;
    psi_u = params.Value(strcmp(params.Names, "psi_u"));
    psi_p = params.Value(strcmp(params.Names, "psi_p"));
    psi_n = params.Value(strcmp(params.Names, "psi_n"));
    psi_T = params.Value(strcmp(params.Names, "psi_T"));
    theta_u = params.Value(strcmp(params.Names, "theta_u"));
    theta_p = params.Value(strcmp(params.Names, "theta_p"));

    eu = params.Value(strcmp(params.Names, "eu"));
    ep = params.Value(strcmp(params.Names, "ep"));
    en = params.Value(strcmp(params.Names, "en"));
    pu = params.Value(strcmp(params.Names, "pu"));
    un = params.Value(strcmp(params.Names, "un"));
    pn = params.Value(strcmp(params.Names, "pn"));

    cov_half_params = [params.Value(7),0,0;params.Value(10),params.Value(8),0;params.Value(11),params.Value(12),params.Value(9)];
    cov_init = cov_half_params * (cov_half_params');
    corr_pu = cov_init(2,1)/sqrt(cov_init(1,1)*cov_init(2,2));
    corr_nu = cov_init(3,1)/sqrt(cov_init(1,1)*cov_init(3,3));
    corr_pn = cov_init(3,2)/sqrt(cov_init(2,2)*cov_init(3,3));
    pi_pu = corr_pu;
    pi_nu = corr_nu;
    pi_pn = (corr_pn-pi_pu*pi_nu)/sqrt((1-pi_pu^2)*(1-pi_nu^2));

    mean_ft = [max(min(log(sqrt(cov_init(1,1))),0.5),-10);...
                max(min(log(sqrt(cov_init(2,2))),0.5),-10);...
                max(min(log(sqrt(cov_init(3,3))),0.5),-10);...
                max(min(atanh(pi_pu),3),-3);...
                max(min(atanh(pi_nu),3),-3);...
                max(min(atanh(pi_pn),3),-3)];

    A = InvMaxZero(0.99, 0.5*ones(6,1)); 
    B = InvMaxZero(0.1, [0.015*ones(6,1)]);
    c = (eye(6)-diag(MaxZero(0.99, A)))*mean_ft;
    kap_hes = InvMaxZero(1,.015);
    mom_ = InvMaxZero(1, 0.5);

    InitialParams = [InvMaxZero(1,psi_u);InvMaxZero(1,psi_p);InvMaxZero(1,psi_n);InvMaxZero(1,psi_T);InvMaxZero(-1,theta_u);InvMaxZero(-1,theta_p);kap_hes;c;A;B;mom_];

    save(strcat('mean_ft/',px_nav.isin{j}), 'mean_ft');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimation
    addpath('../../../Sims_Solver/')
    LossToMinmize = @(vparam) -kf_tvp_mymodel(vparam,y, mean_ft);
    lb = [(-5)*ones(1,6), -Inf, (-2)*ones(1,6), (-5)*ones(1,6), (-3)*ones(1,6), -Inf];
    ub = [(5)*ones(1,6), Inf, (2)*ones(1,6), (5)*ones(1,6), (3)*ones(1,6), Inf];

    options = optimoptions('patternsearch', ...
            'MaxIterations', 1000, ... 
            'UseParallel', true, ...
            'Display', 'iter',...
            'ConstraintTolerance', 1e-6,...
            'MeshTolerance', 1e-10,...
            'StepTolerance', 1e-10);
    [EstimParams, fval] = patternsearch(LossToMinmize, InitialParams, [], [], [], [], lb, ub, @nonlcon, options);
    save(strcat('rawestimates/',px_nav.isin{j}), 'EstimParams');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Examine estimated parameters
    names = ["psi_u","psi_p","psi_n","psi_T","theta_u","theta_p","kap_hes",repmat("c",1,6),repmat("A",1,6),repmat("B",1,6),"mom"]; % lower half of cov matrix
    OutParams = [MaxZero(1,EstimParams(1)),MaxZero(1,EstimParams(2)),MaxZero(1,EstimParams(3)),MaxZero(1,EstimParams(4)),...
    MaxZero(-1,EstimParams(5)),MaxZero(-1,EstimParams(6)),MaxZero(1,EstimParams(7)),EstimParams(8:13)',...
    MaxZero(0.99, EstimParams(14:19)'),MaxZero(0.1, EstimParams(20:25)'),MaxZero(1, EstimParams(26))];
    param_tab = table(names', OutParams', 'VariableNames', {'Names', 'Value'});
    save(strcat('param_est_score/',px_nav.isin{j}), 'param_tab');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Examine ft
    [LL,Res]=kf_tvp_mymodel(EstimParams,y, mean_ft);
    save(strcat('Resall/',px_nav.isin{j}), 'Res');

    fig = figure();
    subplot(6,1,1);
    plot(px_nav.date, squeeze(Res.f_t(1,2:end)));
    subplot(6,1,2);
    plot(px_nav.date, squeeze(Res.f_t(2,2:end)));
    subplot(6,1,3);
    plot(px_nav.date, squeeze(Res.f_t(3,2:end)));
    subplot(6,1,4);
    plot(px_nav.date, squeeze(Res.f_t(4,2:end)));
    subplot(6,1,5);
    plot(px_nav.date, squeeze(Res.f_t(5,2:end)));
    subplot(6,1,6);
    plot(px_nav.date, squeeze(Res.f_t(6,2:end)));
    saveas(fig, strcat('ft/',px_nav.isin{j},'.png'));



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Examine Qt and its se
    color = "g";
    alpha=0.3;
    fig = figure();
    subplot(6,1,1);
    plot(px_nav.date(5:end), squeeze(Res.Q(1,1,6:end)));
    % hold on;
    % confplot(1:(T-1), squeeze(Res.Q(1,1,2:end)),squeeze(mad(Save_Collect.vare,0,2) * 1.96)', color, alpha);
    % hold on;
    % plot(1:(T-1), exp(ft(1,3:end)).^2);
    subplot(6,1,2);
    plot(px_nav.date(5:end), squeeze(Res.Q(2,2,6:end)));
    % hold on;
    % confplot(1:(T-1), squeeze(Res.Q(2,2,2:end)),squeeze(mad(Save_Collect.varw,0,2) * 1.96)', color, alpha);
    % hold on;
    % plot(1:(T-1), exp(ft(2,3:end)).^2);
    subplot(6,1,3);
    plot(px_nav.date(5:end), squeeze(Res.Q(3,3,6:end)));
    % hold on;
    % confplot(1:(T-1), squeeze(Res.Q(3,3,2:end)),squeeze(mad(Save_Collect.varr,0,2) * 1.96)', color, alpha);
    % hold on;
    % plot(1:(T-1), exp(ft(3,3:end)).^2);
    subplot(6,1,4);
    plot(px_nav.date(5:end), squeeze(Res.Q(2,1,6:end)./(Res.Q(1,1,6:end).*Res.Q(2,2,6:end)).^(0.5)));
    % hold on;
    % confplot(1:(T-1), squeeze(Res.Q(2,1,2:end)./(Res.Q(1,1,2:end).*Res.Q(2,2,2:end)).^(0.5)), squeeze(mad(Save_Collect.correw,0,2) * 1.96)', color, alpha);
    % hold on;
    % plot(1:(T-1), tanh(ft(4,3:end)));
    subplot(6,1,5);
    plot(px_nav.date(5:end), squeeze(Res.Q(3,1,6:end)./(Res.Q(1,1,6:end).*Res.Q(3,3,6:end)).^(0.5)));
    subplot(6,1,6);
    plot(px_nav.date(5:end), squeeze(Res.Q(3,2,6:end)./(Res.Q(2,2,6:end).*Res.Q(3,3,6:end)).^(0.5)));
    saveas(fig, strcat('Qt/',px_nav.isin{j},'.png'));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filtered states
    color = "g";
    alpha=0.3;
    fig = figure('Visible', 'off');
    set(fig, 'Position', [100, 100, 1600, 400]);
    plot(px_nav.date(2:end)', y(:,2:end), 'LineWidth', 1);
    hold on
    confplot(px_nav.date(2:end)', Res.alfa_t(1,3:end),squeeze(sqrt(Res.P_t(1,1,3:end)) * 1.96)', color, alpha);
    legend('px', 'nav','state');
    hold off
    saveas(fig, strcat('px_nav_state_plot/',px_nav.isin{j},'.png'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagnostic checks
    residual = cell(length(Res.v_t));
    for i=1:length(Res.v_t)
        residual{i} = chol(Res.invF_t{i})*Res.v_t{i};
    end

    for tt=1:length(residual)
        [Wt,nt]=SelectMatW(y(:,tt));
        if isequal(Wt,[1 0])
            residual{tt} = [residual{tt}; NaN];
        elseif isequal(Wt,[0 1])
            residual{tt} = [NaN; residual{tt}];
        elseif isequal(Wt,eye(2))
            residual{tt} = residual{tt};
        else
            residual{tt} = [NaN; NaN];
        end
    end

    % Diagnostics Check
    resid = [residual{:}];
    [p_jb, p_arch, p_lb]=gen_diagnostics(px_nav, resid);
end