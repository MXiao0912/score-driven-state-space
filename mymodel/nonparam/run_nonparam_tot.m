function run_nonparam_tot(j)
    j = str2double(j);
    cd("/rds/user/mx235/hpc-work/px_nav_ssm/mytvp/mymodel/nonparam")
    px_nav_tot = readtable("../../../mydata/px_nav.csv");
    isin_list = unique(px_nav_tot.isin);

    subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

    for i = 1:length(isin_list)
        subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
    end

    px_nav = subTables{j};
    px_nav.lprice = log(px_nav.price);
    px_nav.lnav = log(px_nav.nav);
    px_nav.dprice = [NaN; diff(px_nav.lprice)];
    px_nav.dnav = [NaN; diff(px_nav.lnav)];
    px_nav.dprice = fillmissing(px_nav.dprice,'previous');
    px_nav.dnav = fillmissing(px_nav.dnav,'previous');

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

    ptilde = px_nav.dprice(4:end)-(psi_u+psi_p)*px_nav.dprice(3:end-1)+psi_u*psi_p*px_nav.dprice(2:end-2);
    ntilde = px_nav.dnav(4:end)+(-psi_n-psi_T-psi_u)*px_nav.dnav(3:end-1)+(psi_n*psi_T + psi_n*psi_u + psi_T*psi_u)*px_nav.dnav(2:end-2)-psi_n*psi_T*psi_u*px_nav.dnav(1:end-3);

    c1 = 2*(1 + psi_u * (-1 + theta_p)^2 - theta_p + theta_p^2 + psi_u^2 * (1 - theta_p + theta_p^2));
    c2 = 1 - 2 * psi_p * theta_u + theta_u^2 + psi_p^2 * (1 + theta_u^2);
    c3 = 2*(1 + (-1 - psi_u + theta_p) * theta_u + psi_p * (1 + theta_p * (-1 + theta_u) + psi_u * (1 + (-1 + theta_p) * theta_u)));
    c4 = 2 - psi_u^2 * (-2 + theta_p) - 2 * psi_u * (-1 + theta_p) - theta_p;
    c5 = 1 + psi_p * (1 - psi_u * (-1 + theta_u)) - (1 + psi_u) * theta_u;
    c6 = (-1 + psi_n) * (1 + (-1 - psi_u + theta_p) * theta_u + psi_T * (1 + theta_p * (-1 + theta_u) + psi_u * (1 + (-1 + theta_p) * theta_u)));
    c7 = (-1 + psi_n) * (1 - psi_T * theta_u + theta_u^2 + psi_p * (-theta_u + psi_T * (1 + theta_u^2)));
    c8 = 2*(1 + psi_u + psi_u^2);
    c9 = 2*((-1 + psi_n) * (-1 + psi_T * (-1 + psi_u * (-1 + theta_u)) + (1 + psi_u) * theta_u));
    c10 = (-1 + psi_n)^2 * (1 - 2 * psi_T * theta_u + theta_u^2 + psi_T^2 * (1 + theta_u^2));
    c11 = -(-1 + theta_p)^2 - psi_u^2 * (-1 + theta_p)^2 + psi_u * (-2 + 3 * theta_p - 2 * theta_p^2);
    c12 = -((-1 + theta_p) * (-1 + theta_u)) - psi_p * (1 + psi_u) * (-1 + theta_p) * (-1 + theta_u) + psi_u * (-1 - (-1 + theta_p) * theta_u);
    c13 = theta_u + psi_p^2 * theta_u - psi_p * (1 + theta_u^2);
    c14 = -1 + 2 * psi_u * (-1 + theta_p) + 2 * theta_p + psi_u^2 * (-1 + 2 * theta_p);
    c15 = theta_u + psi_p * (-1 + (1 + psi_u) * theta_u);
    c16 = (-1 + psi_n) * (1 + theta_p * (-1 - psi_T + theta_u) + psi_u * (1 + psi_T * (1 + theta_p * (-1 + theta_u)) + (-1 + theta_p) * theta_u));
    c17 = (-1 + psi_n) * (-theta_u + psi_p * (1 - psi_T * theta_u + theta_u^2));
    c18 = -1 - psi_u^2 + psi_u * (-2 + theta_p);
    c19 = -1 + psi_u * (-1 - psi_p + theta_u);
    c20 = (-1 + psi_n) * (theta_u + psi_T * (-1 + (1 + psi_u - theta_p) * theta_u));
    c21 = (-1 + psi_n) * (-theta_u + psi_T * (1 - psi_p * theta_u + theta_u^2));
    c22 = -(1+psi_u)^2;
    c23 = -((-1 + psi_n) * (1 + psi_T) * (1 + psi_u) * (-1 + theta_u));
    c24 = (-1 + psi_n)^2 * (theta_u + psi_T^2 * theta_u - psi_T * (1 + theta_u^2));

    t_ = (1/size(ptilde,1)):(1/size(ptilde,1)):1;
    silverman = @(n) 1.06*min(std(t_),iqr(t_)/1.34)/(n^0.2);
    h = silverman(size(ptilde,1))*size(ptilde,1);

    knn = smoothdata(ntilde.^2, 'gaussian',h);
    knp = smoothdata(ntilde.*ptilde, 'gaussian',h);
    kpp = smoothdata(ptilde.^2, 'gaussian',h);
    knln = smoothdata(ntilde.*[NaN;ntilde(1:end-1)], 'gaussian',h);
    knlp = smoothdata(ntilde.*[NaN;ptilde(1:end-1)], 'gaussian',h);
    kpln = smoothdata(ptilde.*[NaN;ntilde(1:end-1)], 'gaussian',h);
    kplp = smoothdata(ptilde.*[NaN;ptilde(1:end-1)], 'gaussian',h);


    res = cell(1,size(knln,1));
    res_fval = cell(1,size(knln,1)); 
    res_flag = cell(1,size(knln,1));
    res_out = cell(1,size(knln,1));
    res_allmins = cell(1,size(knln,1));

    for k=2:size(knln,1)
        disp(k)
        objfunc = @(vparam) moment_rec(vparam,[knn(k,1), kpp(k,1), knp(k,1), knln(k,1),...
            kplp(k,1), knlp(k,1), kpln(k,1)], [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24]);
        x0 = [cov_init(1,1), cov_init(2,2),cov_init(3,3),cov_init(2,1),cov_init(3,1),cov_init(3,2)];
        lb = [1e-10, 1e-10, 1e-10, -Inf, -Inf, -Inf];
        ub = [Inf, Inf, Inf, Inf, Inf, Inf];
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        %[x, fval, exitflag, output]=fmincon(objfunc, x0, A, b, Aeq, beq, lb, ub, nonlcon, optimoptions('fmincon','Algorithm', 'sqp','TolFun',1e-20, 'TolCon',1e-20));
        % [x, fval, exitflag, output]=fmincon(objfunc, x0, A, b, Aeq, beq, lb, ub, @nonlcon, optimoptions('fmincon', 'Algorithm', 'sqp','TolFun',1e-12, 'TolX', 1e-7, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 5000, 'TolCon', 1e-10));

        % [x, fval, exitflag, output]=fmincon(objfunc, x0, A, b, Aeq, beq, lb, ub, @nonlcon, optimoptions('fmincon', 'Algorithm', 'interior-point', 'HessianApproximation', 'lbfgs','TolFun',1e-12, 'TolX', 1e-7, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 5000, 'TolCon', 1e-10));

        problem = createOptimProblem('fmincon', ...
            'objective', objfunc, ...
            'x0', x0, ...
            'lb', lb, ...
            'ub', ub, ...
            'nonlcon', @nonlcon,...
            'options', optimoptions('fmincon', 'Algorithm', 'sqp','TolFun',1e-12, 'TolX', 1e-7, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 5000, 'TolCon', 1e-10));
        %ms = MultiStart('Display', 'iter', 'UseParallel', true, 'StartPointsToRun','bounds-ineqs');
        
        % Run MultiStart with fmincon
        %[x1, x2, x3, x4] = meshgrid(linspace(0.1,1,10), linspace(0.1,1,10), linspace(0.1,1,10), linspace(-0.9,0.9,20));
        %start_points = [x1(:),x2(:),x3(:),x4(:)];
        %start_points_set = CustomStartPointSet(start_points);
        %[x, fval, flag, out, allmins] = run(ms, problem, 100); % 50 runs from different start points
        gs = GlobalSearch;
        %[x, fval, flag, out, allmins] = run(gs, problem);
        [t] = evalc('[x, fval, flag, out, allmins] = run(gs, problem);');

        res{k}=x;
        res_fval{k}=fval;
        % res_flag{k}=flag;
        % res_out{k}=out;
        % res_allmins{k}=allmins;
    end

    M = cell2mat(res);
    M = reshape(M,6,[]).';
    % M(:,1) = smoothdata(M(:,1), 'gaussian',h);
    % M(:,2) = smoothdata(M(:,2), 'gaussian',h);
    % M(:,3) = smoothdata(M(:,3), 'gaussian',h);
    % M(:,4) = smoothdata(M(:,4), 'gaussian',h);
    % M(:,5) = smoothdata(M(:,5), 'gaussian',h);
    % M(:,6) = smoothdata(M(:,6), 'gaussian',h);

    M(:,4) = M(:,4)./sqrt(M(:,1).*M(:,2));
    M(:,5) = M(:,5)./sqrt(M(:,1).*M(:,3));
    M(:,6) = M(:,6)./sqrt(M(:,2).*M(:,3));

    corr_values = M(:,4);
    corr_values(corr_values > 0.9 | corr_values < -0.9) = NaN;
    M(:,4) = corr_values;

    corr_values = M(:,5);
    corr_values(corr_values > 0.9 | corr_values < -0.9) = NaN;
    M(:,5) = corr_values;

    corr_values = M(:,6);
    corr_values(corr_values > 0.9 | corr_values < -0.9) = NaN;
    M(:,6) = corr_values;

    figure;  % Open a new figure window
    subplot(2,3,1);
    %plot(px_nav.date(2:end), abs(M(:,4)), 'o-', MarkerSize=1);  % Change 'o-' to any other plot style as needed
    plot(px_nav.date(5:end), M(:,1), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('\sigma_\mu^2');  % Label y-axis

    subplot(2,3,2);
    plot(px_nav.date(5:end), M(:,2), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('\sigma_p^2');  % Label y-axis

    subplot(2,3,3);
    plot(px_nav.date(5:end), M(:,3), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('\sigma_n^2');  % Label y-axis

    subplot(2,3,4);
    plot(px_nav.date(5:end), M(:,4), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('\rho_{p\mu}');  % Label y-axis

    subplot(2,3,5);
    plot(px_nav.date(5:end), M(:,5), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('\rho_{n\mu}');  % Label y-axis

    subplot(2,3,6);
    plot(px_nav.date(5:end), M(:,6), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('\rho_{pn}');  % Label y-axis

    saveas(gcf, strcat(px_nav.isin{1},'.png'))

    save(strcat('time_varying_cov/',px_nav.isin{j}), 'M');

end

% 109 and 53 are ridiculously slow, check why