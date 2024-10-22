clear all;
cd("mymodel")
run Useful_Transformations.m
addpath("../")
%IE00BYZK4776

px_nav_tot = readtable("../../mydata/px_nav.csv");
isin_list = unique(px_nav_tot.isin);

subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

for i = 1:length(isin_list)
    subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
end

parfor j = 1:length(isin_list)
    % Non-parametric estimation of F1:  and initial parameter values
    px_nav = subTables{j};
    px_nav.lprice = log(px_nav.price);
    px_nav.lnav = log(px_nav.nav);
    px_nav.dprice = [NaN; diff(px_nav.lprice)];
    px_nav.dnav = [NaN; diff(px_nav.lnav)];
    px_nav.dprice = fillmissing(px_nav.dprice,'previous');
    px_nav.dnav = fillmissing(px_nav.dnav,'previous');

    initprice = fillmissing(px_nav.lprice,'previous');
    initnav = fillmissing(px_nav.lnav,'previous');
    initfund = (initprice+initnav)/2;
    dinitfund = diff(initfund);

    Mdl = arima('ARLags', 1, 'MALags', 1);
    EstMdl = estimate(Mdl, dinitfund);
    psi_u = min(max(1e-6,EstMdl.AR{1}),1-(1e-6));
    theta_u = min(max(-1+(1e-6),EstMdl.MA{1}),-(1e-6));
    eu = EstMdl.Variance;

    EstMdl = estimate(Mdl, initprice-initfund);
    psi_p = min(max(1e-6,EstMdl.AR{1}),1-(1e-6));
    theta_p = min(max(-1+1e-6,EstMdl.MA{1}),-(1e-6));
    ep = EstMdl.Variance;

    model = armax(misdata(iddata(diff(initnav),dinitfund)), [2 2 0 0]);
    psi_n = min(max(1e-6,1-model.B(1)),1-(1e-6));
    if psi_n==1
        psi_T = min(max(1e-6,model.A(3)),1-(1e-6));
    else
        psi_T = min(max(1e-6,model.B(2)/(-model.B(1))),1-(1e-6));
    end
    en = model.NoiseVariance;

    y = table2array(px_nav(:,{'lprice','lnav'}))';
    InitialParams = [InvMaxZero(1,psi_u),InvMaxZero(1,psi_p),InvMaxZero(1,psi_n),InvMaxZero(1,psi_T),InvMaxZero(-1,theta_u),InvMaxZero(-1,theta_p),sqrt([eu,ep,en]),0,0,0];
    LossToMinmize = @(vparam) -kf_simple_mymodel(vparam,y)/(size(y,2)-1);

    optionsIVAN = optimset('Display', 'iter-detailed','LargeScale', 'off','MaxFunEvals',5000);
    [EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);
    
    % [fhat,EstimParams] = csminwel(LossToMinmize,InitialParams,eye(size(InitialParams,1))*.5,[] ,1e-14,100);
    % options = optimoptions('patternsearch', ...
    %         'MaxIterations', 1000, ... 
    %         'UseParallel', true, ...
    %         'Display', 'iter',...
    %         'ConstraintTolerance', 1e-6,...
    %         'MeshTolerance', 1e-10,...
    %         'StepTolerance', 1e-10);

    % lb = [(-5)*ones(1,4), (-5)*ones(1,2), -Inf(1,3), zeros(1,3),-Inf(1,3)];
    % ub = [(5)*ones(1,4), (5)*ones(1,2), Inf(1,3), Inf(1,3),Inf(1,3)];

    % [EstimParams, fval] = patternsearch(LossToMinmize, InitialParams, [], [], [], [], lb, ub, options);

    [LL,Res] = kf_simple_mymodel(EstimParams,y);
    % [a_sm, v_sm] = kf_smooth_original(Res.v_t, Res.invF_t, Res.L_t, Res.alfa_t(:,1:(end-1)), Res.P_t(:,:,1:(end-1)), Res.Z);
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
    names = {"isin","p_jb", "p_arch", "p_lb"}; 
    test_tab = table(names(:), {isin_list{j}, p_jb, p_arch, p_lb}', 'VariableNames', {'Names', 'Value'});
    save(strcat('test/',px_nav.isin{j}), 'test_tab');
end