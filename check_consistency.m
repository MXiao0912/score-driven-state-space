all_var = NaN(14, 10);
for i=1:10

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate dgp
    T = i+100;
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
    A = InvMaxZero(0.9, 0.8*ones(4,1)); 
    B = InvMaxZero(0.1, [0.015*ones(4,1)]);
    c = (eye(4)-diag(MaxZero(0.99, A)))*mean_ft;
    kap_hes = InvMaxZero(.5,.015);
    mom_ = InvMaxZero(1, 0.5);
    
    InitialParams = [kap_hes;A;B;mom_];
    fixedParams = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1); NaN(1); zeros(4,1); NaN(4,1); NaN(4,1); NaN(1)];
    
    y = table2array(obs(:,{'price','nav'}))';
    [EstimParams, FullParams] = sim_est_corr1(y, mean_ft, "fminunc", false, InitialParams, fixedParams);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % estimate variance of parameters
    ContrLossToMinmize = @(vparam) ContributionLogLik(vparam, fixedParams, y, mean_ft);
    
    % HHessian = fdhess2(LossToMinmize, EstimParamsNew);
    jac = fdjacob(ContrLossToMinmize, EstimParams);
    HHessianSand = jac'*jac; 

    VARIANZA_Coeffs =inv(HHessianSand);

    all_var(:,i) = diag(VARIANZA_Coeffs);
end

