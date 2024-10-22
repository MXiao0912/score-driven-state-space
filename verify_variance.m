function EstimParams = verify_variance(Simu)
    Simu = strtrim(Simu);
    seed = str2double(Simu); % Assuming Simu is a string containing the array task ID
    rng(seed); % Set the random seed to ensure different random numbers

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

    InitialParams = [kap_hes;c;A;B;mom_];
    fixedParams = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1); NaN(1); NaN(4,1); NaN(4,1); NaN(4,1); NaN(1)];

    y = table2array(obs(:,{'price','nav'}))';

    [EstimParams, FullParams] = sim_est_corr1(y, mean_ft, "patternsearch", false, InitialParams, fixedParams);


    save(strcat('verify_var_res/EstimParams_res_', Simu, '.mat'), "EstimParams");
end