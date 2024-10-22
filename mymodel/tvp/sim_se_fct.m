function Save_Collect = sim_se_fct(y, mean_ft, EstimParams, Res,isin_)
    ContrLossToMinmize = @(vparam) ContributionLogLik(vparam, y, mean_ft);

    LossToMinmize = @(vparam) -kf_tvp_mymodel(vparam,y, mean_ft);

    TT = size(y,2);

    HHessian = (1/TT)*fdhess2(LossToMinmize, EstimParams);
    jac = fdjacob(ContrLossToMinmize, EstimParams);
    HHessianSand = (1/TT)*jac'*jac; 

    % WHY IS JAC 0 FOR SOME VARIABLES FOR ALL PERIODS?

    [TransParam,Jacob] = TransformParametersLik(EstimParams);
    TransParam = TransParam';

    VARIANZA_Coeffs =(1/TT)*inv(HHessianSand+1e-10*eye(size(HHessianSand)));
    VARIANZA_Coeffs = (VARIANZA_Coeffs+VARIANZA_Coeffs')/2;
    VARIANZARob_Coeffs =(1/TT)*inv(HHessian+1e-6*eye(size(HHessianSand)))*HHessianSand*inv(HHessian+1e-6*eye(size(HHessianSand)));
    VARIANZARob_Coeffs=(VARIANZARob_Coeffs+VARIANZARob_Coeffs')/2;

    VARIANZA_TransfCoeffs =diag(Jacob)*inv(HHessianSand+1e-6*eye(size(HHessianSand)))*diag(Jacob);

    % VARIANZARob_TransfCoeffs =(1/TT)*diag(Jacob)*inv(HHessian+1e-6*eye(size(HHessianSand)))*HHessianSand*inv(HHessian+1e-6*eye(size(HHessianSand)))*diag(Jacob);

    PickNoSimulations = 1000; 

    % Pre-allocate Save_Collect structure with NaNs
    Save_Collect.varu = NaN(TT, PickNoSimulations);
    Save_Collect.varp = NaN(TT, PickNoSimulations);
    Save_Collect.varn = NaN(TT, PickNoSimulations);
    Save_Collect.corrpu = NaN(TT, PickNoSimulations);
    Save_Collect.corrnu = NaN(TT, PickNoSimulations);
    Save_Collect.corrpn = NaN(TT, PickNoSimulations);
    Save_Collect.TransfParam = NaN(size(TransParam, 1), PickNoSimulations);
    Save_Collect.Param = NaN(size(EstimParams, 1), PickNoSimulations);

    % Fill the first column with the initial results
    Save_Collect.varu(:, 1) = squeeze(Res.Q(1, 1, 2:end));
    Save_Collect.varp(:, 1) = squeeze(Res.Q(2, 2, 2:end));
    Save_Collect.varn(:, 1) = squeeze(Res.Q(3, 3, 2:end));
    Save_Collect.corrnu(:, 1) = squeeze(Res.Q(3, 1, 2:end) ./ sqrt(Res.Q(1, 1, 2:end) .* Res.Q(3, 3, 2:end)));
    Save_Collect.corrpn(:, 1) = squeeze(Res.Q(3, 2, 2:end) ./ sqrt(Res.Q(3, 3, 2:end) .* Res.Q(2, 2, 2:end)));
    Save_Collect.corrpu(:, 1) = squeeze(Res.Q(2, 1, 2:end) ./ sqrt(Res.Q(1, 1, 2:end) .* Res.Q(2, 2, 2:end)));
    Save_Collect.TransfParam(:, 1) = TransParam;
    Save_Collect.Param(:, 1) = EstimParams;

    % Fix seed 
    rng(1000);

    % Pre-allocate arrays to store results temporarily
    varu_all = NaN(TT, PickNoSimulations); 
    varp_all = NaN(TT, PickNoSimulations);
    varn_all = NaN(TT, PickNoSimulations);
    corrpu_all = NaN(TT, PickNoSimulations);
    corrnu_all = NaN(TT, PickNoSimulations);
    corrpn_all = NaN(TT, PickNoSimulations);
    param_all = NaN(size(EstimParams, 1), PickNoSimulations);
    transf_param_all = NaN(size(TransParam, 1), PickNoSimulations);

    parfor Simu = 2:PickNoSimulations
        disp(Simu);
        CheckFinite = 0;
        
        while CheckFinite == 0
            % Generate parameters
            Paramter_i = mvnrnd(EstimParams, VARIANZARob_Coeffs);
            
            % Run the filter
            [LL_i, Res_i] = kf_tvp_mymodel(Paramter_i,y, mean_ft);
            
            if LL_i ~= -1.0000e+10
                % Extract results
                varu_i = squeeze(Res_i.Q(1, 1, 2:end));
                varp_i = squeeze(Res_i.Q(2, 2, 2:end));
                varn_i = squeeze(Res_i.Q(3, 3, 2:end));
                corrnu_i = squeeze(Res_i.Q(3, 1, 2:end) ./ sqrt(Res_i.Q(1, 1, 2:end) .* Res_i.Q(3, 3, 2:end)));
                corrpn_i = squeeze(Res_i.Q(2, 3, 2:end) ./ sqrt(Res_i.Q(3, 3, 2:end) .* Res_i.Q(2, 2, 2:end)));
                corrpu_i = squeeze(Res_i.Q(2, 1, 2:end) ./ sqrt(Res_i.Q(1, 1, 2:end) .* Res_i.Q(2, 2, 2:end)));
                
                [TransParam, Jacob] = TransformParametersLik(Paramter_i);
                paramTransform_i = TransParam';
                
                % Check the conditions
                if all(varu_i < 10) && all(varp_i < 10) && all(varn_i < 10)
                    % Store results in temporary arrays
                    varu_all(:, Simu) = varu_i;
                    varp_all(:, Simu) = varp_i;
                    varn_all(:, Simu) = varn_i;
                    corrpu_all(:, Simu) = corrpu_i;
                    corrnu_all(:, Simu) = corrnu_i;
                    corrpn_all(:, Simu) = corrpn_i;
                    param_all(:, Simu) = Paramter_i;
                    transf_param_all(:, Simu) = paramTransform_i;
                    
                    CheckFinite = 1;
                end
            end
        end
    end

    % Collect results into Save_Collect after the parfor loop
    Save_Collect.varu(:, 2:end) = varu_all(:, 2:end);
    Save_Collect.varp(:, 2:end) = varp_all(:, 2:end);
    Save_Collect.varn(:, 2:end) = varn_all(:, 2:end);
    Save_Collect.corrpu(:, 2:end) = corrpu_all(:, 2:end);
    Save_Collect.corrnu(:, 2:end) = corrnu_all(:, 2:end);
    Save_Collect.corrpn(:, 2:end) = corrpn_all(:, 2:end);

    Save_Collect.Param(:, 2:end) = param_all(:, 2:end);
    Save_Collect.TransfParam(:, 2:end) = transf_param_all(:, 2:end);

    save(strcat('ft_se/',isin_), "Save_Collect");

end