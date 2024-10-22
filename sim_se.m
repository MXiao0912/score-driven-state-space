%% CALCULATE SE ESTIMATES WITH SIMULATIONS (coded up only 'AllCoeffs' specification)c
% LossToMinmize = @(vparam) TVKF_pd_Pd_Model([vparam],y',initialF1,[],SettingsTransfor);

ContrLossToMinmize = @(vparam) ContributionLogLik(vparam, fixedParams, y, mean_ft);

LossToMinmize = @(vparam) -kf_tvp_corr1(insertParams(vparam, fixedParams),y, mean_ft,[]);

TT = size(y,2)

HHessian = (1/TT)*fdhess2(LossToMinmize, EstimParams);

jac = fdjacob(ContrLossToMinmize, EstimParams);
HHessianSand = (1/TT)*jac'*jac; 

% WHY IS JAC 0 FOR SOME VARIABLES FOR ALL PERIODS?

[TransParam,Jacob] = TransformParametersLik(EstimParams);
TransParam = TransParam';

% VARIANZA_Coeffs =(1/TT)*inv(HHessianSand+1e-6*eye(size(HHessianSand)));
VARIANZARob_Coeffs =(1/TT)*inv(HHessian+1e-6*eye(size(HHessianSand)))*HHessianSand*inv(HHessian+1e-6*eye(size(HHessianSand)));
 
VARIANZA_TransfCoeffs =diag(Jacob)*inv(HHessianSand+1e-6*eye(size(HHessianSand)))*diag(Jacob);

VARIANZARob_TransfCoeffs =(1/TT)*diag(Jacob)*inv(HHessian+1e-6*eye(size(HHessianSand)))*HHessianSand*inv(HHessian+1e-6*eye(size(HHessianSand)))*diag(Jacob);

PickNoSimulations = 1000; 

% Pre-allocate Save_Collect structure with NaNs
Save_Collect.vare = NaN(TT-1, PickNoSimulations);
Save_Collect.varw = NaN(TT-1, PickNoSimulations);
Save_Collect.varr = NaN(TT-1, PickNoSimulations);
Save_Collect.correw = NaN(TT-1, PickNoSimulations);
Save_Collect.TransfParam = NaN(size(TransParam, 1), PickNoSimulations);
Save_Collect.Param = NaN(size(EstimParams, 1), PickNoSimulations);

% Fill the first column with the initial results
Save_Collect.vare(:, 1) = squeeze(Res.Q(1, 1, 2:end));
Save_Collect.varw(:, 1) = squeeze(Res.Q(2, 2, 2:end));
Save_Collect.varr(:, 1) = squeeze(Res.Q(3, 3, 2:end));
Save_Collect.correw(:, 1) = squeeze(Res.Q(2, 1, 2:end) ./ sqrt(Res.Q(1, 1, 2:end) .* Res.Q(2, 2, 2:end)));
Save_Collect.TransfParam(:, 1) = TransParam;
Save_Collect.Param(:, 1) = EstimParams;

% Fix seed 
rng(1000);

% Pre-allocate arrays to store results temporarily
vare_all = NaN(T-1, PickNoSimulations); 
varw_all = NaN(T-1, PickNoSimulations);
varr_all = NaN(T-1, PickNoSimulations);
correw_all = NaN(T-1, PickNoSimulations);
param_all = NaN(size(EstimParams, 1), PickNoSimulations);
transf_param_all = NaN(size(TransParam, 1), PickNoSimulations);

parfor Simu = 2:PickNoSimulations
    disp(Simu);
    CheckFinite = 0;
    
    while CheckFinite == 0
        % Generate parameters
        Paramter_i = mvnrnd(EstimParams, VARIANZARob_Coeffs);
        
        % Run the filter
        [LL_i, Res_i] = kf_tvp_corr1(insertParams(Paramter_i, fixedParams), y, mean_ft, []);
        
        if LL_i ~= -1.0000e+10
            % Extract results
            vare_i = squeeze(Res_i.Q(1, 1, 2:end));
            varw_i = squeeze(Res_i.Q(2, 2, 2:end));
            varr_i = squeeze(Res_i.Q(3, 3, 2:end));
            correw_i = squeeze(Res_i.Q(2, 1, 2:end) ./ sqrt(Res_i.Q(1, 1, 2:end) .* Res_i.Q(2, 2, 2:end)));
            
            [TransParam, Jacob] = TransformParametersLik(Paramter_i);
            paramTransform_i = TransParam';
            
            % Check the conditions
            if all(vare_i < 10) && all(varw_i < 10) && all(varr_i < 10)
                % Store results in temporary arrays
                vare_all(:, Simu) = vare_i;
                varw_all(:, Simu) = varw_i;
                varr_all(:, Simu) = varr_i;
                correw_all(:, Simu) = correw_i;
                param_all(:, Simu) = Paramter_i;
                transf_param_all(:, Simu) = paramTransform_i;
                
                CheckFinite = 1;
            end
        end
    end
end

% Collect results into Save_Collect after the parfor loop
Save_Collect.vare(:, 2:end) = vare_all(:, 2:end);
Save_Collect.varw(:, 2:end) = varw_all(:, 2:end);
Save_Collect.varr(:, 2:end) = varr_all(:, 2:end);
Save_Collect.correw(:, 2:end) = correw_all(:, 2:end);
Save_Collect.Param(:, 2:end) = param_all(:, 2:end);
Save_Collect.TransfParam(:, 2:end) = transf_param_all(:, 2:end);

save("savecollect.mat", "Save_Collect");
save("EstimParams.mat", "EstimParams");