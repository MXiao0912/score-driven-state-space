function [EstimParams, FullParams] = sim_est_corr1(y, mean_ft, method, verbose, InitialParams, fixedParams, pred_ft)

addpath('../Sims_Solver/')
LossToMinmize = @(vparam) -kf_tvp_corr1(insertParams(vparam, fixedParams),y, mean_ft, pred_ft);

function [stop,options,optchanged] = displayXValues(optimValues,options,flag)
    % This function is called at each iteration
    x = optimValues.x; % Current point (x values)

    % Convert the vector to a string using arrayfun
    xStr = strjoin(arrayfun(@(xi) sprintf('%.6f', xi), x, 'UniformOutput', false), ' ');

    % Display the iteration number and x values
    fprintf('Iteration %d: x = [%s]\n', optimValues.iteration, xStr);

    % disp(['Iteration ', num2str(optimvalues.iteration), ': x = ', num2str(x)]);
    stop = false; % Return false so the algorithm continues
    optchanged = false;
end

function [stop,options,optchanged] = displayXValues_fminunc(x,options,flag)
    % This function is called at each iteration

    % Convert the vector to a string using arrayfun
    xStr = strjoin(arrayfun(@(xi) sprintf('%.6f', xi), x, 'UniformOutput', false), ' ');

    % Display the iteration number and x values
    fprintf(xStr);

    % disp(['Iteration ', num2str(optimvalues.iteration), ': x = ', num2str(x)]);
    stop = false; % Return false so the algorithm continues
    optchanged = false;
end

if method=="csminwel"

    % optionsIVAN = optimoptions('fminunc', 'Display', 'iter-detailed','MaxFunEvals',1000, 'StepTolerance', 1e-10,'OutputFcn', @displayXValues_fminunc);
    optionsIVAN = optimoptions('fminunc', 'Display', 'iter-detailed','MaxFunEvals',1000, 'StepTolerance', 1e-10);
    [fhat,EstimParams] = csminwel(LossToMinmize,InitialParams,eye(size(InitialParams,1))*.5,[] ,1e-20,1000);
    
elseif method=="fminunc"

    if verbose
        optionsIVAN = optimoptions('fminunc', 'Display', 'iter-detailed','MaxFunEvals',5000, 'StepTolerance', 1e-10,'OutputFcn', @displayXValues_fminunc);
    else
        optionsIVAN = optimoptions('fminunc', 'Display', 'iter-detailed','MaxFunEvals',5000, 'StepTolerance', 1e-10);
    end
    [EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);

elseif method =="patternsearch"
    
    if verbose
        options = optimoptions('patternsearch', ...
        'MaxIterations', 1000, ... 
        'UseParallel', true, ...
        'Display', 'iter',...
        'OutputFcn', @displayXValues, ...
        'ConstraintTolerance', 1e-6,...
        'MeshTolerance', 1e-10,...
        'StepTolerance', 1e-10);
    else
        options = optimoptions('patternsearch', ...
        'MaxIterations', 1000, ... 
        'UseParallel', true, ...
        'Display', 'iter',...
        'ConstraintTolerance', 1e-6,...
        'MeshTolerance', 1e-10,...
        'StepTolerance', 1e-10);
    end
    lb = [-Inf(1, 5), -Inf, (-2)*ones(1,4), (-5)*ones(1,4), (-3)*ones(1,4), -Inf];
    ub = [Inf(1, 5), Inf, (2)*ones(1,4), (5)*ones(1,4), (3)*ones(1,4), Inf];

    [EstimParams, fval] = patternsearch(LossToMinmize, InitialParams, [], [], [], [], lb, ub, options);
    % [EstimParams, fval] = patternsearch(LossToMinmize, InitialParams, [], [], [], [], [], [], options);

elseif method=="ga"
     [EstimParams, fval] = ga(LossToMinmize, 14);
end   

FullParams = insertParams(EstimParams, fixedParams);

% solution 1: use score without scaling : try with the original scaling (more sensitive to later scores)
% solution 2: reduce the number of tvp to 3 if Fdot only have 3 unique rows
% solution 3: more smoothing
% question: why originally we think all 5 params could be identified?


end