% [psi, phi, initF1] = nonparam_tvp(obs);

 [PositiveTrans, InvPositiveTrans, DerPositiveTrans, ...
     UnoMenoUno, InvUnoMenoUno, DerUnoMenoUno, ...
     LogisticFun, MaxZero, InvMaxZero, DerMaxZero, vec] = Useful_Transformations();


y = table2array(obs(2:end,{'price','nav'}))';

% user-specify A and B
A = InvMaxZero(0.99, 0.9*ones(5,1)); 
B = InvMaxZero(0.1, [0.015*ones(5,1)]);
c = (eye(5)-diag(MaxZero(0.99, A)))*mean_ft;
kap_hes = InvMaxZero(.5,.015);
mom_ = InvMaxZero(1, 0.5);

% Define the initial parameter sets
InitialParams1 = [kap_hes; c(1:3); A(1:3); B(1:3); mom_];  % Example: First set of parameters to be optimized
InitialParams2 = [c(4:5); A(4:5); B(4:5)];        % Example: Second set of parameters to be optimized

% Define fixed parameter structures with NaN for positions that will be optimized
fixedParams1 = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1); NaN(1); NaN(3,1); InitialParams2(1:2); NaN(3,1); InitialParams2(3:4); NaN(3,1); InitialParams2(5:6); NaN(1)];
fixedParams2 = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1);InitialParams1(1); InitialParams1(2:4); NaN(2,1); InitialParams1(5:7); NaN(2,1); InitialParams1(8:10); NaN(2,1); InitialParams1(11)];

optionsPatternSearch = optimoptions('patternsearch', ...
    'UseParallel', true, ...
    'Display', 'iter', ...
    'MaxIterations', 300, ... 
    'ConstraintTolerance', 1e-6);


% Define convergence criteria and maximum number of iterations
maxIterations = 10;
tolerance = 1e-6;

% Initialize previous parameter values for convergence checking
prevParams1 = InitialParams1;
prevParams2 = InitialParams2;

for iter = 1:maxIterations
    % Step 1: Estimate the first set of parameters while keeping the second set fixed
    fixedParams1 = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1); NaN(1); NaN(3,1); prevParams2(1:2); NaN(3,1); prevParams2(3:4); NaN(3,1); prevParams2(5:6); NaN(1)];; % Update fixedParams1 with latest estimates of params2
    wrappedObjectiveFunction1 = @(variableParams1) ...
        kf_tvp_variance_only(insertParams(variableParams1, fixedParams1), y, mean_ft(1:3)) / (size(y,2)-1);

    LossToMinimize1 = @(vparam1) -wrappedObjectiveFunction1(vparam1);
    [EstimParams1, fval1] = patternsearch(LossToMinimize1, prevParams1, [], [], [], [], [], [], optionsPatternSearch);

    % Step 2: Estimate the second set of parameters while keeping the first set fixed
    fixedParams2 = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1);prevParams1(1); prevParams1(2:4); NaN(2,1); prevParams1(5:7); NaN(2,1); prevParams1(8:10); NaN(2,1); prevParams1(11)]; % Update fixedParams2 with latest estimates of params1
    wrappedObjectiveFunction2 = @(variableParams2) ...
        kf_tvp_variance_only(insertParams(variableParams2, fixedParams2), y, mean_ft(1:3)) / (size(y,2)-1);

    LossToMinimize2 = @(vparam2) -wrappedObjectiveFunction2(vparam2);
    [EstimParams2, fval2] = patternsearch(LossToMinimize2, prevParams2, [], [], [], [], [], [], optionsPatternSearch);

    % Check convergence by comparing parameter changes
    deltaParams1 = max(abs(EstimParams1 - prevParams1));
    deltaParams2 = max(abs(EstimParams2 - prevParams2));

    if deltaParams1 < tolerance && deltaParams2 < tolerance
        disp(['Converged at iteration ', num2str(iter)]);
        break;
    end

    % Update previous parameters for the next iteration
    prevParams1 = EstimParams1;
    prevParams2 = EstimParams2;

    disp(['Iteration ', num2str(iter), ' completed.']);
end

% InitialParams = [kap_hes;c;A;B;mom_];
% fixedParams = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1); NaN(1); mean_ft(1:3); NaN(2,1); InvMaxZero(0.99, zeros(3,1)); NaN(2,1); InvMaxZero(0.99, zeros(3,1)); NaN(2,1);NaN(1)];
% % fixedParams = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1); NaN(1); NaN(5,1); NaN(5,1); NaN(5,1);NaN(1)];
% % fixedIndices = [3:5,7:11,12:21];

% Combine the final estimated parameters
EstimParams = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1);EstimParams1(1); EstimParams1(2:4); EstimParams2(1:2); EstimParams1(5:7); EstimParams2(3:4); EstimParams1(8:10); EstimParams2(5:6); EstimParams1(11)];


% wrappedObjectiveFunction = @(variableParams) ...
%     kf_tvp(insertParams(variableParams, fixedParams),y, mean_ft)/(size(y,2)-1);

% % TVP
% LossToMinmize = @(vparam) -wrappedObjectiveFunction(vparam);
% optionsIVAN = optimoptions('fminunc', 'Display', 'iter-detailed','MaxFunEvals',1000, 'StepTolerance', 1e-10);

% options = optimoptions('patternsearch', ...
%     'UseParallel', true, ...
%     'Display', 'iter',...
%     'ConstraintTolerance', 1e-6);

% % [EstimParams, fval] = patternsearch(LossToMinmize, InitialParams, [], [], [], [], [], [], options);
% % [fhat,EstimParams] = csminwel(LossToMinmize,InitialParams,eye(size(InitialParams,1))*.5,[] ,1e-20,100);
% [EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);

% EstimParams = insertParams(EstimParams, fixedParams);

psi = UnoMenoUno(EstimParams(1))
phi = UnoMenoUno(EstimParams(2))
err_mean = EstimParams(3:5)
kap_hes = MaxZero(.5,EstimParams(6))
c = EstimParams(7:11)
A = MaxZero(0.99, EstimParams(12:16))
B = MaxZero(0.1, EstimParams(17:21))
mom = MaxZero(1,EstimParams(22))
% c = (eye(5)-diag(A))*mean_ft

% solution 1: use score without scaling : try with the original scaling (more sensitive to later scores)
% solution 2: reduce the number of tvp to 3 if Fdot only have 3 unique rows
% solution 3: more smoothing
% question: why originally we think all 5 params could be identified?




% nonlinear condition
% inverse of information
