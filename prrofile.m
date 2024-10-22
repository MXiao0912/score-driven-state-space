[psi, phi, initF1] = nonparam_tvp(obs);

 [PositiveTrans, InvPositiveTrans, DerPositiveTrans, ...
     UnoMenoUno, InvUnoMenoUno, DerUnoMenoUno, ...
     LogisticFun, MaxZero, InvMaxZero, DerMaxZero, vec] = Useful_Transformations();

% user-specify A and B
A = InvUnoMenoUno(0.9*ones(5,1)); 
B = InvMaxZero(0.1, [0.015*ones(5,1)]);
kap_hes = InvMaxZero(.5,.015);
c = (eye(5)-diag(UnoMenoUno(A)))*initF1;

InitialParams = [InvUnoMenoUno(psi);InvUnoMenoUno(phi);zeros(3,1); kap_hes; c; A;B];
% fixedParams = [NaN(2,1);NaN(3,1); NaN(1); NaN(5,1); NaN(5,1); NaN(5,1)];
fixedParams = InitialParams;
% fixedIndices = [3:5,7:11,12:21];

y = table2array(obs(2:end,{'price','nav'}))';

wrappedObjectiveFunction = @(variableParams) ...
    kf_tvp(insertParams(variableParams, fixedParams),y, initF1)/(size(y,2)-1);

% % TVP
LossToMinimize = @(vparam) -wrappedObjectiveFunction(vparam);
% optionsIVAN = optimoptions('fminunc', 'Display', 'iter-detailed','MaxFunEvals',1000, 'StepTolerance', 1e-10);

options = optimoptions('patternsearch', ...
    'UseParallel', true, ...
    'Display', 'iter');

% [EstimParams, fval] = patternsearch(LossToMinimize, InitialParams, [], [], [], [], [], [], [], options);


% [fhat,EstimParams] = csminwel(LossToMinmize,InitialParams,eye(size(InitialParams,1))*.5,[] ,1e-20,1000);
% [EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);

% profile likelihood
tolerance = 1e-6;
maxIterations = 100;
converged = false;
iteration = 0;

currentParams = InitialParams;
previousParams = currentParams;

group1Indices = 1:5; % Indices for mean variables
group2Indices = 6:21; % Indices for zeros(3,1)

while ~converged && iteration < maxIterations
    iteration = iteration + 1;
    fprintf('Iteration %d\n', iteration);
    
    % Optimize group 1
    fixedParams(group1Indices) = NaN;
    [currentParams(group1Indices), fval] = patternsearch(LossToMinimize, currentParams(group1Indices), [], [], [], [], [], [], [], options);
    fixedParams(group1Indices) = currentParams(group1Indices);
    
    % Optimize group 2
    fixedParams(group2Indices) = NaN;
    [currentParams(group2Indices), fval] = patternsearch(LossToMinimize, currentParams(group2Indices), [], [], [], [], [], [], [], options);
    fixedParams(group2Indices) = currentParams(group2Indices);

    paramChange = norm(currentParams - previousParams);
    if paramChange < tolerance
        converged = true;
    else
        previousParams = currentParams;
    end
end

EstimParams = currentParams;


EstimParams = insertParams(EstimParams, fixedParams);


psi = UnoMenoUno(EstimParams(1))
phi = UnoMenoUno(EstimParams(2))
mean = EstimParams(3:5)
kap_hes = MaxZero(.5,EstimParams(6))
c = EstimParams(7:11)
A = UnoMenoUno(EstimParams(12:16))
B = MaxZero(0.1, EstimParams(17:21))


% solution 1: use score without scaling : try with the original scaling (more sensitive to later scores)
% solution 2: reduce the number of tvp to 3 if Fdot only have 3 unique rows
% solution 3: more smoothing
% question: why originally we think all 5 params could be identified?
