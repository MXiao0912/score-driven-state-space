% [psi, phi, initF1] = nonparam_tvp(obs);

 [PositiveTrans, InvPositiveTrans, DerPositiveTrans, ...
     UnoMenoUno, InvUnoMenoUno, DerUnoMenoUno, ...
     LogisticFun, MaxZero, InvMaxZero, DerMaxZero, vec] = Useful_Transformations();

% user-specify A and B
A = InvMaxZero(0.99, 0.9*ones(3,1)); 
B = InvMaxZero(0.1, [0.015*ones(3,1)]);
c = (eye(3)-diag(MaxZero(0.99, A)))*mean_ft(1:3);
kap_hes = InvMaxZero(.5,.015);
mom_ = InvMaxZero(1, 0.5);

InitialParams = [kap_hes;c;A;B;mom_];
fixedParams = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1); NaN(1); NaN(3,1);NaN(3,1); NaN(3,1);NaN(1)];
% fixedIndices = [3:5,7:11,12:21];

y = table2array(obs(2:end,{'price','nav'}))';

wrappedObjectiveFunction = @(variableParams) ...
    kf_tvp_variance_only(insertParams(variableParams, fixedParams),y, mean_ft(1:3))/(size(y,2)-1);

% TVP
LossToMinmize = @(vparam) -wrappedObjectiveFunction(vparam);
optionsIVAN = optimoptions('fminunc', 'Display', 'iter-detailed','MaxFunEvals',1000, 'StepTolerance', 1e-10);

options = optimoptions('patternsearch', ...
    'UseParallel', true, ...
    'Display', 'iter',...
    'ConstraintTolerance', 1e-6);

% [EstimParams, fval] = patternsearch(LossToMinmize, InitialParams, [], [], [], [], [], [], options);
% [fhat,EstimParams] = csminwel(LossToMinmize,InitialParams,eye(size(InitialParams,1))*.5,[] ,1e-20,100);
[EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);

EstimParams = insertParams(EstimParams, fixedParams);

psi = UnoMenoUno(EstimParams(1))
phi = UnoMenoUno(EstimParams(2))
err_mean = EstimParams(3:5)
kap_hes = MaxZero(.5,EstimParams(6))
c = EstimParams(7:9)
A = MaxZero(0.99, EstimParams(10:12))
B = MaxZero(0.1, EstimParams(13:15))
mom = MaxZero(1,EstimParams(16))
% c = (eye(5)-diag(A))*mean_ft

% solution 1: use score without scaling : try with the original scaling (more sensitive to later scores)
% solution 2: reduce the number of tvp to 3 if Fdot only have 3 unique rows
% solution 3: more smoothing
% question: why originally we think all 5 params could be identified?




% nonlinear condition
% inverse of information
