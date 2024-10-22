 [PositiveTrans, InvPositiveTrans, DerPositiveTrans, ...
    UnoMenoUno, InvUnoMenoUno, DerUnoMenoUno, ...
    LogisticFun, MaxZero, InvMaxZero, DerMaxZero, vec] = Useful_Transformations();


y = table2array(obs(2:end,{'price','nav'}))';
A = InvMaxZero(0.99, [0.9*ones(3,1)]); 
B = InvMaxZero(0.1, [0.015*ones(3,1)]);
kap_hes = 0.02;
c = f_intercept(1:3);

true_param = [InvUnoMenoUno(0.6);InvUnoMenoUno(0.8);zeros(3,1); kap_hes; c; A; B];

[LL, Res]=kf_tvp_variance_only(true_param,y, mean_ft(1:3));



% restrict phi further
% seperate the check into the mean and covariance parts