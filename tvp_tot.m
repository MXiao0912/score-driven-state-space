clear;
cd '/Users/mingmeixiao/code/fca_price_nav/mytvp'
run Useful_Transformations;
addpath("../mydata/")
addpath("../Sims Solver/")

px_nav_tot = readtable("../mydata/px_nav.csv");
isin_list = unique(px_nav_tot.isin);

subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

for i = 1:length(isin_list)
    subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
end

j=1;
% Non-parametric estimation of F1:  and initial parameter values
px_nav = subTables{j};
px_nav.lprice = log(px_nav.price);
px_nav.lnav = log(px_nav.nav);
px_nav.dprice = [NaN; diff(px_nav.lprice)];
px_nav.dnav = [NaN; diff(px_nav.lnav)];
px_nav.dprice = fillmissing(px_nav.dprice,'previous');
px_nav.dnav = fillmissing(px_nav.dnav,'previous');

% initialize with kernel method
[psi, phi, initF1] = nonparam_tvp(px_nav);


% user-specify A and B
A = InvMaxZero(1,.9*ones(5,1)); 
B = InvMaxZero(0.1, 0.015*ones(5,1));
kap_hes = InvMaxZero(.5,.015);
c = (eye(5)-diag(MaxZero(1,A)))*initF1;

% Init TVP
y = table2array(px_nav(:,{'lprice','lnav'}))';
InitialParams = [InvUnoMenoUno(psi);InvUnoMenoUno(phi);0;0;0;kap_hes;c;A;B];

% TVP
LossToMinmize = @(vparam) -kf_tvp(vparam,y, initF1)/(size(y,2)-1);
optionsIVAN = optimset('Display', 'iter-detailed','LargeScale', 'off','MaxFunEvals',5000);
[EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);


[fhat,EstimParams] = csminwel(LossToMinmize,InitialParams,eye(size(InitialParams,1))*.5,[] ,1e-14,100);




