j=12;

cd("/rds/user/mx235/hpc-work/px_nav_ssm/mytvp/mymodel/tvp/")
addpath("../")
run Useful_Transformations.m

px_nav_tot = readtable("../../../mydata/px_nav.csv");
isin_list = unique(px_nav_tot.isin);

subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

for i = 1:length(isin_list)
    subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
end

px_nav = subTables{j};

y = load(strcat('yall/',px_nav.isin{j})).y;
mean_ft = load(strcat('mean_ft/',px_nav.isin{j})).mean_ft;
EstimParams = load(strcat('rawestimates/',px_nav.isin{j})).EstimParams;
Res = load(strcat('Resall/',px_nav.isin{j})).Res;

sim_se_fct(y, mean_ft, EstimParams, Res, px_nav.isin{j})；

Save_Collect = load(strcat('ft_se/',px_nav.isin{j})).Save_Collect;
color = "g";
alpha=0.3;
fig = figure();
subplot(6,1,1);
plot(px_nav.date(5:end), squeeze(Res.Q(1,1,6:end)));
hold on;
confplot(px_nav.date(5:end), squeeze(Res.Q(1,1,6:end)),squeeze(mad(Save_Collect.varu(5:end,:),0,2) * 1.96)', color, alpha);
subplot(6,1,2);
plot(px_nav.date(5:end), squeeze(Res.Q(2,2,6:end)));
hold on;
confplot(px_nav.date(5:end), squeeze(Res.Q(2,2,6:end)),squeeze(mad(Save_Collect.varp(5:end,:),0,2) * 1.96)', color, alpha);
subplot(6,1,3);
plot(px_nav.date(5:end), squeeze(Res.Q(3,3,6:end)));
hold on;
confplot(px_nav.date(5:end), squeeze(Res.Q(3,3,6:end)),squeeze(mad(Save_Collect.varn(5:end,:),0,2) * 1.96)', color, alpha);
subplot(6,1,4);
plot(px_nav.date(5:end), squeeze(Res.Q(2,1,6:end)./(Res.Q(1,1,6:end).*Res.Q(2,2,6:end)).^(0.5)));
hold on;
confplot(px_nav.date(5:end), squeeze(Res.Q(2,1,6:end)./(Res.Q(1,1,6:end).*Res.Q(2,2,6:end)).^(0.5)), squeeze(mad(Save_Collect.corrpu(5:end,:),0,2) * 1.96)', color, alpha);
subplot(6,1,5);
plot(px_nav.date(5:end), squeeze(Res.Q(3,1,6:end)./(Res.Q(1,1,6:end).*Res.Q(3,3,6:end)).^(0.5)));
hold on;
confplot(px_nav.date(5:end), squeeze(Res.Q(3,1,6:end)./(Res.Q(1,1,6:end).*Res.Q(3,3,6:end)).^(0.5)), squeeze(mad(Save_Collect.corrnu(5:end,:),0,2) * 1.96)', color, alpha);
subplot(6,1,6);
plot(px_nav.date(5:end), squeeze(Res.Q(3,2,6:end)./(Res.Q(2,2,6:end).*Res.Q(3,3,6:end)).^(0.5)));
hold on;
confplot(px_nav.date(5:end), squeeze(Res.Q(2,3,6:end)./(Res.Q(3,3,6:end).*Res.Q(2,2,6:end)).^(0.5)), squeeze(mad(Save_Collect.corrpn(5:end,:),0,2) * 1.96)', color, alpha);

saveas(fig, strcat('q_se_plot/',px_nav.isin{j},'.png'));