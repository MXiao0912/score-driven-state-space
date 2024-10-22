cd("/rds/user/mx235/hpc-work/px_nav_ssm/mytvp/mymodel/tvp")
px_nav_tot = readtable("../../../mydata/px_nav.csv");
isin_list = unique(px_nav_tot.isin);

subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

for i = 1:length(isin_list)
    subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
end


tot_fund = table();
for i=1:length(isin_list)
    px_nav = subTables{i};
    dates = px_nav.date(2:end);
    Res = load(strcat("Resall/",isin_list(i), ".mat")).Res;
    state = squeeze(Res.alfa_t(1,3:end))';
    add_fund = table(dates, state);
    add_fund.isin = repmat(isin_list{i},size(state,1),1);
    tot_fund = vertcat(tot_fund, add_fund);
end

writetable(tot_fund, "fundvalue.csv");