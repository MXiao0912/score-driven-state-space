cd("/rds/user/mx235/hpc-work/px_nav_ssm/mytvp/mymodel/tvp/")

px_nav_tot = readtable("../../../mydata/px_nav.csv");
isin_list = unique(px_nav_tot.isin);

subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

for i = 1:length(isin_list)
    subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
end

j=1
px_nav = subTables{j};
Res = load(strcat('Resall/',px_nav.isin{j}), 'Res').Res;
Qu = squeeze(Res.Q(1,1,6:end));
Qp = squeeze(Res.Q(2,2,6:end));
Qn = squeeze(Res.Q(3,3,6:end));
Qpu = squeeze(Res.Q(2,1,6:end));
Qnu = squeeze(Res.Q(3,1,6:end));
Qpn = squeeze(Res.Q(3,2,6:end));
isin_col = repmat(px_nav.isin{1}, size(Qu , 1), 1);
isin_tab = table(isin_col,px_nav.date(5:end), Qu,Qp,Qn,Qpu,Qnu,Qpn);



T=table()
for j=1:137
    disp(j)
    px_nav = subTables{j};
    Res = load(strcat('Resall/',px_nav.isin{j}), 'Res').Res;
    Qu = squeeze(Res.Q(1,1,6:end));
    Qp = squeeze(Res.Q(2,2,6:end));
    Qn = squeeze(Res.Q(3,3,6:end));
    Qpu = squeeze(Res.Q(2,1,6:end))./sqrt((squeeze(Res.Q(2,2,6:end)).*squeeze(Res.Q(1,1,6:end))));
    Qnu = squeeze(Res.Q(3,1,6:end))./sqrt((squeeze(Res.Q(3,3,6:end)).*squeeze(Res.Q(1,1,6:end))));
    Qpn = squeeze(Res.Q(3,2,6:end))./sqrt((squeeze(Res.Q(2,2,6:end)).*squeeze(Res.Q(3,3,6:end))));
    isin_col = repmat(px_nav.isin{1}, size(Qu , 1), 1);
    isin_tab = table(isin_col,px_nav.date(5:end), Qu,Qp,Qn,Qpu,Qnu,Qpn);

    % Optionally, give appropriate column names
    isin_tab.Properties.VariableNames ={'ISIN', 'Date', '\sigma_\mu^2', '\sigma_p^2', '\sigma_n^2', '\rho_{p\mu}', '\rho_{n\mu}', '\rho_{pn}'};

    T = [T; isin_tab];
end


% Get the column names to include for the mean (excluding 'ISIN')
columns_to_mean = T.Properties.VariableNames(~ismember(T.Properties.VariableNames, {'ISIN', 'Date'}));

% Perform group-by on 'Date' and calculate the mean for the selected numeric columns
grouped_means = groupsummary(T, 'Date', 'median', columns_to_mean);


t_ = (1/size(grouped_means,1)):(1/size(grouped_means,1)):1;
silverman = @(n) 1.06*min(std(t_),iqr(t_)/1.34)/(n^0.2);
h = silverman(size(grouped_means,1))*size(grouped_means,1);

grouped_means{:,3} = smoothdata(grouped_means{:,3}, 'gaussian',h/2);
grouped_means{:,4} = smoothdata(grouped_means{:,4}, 'gaussian',h/2);
grouped_means{:,5} = smoothdata(grouped_means{:,5}, 'gaussian',h/2);
grouped_means{:,6} = smoothdata(grouped_means{:,6}, 'gaussian',h/2);
grouped_means{:,7} = smoothdata(grouped_means{:,7}, 'gaussian',h/2);
grouped_means{:,8} = smoothdata(grouped_means{:,8}, 'gaussian',h/2);



figure;  % Open a new figure window
subplot(2,3,1);
%plot(px_nav.date(2:end), abs(M(:,4)), 'o-', MarkerSize=1);  % Change 'o-' to any other plot style as needed
plot(grouped_means.Date, grouped_means{:,3}, '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
xlabel('Date');  % Label x-axis
ylabel('\sigma_\mu^2');  % Label y-axis

subplot(2,3,2);
plot(grouped_means.Date, grouped_means{:,4}, '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
xlabel('Date');  % Label x-axis
ylabel('\sigma_p^2');  % Label y-axis

subplot(2,3,3);
plot(grouped_means.Date, grouped_means{:,5}, '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
xlabel('Date');  % Label x-axis
ylabel('\sigma_n^2');  % Label y-axis

subplot(2,3,4);
plot(grouped_means.Date, grouped_means{:,6}, '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
xlabel('Date');  % Label x-axis
ylabel('\rho_{p\mu}');  % Label y-axis

subplot(2,3,5);
plot(grouped_means.Date, grouped_means{:,7}, '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
xlabel('Date');  % Label x-axis
ylabel('\rho_{n\mu}');  % Label y-axis

subplot(2,3,6);
plot(grouped_means.Date, grouped_means{:,8}, '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
xlabel('Date');  % Label x-axis
ylabel('\rho_{pn}');  % Label y-axis

saveas(gcf, strcat("mean_score",'.png'))