cd("/rds/user/mx235/hpc-work/px_nav_ssm/mytvp/mymodel/nonparam")

px_nav_tot = readtable("../../../mydata/px_nav.csv");
isin_list = unique(px_nav_tot.isin);

subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

for i = 1:length(isin_list)
    subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
end

% j=1
% px_nav = subTables{j};
% isin_data = load(strcat('time_varying_cov/',px_nav.isin{j}), 'M').M;
% sum_matrix = zeros(size(isin_data));
T=table()
for j=1:137
    disp(j)
    px_nav = subTables{j};
    isin_data = load(strcat('time_varying_cov/',px_nav.isin{j}), 'M').M;
    % 1. Repeat ISIN for all rows (ensure it's a character or string array)
    isin_col = repmat(px_nav.isin{1}, size(isin_data, 1), 1);  % Repeat ISIN for each row

    % 2. Extract and keep px_nav.date as it is (assuming it's a datetime array)
    date_col = px_nav.date(5:end);  % Subset of px_nav.date starting from the 5th element

    % 3. Concatenate the ISIN, date, and isin_data into a table
    new_row  = table(isin_col, date_col, isin_data);

    new_row = [new_row(:, {'isin_col', 'date_col'}), array2table(new_row.isin_data, 'VariableNames', {'\sigma_\mu^2', '\sigma_p^2', '\sigma_n^2', '\rho_{p\mu}', '\rho_{n\mu}', '\rho_{pn}'})];

    % Optionally, give appropriate column names
    new_row.Properties.VariableNames ={'ISIN', 'Date', '\sigma_\mu^2', '\sigma_p^2', '\sigma_n^2', '\rho_{p\mu}', '\rho_{n\mu}', '\rho_{pn}'};

    T = [T; new_row];
end


% Get the column names to include for the mean (excluding 'ISIN')
columns_to_mean = T.Properties.VariableNames(~ismember(T.Properties.VariableNames, {'ISIN', 'Date'}));

% Perform group-by on 'Date' and calculate the mean for the selected numeric columns
grouped_means = groupsummary(T, 'Date', 'mean', columns_to_mean);


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

saveas(gcf, strcat("mean",'.png'))