cd("/rds/user/mx235/hpc-work/px_nav_ssm/mytvp/mymodel/tvp/")
run Useful_Transformations.m
addpath("../")

px_nav_tot = readtable("../../../mydata/px_nav.csv");
isin_list = unique(px_nav_tot.isin);

subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

for i = 1:length(isin_list)
    subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
end

test_res = table();
% for i = 1:length(isin_list)
for i = 1:length(isin_list)
    % Non-parametric estimation of F1:  and initial parameter values
    px_nav = subTables{i};
    data = load(strcat('test/',px_nav.isin{i})).test_tab;
    new_data = array2table(data.Value', 'VariableNames', string(data.Names));
    test_res = [test_res;new_data];
end

% Loop through each column of the table
for col = 1:width(test_res)
    % Check if the first element in the column is numeric
    if isnumeric(test_res{1, col}{1})
        test_res.(col) = cell2mat(test_res{:, col});  % Convert cells to numeric array
    elseif ischar(test_res{1, col}{1}) || isstring(test_res{1, col}{1})
        test_res.(col) = string(test_res{:, col});  % Convert cells to string array
    end
end

summary(test_res)         

% Initialize a string to store LaTeX code
latex_str = '\begin{table}[ht]\centering\n\begin{tabular}{|c|c|c|c|c|c|}\n';
latex_str = [latex_str 'Column & Min (0%) & Q1 (25%) & Median (50%) & Q3 (75%) & Max (100%) \\ \hline\n'];

% Now compute detailed quartiles for each numeric column
for i = 1:width(test_res)
    if isnumeric(test_res{:, i})
        quartiles = quantile(test_res{:, i}, [0 0.25 0.5 0.75 1]);
        % Add each row to LaTeX string
        latex_str = [latex_str sprintf('%s & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n', ...
            test_res.Properties.VariableNames{i}, quartiles(1), quartiles(2), quartiles(3), quartiles(4), quartiles(5))];
    end
end

% End the LaTeX table
latex_str = [latex_str '\end{tabular}\n\caption{Detailed Quartiles}\n\end{table}'];

% Display the LaTeX code (or write it to a file)
disp(latex_str);

