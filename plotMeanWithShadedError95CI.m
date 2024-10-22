function plotMeanWithShadedError95CI(data)
    % Function to calculate the mean and 1.96 * standard error of each row
    % and plot the mean with a shaded area representing the 95% confidence interval
    %
    % INPUT:
    % data - matrix where each row represents a time point and each column
    %        represents a different sample or measurement at that time point
    
    % Calculate the mean of each row
    row_mean = mean(data, 2);
    
    % Calculate the standard error of each row
    row_std_error = std(data, 0, 2);
    
    % Calculate the 95% confidence interval (1.96 * standard error)
    % ci95 = 1.96 * row_std_error;
    ci95=0;
    
    % Create a vector for the time points (assuming time points are row indices)
    time_points = 1:size(data, 1);

    figure;
    
    % Plot the mean
    plot(time_points, row_mean, 'LineWidth', 2);
    hold on;
    
    color = "g";
    alpha=0.3;
    confplot(time_points, row_mean, ci95, color, alpha);
    
    % Label the axes
    xlabel('Time Points');
    ylabel('Mean Value');
    
    % Add a title to the plot
    title('Mean with Shaded 95% Confidence Interval over Time');
    
    % Optionally, add grid lines for better visualization
    grid on;
    
    % Turn off the hold to prevent affecting subsequent plots
    hold off;
end