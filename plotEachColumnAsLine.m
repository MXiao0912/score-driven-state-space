function plotEachColumnAsLine(data)
    % Function to plot each column as a separate line
    % The number of lines plotted will equal the number of rows in the data
    %
    % INPUT:
    % data - matrix where each row represents a time point and each column
    %        represents a different sample or measurement at that time point
    
    % Create a vector for the time points (assuming time points are row indices)
    time_points = 1:size(data, 1);
    
    % Plot each column as a line
    figure; % Create a new figure
    plot(time_points, data, 'LineWidth', 1.5);
    
    % Label the axes
    xlabel('Time Points');
    ylabel('Measurement Value');
    
    % Add a title to the plot
    title('Plot of Each Column as a Separate Line');
    
    % Optionally, add grid lines for better visualization
    grid on;
end