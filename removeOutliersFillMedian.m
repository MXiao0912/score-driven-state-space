function cleaned_data = removeOutliersFillMedian(data)
    % removeOutliersFillMedian removes outliers from the data using the IQR method
    % and replaces them with the median of the non-outlier data.
    %
    % Input:
    % - data: A numeric vector containing the time series data.
    %
    % Output:
    % - cleaned_data: The data with outliers replaced by the median.

    % Calculate the first and third quartiles and IQR
    Q1 = prctile(data, 10);
    Q3 = prctile(data, 90);
    IQR = Q3 - Q1;

    % Determine the bounds for identifying outliers
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;

    % Identify the outliers
    outliers = (data < lower_bound) | (data > upper_bound);

    % Compute the median of non-outlier data
    non_outlier_data = data(~outliers);
    median_value = median(non_outlier_data);

    % Replace outliers with the median
    cleaned_data = data;
    cleaned_data(outliers) = median_value;

    % Optionally, display the number of outliers detected
    fprintf('Detected and replaced %d outliers.\n', sum(outliers));
end