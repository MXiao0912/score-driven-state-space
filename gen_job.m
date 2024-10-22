% Define the range of numbers
numbers = (1:100)';

% Define the filename
filename = 'job.csv';

% Write the numbers to the CSV file
writematrix(numbers, filename);

% Display a message to confirm the file was created
disp(['File ', filename, ' created successfully.']);