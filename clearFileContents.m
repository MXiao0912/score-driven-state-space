function clearFileContents(filename)
    % Open the file in write mode, which clears its contents
    fid = fopen(filename, 'w');
    if fid == -1
        error('Failed to open file: %s', filename);
    end
    fclose(fid);
end