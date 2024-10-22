function A_pinv = svd_pseudoinverse(A, tol)
    % SVD_PSEUDOINVERSE Calculates the pseudoinverse of a matrix using SVD.
    % 
    % Usage: 
    %   A_pinv = svd_pseudoinverse(A)
    %   A_pinv = svd_pseudoinverse(A, tol)
    %
    % Inputs:
    %   A - The input matrix
    %   tol - Tolerance for singular values (optional, default: 1e-9)
    %
    % Output:
    %   A_pinv - The pseudoinverse of matrix A
    
    if nargin < 2
        tol = 1e-9; % Default tolerance if not provided
    end

    % Perform SVD
    [U, S, V] = svd(A);

    % Initialize the inverse of S
    S_inv = zeros(size(S'));

    % Invert the diagonal elements of S above the tolerance
    for i = 1:min(size(S))
        if S(i, i) > tol
            S_inv(i, i) = 1 / S(i, i);
        end
    end

    % Form the pseudoinverse
    A_pinv = V * S_inv * U';
end