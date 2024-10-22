function [c, ceq] = myNonlinearConstraints(x)
    % Inequality constraints
    c_ = x(2:6);
    A_ = x(7:11);
    mean_ft_ = (eye(size(A_,1))-diag(A_))\c_
    c = [mean_ft_]; % Example: 1 - x1*x2 <= 0

    % Equality constraints
    ceq = []; % No equality constraints in this example
end