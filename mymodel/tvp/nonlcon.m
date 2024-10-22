function [c, ceq] = nonlcon(x)
    
    run Useful_Transformations.m
    c_ = x(8:13);
    A_ = MaxZero(0.99,x(14:19));
    % Inequality constraints (c(x) <= 0)
    c(1) = [c_(1)/(1-A_(1))-0.5];  % Example: x(1)^2 + x(2)^2 - 1; % A circle constraint (x1^2 + x2^2 ≤ 1)
    c(2) = [-c_(1)/(1-A_(1))-10]; 
    c(3) = [c_(2)/(1-A_(2))-0.5];  % Example: x(1)^2 + x(2)^2 - 1; % A circle constraint (x1^2 + x2^2 ≤ 1)
    c(4) = [-c_(2)/(1-A_(2))-10]; 
    c(5) = [c_(3)/(1-A_(3))-0.5];  % Example: x(1)^2 + x(2)^2 - 1; % A circle constraint (x1^2 + x2^2 ≤ 1)
    c(6) = [-c_(3)/(1-A_(3))-10]; 

    c(7) = [c_(4)/(1-A_(4))-3];
    c(8) = [-c_(4)/(1-A_(4))-3];
    c(9) = [c_(5)/(1-A_(5))-3];
    c(10) = [-c_(5)/(1-A_(5))-3];
    c(11) = [c_(6)/(1-A_(6))-3];
    c(12) = [-c_(6)/(1-A_(6))-3];
    
    % Equality constraints (ceq(x) == 0)
    ceq = [];  % Example: x(1)^2 + x(2) - 1; % Equality constraint (x1^2 + x2 = 1)
end