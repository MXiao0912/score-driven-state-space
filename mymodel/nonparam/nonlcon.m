function [c, ceq] = nonlcon(x)
    cov = [x(1) x(4) x(5); x(4) x(2) x(6); x(5) x(6) x(3)];
    c = 1e-10-eig(cov);
    ceq = [];
end