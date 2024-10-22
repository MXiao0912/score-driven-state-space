function [c, ceq] = nonlcon_tvp(x)
    cov = [x(1) 0 x(4); 0 x(2) x(5); x(4) x(5) x(3)];
    c = 1e-10-eig(cov);
    ceq = [];
end