function [c, ceq] = myNonlinearConstraints(x)

     [PositiveTrans, InvPositiveTrans, DerPositiveTrans, ...
     UnoMenoUno, InvUnoMenoUno, DerUnoMenoUno, ...
     LogisticFun, MaxZero, InvMaxZero, DerMaxZero, vec] = Useful_Transformations();

    % Inequality constraints
    c_ = x(2:6);
    A_ = x(7:11);
    mean_ft_ = (eye(size(A_,1))-diag(MaxZero(0.99, A_)))\c_;
    c = [exp(mean_ft_(1))^2-0.2^2;
         exp(mean_ft_(2))^2-0.2^2;
         exp(mean_ft_(3))^2-0.2^2]; % Example: 1 - x1*x2 <= 0

    % Equality constraints
    ceq = []; % No equality constraints in this example
end