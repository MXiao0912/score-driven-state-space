function [TransParam,Jacob] = TransformParametersLik(Vpar);

run Useful_Transformations

%% 
if size(Vpar,1) ==1 
    Vpar = Vpar'; 
end

%% SET UP THE DYNAMICS OF THE SCORE DRIVEN PARAMTERS 

% psi = UnoMenoUno(EstimParams(1))
% phi = UnoMenoUno(EstimParams(2))

% err_mean = EstimParams(3:5)

% kap_hes = MaxZero(.5,EstimParams(6))

% c = EstimParams(7:10)

% A = MaxZero(0.99, EstimParams(11:14))
% B = MaxZero(0.1, EstimParams(15:18))
% mom = MaxZero(1,EstimParams(19))

% a) Parameters of the score driven process 

TransParam(1) = MaxZero(1,Vpar(1)); %% psi, phi
Jacob(1) =DerMaxZero(1,Vpar(1));
TransParam(2) = MaxZero(1,Vpar(2)); %% psi, phi
Jacob(2) =DerMaxZero(1,Vpar(2));
TransParam(3) = MaxZero(1,Vpar(3)); %% psi, phi
Jacob(3) =DerMaxZero(1,Vpar(3));
TransParam(4) = MaxZero(1,Vpar(4)); %% psi, phi
Jacob(4) =DerMaxZero(1,Vpar(4));
TransParam(5) = MaxZero(-1,Vpar(5)); %% psi, phi
Jacob(5) =DerMaxZero(-1,Vpar(5));
TransParam(6) = MaxZero(-1,Vpar(6)); %% psi, phi
Jacob(6) =DerMaxZero(-1,Vpar(6));
TransParam(7) = MaxZero(1,Vpar(7)); %% psi, phi
Jacob(7) =DerMaxZero(1,Vpar(7));

TransParam(8:13) = Vpar(8:13); %% err_mean
Jacob(8:13) = ones(6,1); 

TransParam(14:19) = MaxZero(0.99, Vpar(14:19)); %% err_mean
Jacob(14:19) = DerMaxZero(0.99,Vpar(14:19));

TransParam(20:25) = MaxZero(0.1, Vpar(20:25)); %% err_mean
Jacob(20:25) = DerMaxZero(0.1,Vpar(20:25));

TransParam(26) = MaxZero(1, Vpar(26)); %% err_mean
Jacob(26) = DerMaxZero(1,Vpar(26));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TransParam(1) = MaxZero(0.5,Vpar(1)); %% kap_hes 
% Jacob(1) = DerMaxZero(0.5,Vpar(1)); 

% TransParam(2:5) = Vpar(2:5); %% c
% Jacob(2:5) = ones(4,1);

% TransParam(6:9) =  MaxZero(0.99, Vpar(6:9)); %% A
% Jacob(6:9) = DerMaxZero(0.99, Vpar(6:9));

% TransParam(10:13) = MaxZero(0.1, Vpar(10:13)); %% B
% Jacob(10:13) = DerMaxZero(0.1, Vpar(10:13));

% TransParam(14) = MaxZero(1,Vpar(14)); %% mom
% Jacob(14) = DerMaxZero(1,Vpar(14));


