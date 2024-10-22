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

TransParam(1) = UnoMenoUno(Vpar(1)); %% psi, phi
Jacob(1) =DerUnoMenoUno(Vpar(1));
TransParam(2) = UnoMenoUno(Vpar(2)); %% psi, phi
Jacob(2) =DerUnoMenoUno(Vpar(2));

TransParam(3:5) = Vpar(3:5); %% err_mean
Jacob(3:5) = ones(3,1); 


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


