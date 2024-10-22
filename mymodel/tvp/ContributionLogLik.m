function [Contribution] = ContributionLogLik(Vpar, y, mean_ft)

[LL,Res] = kf_tvp_mymodel(Vpar,y, mean_ft); 
Contribution = Res.ContributionLogLik;
