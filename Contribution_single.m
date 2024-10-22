function [Contribution_single] = Contribution_single(Vpar,fixedParams, y, mean_ft, m)

[LL,Res] = kf_tvp_corr1(insertParams(Vpar, fixedParams),y, mean_ft,[]); 
Contribution = Res.ContributionLogLik;
Contribution_single = Contribution(m);
