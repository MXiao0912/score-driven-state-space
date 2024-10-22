function [Contribution] = ContributionLogLik(Vpar,fixedParams, y, mean_ft)

[LL,Res] = kf_tvp_corr1(insertParams(Vpar, fixedParams),y, mean_ft,[]); 
Contribution = Res.ContributionLogLik;
