function res_obj = moment_rec(vpar,k_est, params)
sigma_u2 = vpar(1);
sigma_p2 = vpar(2);
sigma_n2 = vpar(3);
cov_up = vpar(4);
cov_un = vpar(5);
cov_pn = vpar(6);


knn = k_est(1);
kpp = k_est(2);
knp = k_est(3);
knln = k_est(4);
kplp = k_est(5);
knlp = k_est(6);
kpln = k_est(7);

c1 = params(1);
c2 = params(2);
c3 = params(3);
c4 = params(4);
c5 = params(5);
c6 = params(6);
c7 = params(7);
c8 = params(8);
c9 = params(9);
c10 = params(10);
c11 = params(11);
c12 = params(12);
c13 = params(13);
c14 = params(14);
c15 = params(15);
c16 = params(16);
c17 = params(17);
c18 = params(18);
c19 = params(19);
c20 = params(20);
c21 = params(21);
c22 = params(22);
c23 = params(23);
c24 = params(24);


eq1 = kpp-(c1*sigma_p2+c2*sigma_u2+c3*cov_up);

eq2 = knp-(c4*cov_pn+c5*cov_un+c6*cov_up+c7*sigma_u2);

eq3 = knn-(c8*sigma_n2+c9*cov_un+c10*sigma_u2);

eq4 = kplp-(c11*sigma_p2+c12*cov_up+c13*sigma_u2);

eq5 = kpln-(c14*cov_pn+c15*cov_un+c16*cov_up+c17*sigma_u2);

eq6 = knlp-(c18*cov_pn+c19*cov_un+c20*cov_up+c21*sigma_u2);

eq7 = knln-(c22*sigma_n2+c23*cov_un+c24*sigma_u2);

res_obj = eq1^2+eq2^2+eq3^2+eq4^2+eq5^2+eq6^2+eq7^2;

end