function res_obj = moment_rec(vpar,k_est, params)
sigma_ep2 = vpar(2);
sigma_r2 = vpar(1);
sigma_omega2 = vpar(3);
cov_eo = vpar(6);
cov_ro = vpar(5);
cov_re = vpar(4);

knn = k_est(1);
kpp = k_est(2);
knp = k_est(3);
knln = k_est(4);
kplp = k_est(5);
knlp = k_est(6);
kpln = k_est(7);

psi = params(1);
phi = params(2);

% eq1 = kpp-((1+psi^2)*sigma_r2+2*sigma_ep2);
% eq2 = knp-((1-phi)*sigma_r2+2*cov_eo);
% eq3 = knn-(((1-phi)^2)*sigma_r2+2*sigma_omega2);
% eq4 = kplp-(-sigma_ep2-psi*sigma_r2);
% eq5 = kpln-(-psi*(1-phi)*(sigma_r2)-cov_eo);
% eq6 = knlp-(-cov_eo);
% eq7 = knln-(-sigma_omega2);

eq1 = kpp-((1+psi^2)*sigma_r2+2*sigma_ep2+(1+psi)*cov_re);
eq2 = knp-((1-phi)*sigma_r2+2*cov_eo+(1+psi)*cov_ro+(1-phi)*cov_re);
eq3 = knn-(((1-phi)^2)*sigma_r2+2*sigma_omega2+(1-phi)*cov_ro);
eq4 = kplp-(-sigma_ep2-psi*sigma_r2-(psi+1)*cov_re);
eq5 = kpln-(-psi*(1-phi)*(sigma_r2)-cov_eo-(1-phi)*cov_re);
eq6 = knlp-(-cov_eo-cov_ro);
eq7 = knln-(-sigma_omega2-(1-phi)*cov_ro);

res_obj = eq1^2+eq2^2+eq3^2+eq4^2+eq5^2+eq6^2+eq7^2;

end