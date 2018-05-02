Rin=0.22e-3;
Rout=0.9e-3;
epsilon_re=2.79;
mu0=4e-7*pi;
epsilon0=8.85e-12;
alpha_total0 = 0.172;
f0 = 1e8;

Zc0 = 138*sqrt(1./epsilon_re)*log10(Rout/Rin);

T = 20:1:600;

T = T';

res_cond = 1.72e-8 * (1 + 3.8e-3*(T-20));

delta = sqrt(2.*res_cond/(mu0.*2.*pi.*f0));

R = (res_cond./(delta.*Rin) + res_cond./(delta.*Rout)) ./ (2.*pi);

alpha_cond = R.*8.686./(2.*Zc0);

res_diel = 2.2e13.*exp(-0.035.*T);

sigma_diel = 1./res_diel;

alpha_total = alpha_total0 + 1e-3.*(T-20);

alpha_diel = (alpha_total - alpha_cond)./8.686;

G = 2.*alpha_diel./Zc;

epsilon_im = G.*log(Rout/Rin)./(4.*pi^2.*f0.*epsilon0);

epsilon_im_cond = 1./(res_diel.*2.*pi.*f0.*epsilon0);

tan_delta = (epsilon_im-epsilon_im_cond)./epsilon_re;

L = mu0./(2*pi).*(log(Rout./Rin)+delta./2.*(1./Rin+1./Rout));

C = 2.*pi.*epsilon0.*epsilon_re./log(Rout./Rin);

Zc = sqrt(L./C);

sigma_diel_T = [T, sigma_diel];
tan_delta_T = [T, tan_delta];
epsilon_im_T = [T, epsilon_im];
epsilon_im_cond_T = [T, epsilon_im_cond];

save('sigma_diel_T.txt', 'sigma_diel_T', '-ascii');
save('tan_delta_T.txt', 'tan_delta_T', '-ascii');
save('epsilon_im_T.txt', 'epsilon_im_T', '-ascii');
save('epsilon_im_cond_T.txt', 'epsilon_im_cond_T', '-ascii');
