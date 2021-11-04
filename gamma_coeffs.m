function gam_vec = gamma_coeffs(n,gam)
% % FUNCTION: gamma_coeffs
% % Author: Fabian Santiago 
% % Description: 
% %     Generates the coefficients (1/gamm)^n*(1-(1/gamma))^(n-i).
I1 = flip(0:n);
I2 = 0:n; 
P1 = 1/gam;
P2 = 1/(1-gam);

gam_vec = ((P1.^I1).*(P2.^I2))';


