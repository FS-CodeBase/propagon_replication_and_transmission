function [Yn,yn,cn,log_cn,log_Yn] = load_model_solutions(mu,sig)
% % FUNCTION load_model_solutions 
% % Author: Fabian Santiago
% % E-mail: FabianSantiago707@gmail.com
% % Description: loads structured population model solutions with linear
% %              intracellular aggregate dynamics model.
% % 
% % Parameters:
% % CELL DYNAMICS
% % alp - cell division rate
% % bet - cell death rate
% %
% % INTRACELLULAR AGGREGATE DYNAMICS
% % lam - aggregate replication rate
% % gam - Controls the fraction of aggregates that are split between the
% %       mother and daughter cells. gam*Aggregates -> Mother 
% %                                    and (1-gam)*Aggregates-> Daughter
% % 
% % Variables:
% % a(t)   - aggregates at time t
% % da/dt  - aggregate dynamics [ lam*a ] Linear
% % Y(t,a) - aggregate distribution at time t 
% % y(t,a) - normalized aggregate distribution at time t (a PDF) 
% % c(t)   - number of cells at time t
% % U(x)   - initial distribution of aggregates
% % n      - the generation, n = 0, 1, 2, ..., M (maximum generation_
% %
% % PDE System:
% % Y0_t+(da/dt*Y0)_a = -(alp+bet)*Y0                            - 0th gen
% % Yn_t+(da/dt*Yn)_a = -(alp+bet)*Yn+2*alp*gam*Y_[n-0](t,gam*a) - ith gen
% %
% % Note that: Yn(t,a) = cn(t)*yi(t,a)

% Initial Number of Cells
% For ODE model
n00 = 10;

% Population Growth Rate Parameters
alp  = log(2)/(90/60); % Yeast cell division rate (~0.46/hr)
bet  = 0;              % Yeast cell death rate (~0,negligible)

% % INITIALIZE THE DISTRIBUTION OF AGGREGATES IN THE COLONY
% Truncated normal U(x) on (0,inf)
Phi =@(x) 1/2*(1+erf(x/sqrt(2))); phi = @(x) 1/sqrt(2*pi)*exp(-1/2*x.^2);
U   =@(x) 1/sig*phi((x-mu)/sig)/(1-Phi((1-mu)/sig)).*heaviside(x-1);

% Load Pascal's Triangle Factors function
pt_fac = pascals_triangle(25); % 25 rows of Pascal's triangle (can be less) 

% % DEFINE MODEL SOLUTIONS ODEs/PDEs
% Number of cells in generation n at time t
cn =@(n,t) (alp*t).^n./factorial(n).*exp(-(alp+bet)*t).*n00;

% Normalized distribution of aggregate in generation n at time t
yn =@(n,t,a,lam,rho) exp(-lam*t).*...
            sum(pt_fac(n).*gamma_coeffs(n,rho).*...
            U(gamma_coeffs(n,rho).*repmat(a,[n+1,1]).*exp(-lam*t)),1);

% Distribution of aggregate in generation n at time t
Yn =@(n,t,a,lam,rho) cn(n,t).*yn(n,t,a,lam,rho);

% % LOG VERSION OF GENERAL SOLUTION
% Log version of the number of cells in generation n at time t
log_cn =@(n,t) n*log(2*alp*t)-log(factorial(n))+log(n00)-(alp+bet)*t;
% Log version of the initial distribution of aggregates ln(U(x))
log_U  =@(n,t,a,lam,rho) -1/(2*sig^2)*...
                    (rho.^((0:n)-n).*(1-rho).^(-(0:n))...
                                        .*exp(-lam*t).*a-mu).^2;                                 
% Log version of the distribution of aggregate in generation n at time t
log_Yn =@(n,t,a,lam,rho) n*log(alp*t)-log(factorial(n))-(alp+lam)*t... 
                            +log(pt_fac(n)')+((0:n)-n)*log(rho)-...
                                (0:n)*log(1-rho)+log_U(n,t,a,lam,rho)...
                                -1/2*log(2*pi*sig^2);