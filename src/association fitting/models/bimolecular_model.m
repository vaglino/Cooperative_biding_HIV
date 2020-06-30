function [Pa_at_exp,t,Pa,n,n_tot] = bimolecular_model(t_exp,type,densities,indep_vars,k)

% BIMOLECULAR MODEL, solution of bimolecular binding model (used for either
% CD4 only or X4 only bonds
% ------------------------------------------------------------------------
% INPUTS: 
% t_exp - query times
% type - empty (in case analytical solution is provided)
% densities - as receptor, ligand, coreceptor, contact area
% k - rate constants
% OUTPUTS:
% Pa_at_exp - adhesion frequency at experimental time pts
% t - time points at which ode is evaluated
% Pa - adhesion frequency at t
% n - <n> at t
% n_tot = total <n> at t

% Stefano Travaglino, ZhuLab, 2020
% ------------------------------------------------------------------------

mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

kon1=k(1);
koff1=k(2);

t=0:0.1:max(t_exp);
% t=0:0.01:max(t_exp);


Pa = 1-exp(-mr*ml*Ac*(kon1/koff1)*(1-exp(-koff1.*t)));
Pa_at_exp = 1-exp(-mr*ml*Ac*(kon1/koff1)*(1-exp(-koff1.*t_exp)));

n = mr*ml*Ac*(kon1/koff1)*(1-exp(-koff1.*t));
% n_at_5s = mr*ml*Ac*(kon1/koff1)*(1-exp(-koff1.*5));
n_tot = n;

end