function [Pa_at_exp,t,Pa,n,n_tot] = CD4_X4_HIV_association(t_exp,type,densities,indep_vars,k)

% CD4_X4_HIV_association, solution of trimolecular binding model (used for
% fitting CD4-Env-X4 adhesion frequency
% ------------------------------------------------------------------------
% INPUTS: 
% t_exp - query times
% type - empty (in case analytical solution is provided)
% densities - as receptor, ligand, coreceptor, contact area
% indep_vars = constant rates
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

kon3=k(1);
koff3=k(2);
kon4=k(3);
koff4=k(4);

kon1=indep_vars(1);
koff1=indep_vars(2);
kon2=indep_vars(3);
koff2=indep_vars(4);

p0 = [0 0 0];


tspan = [0:0.1:max(t_exp)];
% tspan = [0:0.01:max(t_exp)];
% p0 = [0.8 0 0.2]; %starting distribution of species
[t,n] = ode23s(@(t,n) transient_model(t,n,kon1,kon2,kon3,kon4,koff1,koff2,koff3,koff4,mr,ml,mx,Ac), tspan, p0);
n_tot = sum(n,2);
Pa = 1-exp(-n_tot); 
Pa_at_exp = zeros(length(t_exp),1); 
for i=1:length(t_exp)
    diff = abs(t-t_exp(i));
    [tmin,imin] = min(diff);
    Pa_at_exp(i) = Pa(imin);   
end 
n = n';
end


function dndt = transient_model(t,n,kon1,kon2,kon3,kon4,koff1,koff2,koff3,koff4,mr,ml,mx,Ac)
nr = n(1);
nx = n(2);
n3 = n(3);
dndt = zeros(3,1);
dndt(1) = -mx.*kon4.*nr - koff1.*nr + koff4.*n3 + mr.*ml.*Ac.*kon1;
dndt(2) = -mr.*kon3.*nx - koff2.*nx + koff3.*n3 + mx.*ml.*Ac.*kon2;
dndt(3) = mx.*kon4.*nr + mr.*kon3.*nx - koff3.*n3 - koff4.*n3;
end
