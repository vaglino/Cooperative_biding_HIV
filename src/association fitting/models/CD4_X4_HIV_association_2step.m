function [Pa_at_exp,t,Pa,n,n_tot] = CD4_X4_HIV_association_2step(t_exp,type,densities,indep_vars,k)

% CD4_X4_HIV_association_2_step, solution of trimolecular binding model (used for
% fitting CD4-Env-X4 adhesion frequency to Pathway I or II only
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

kon2=k(1);
koff2=k(2);
kon1=indep_vars(1);
koff1=indep_vars(2);

tspan = [0:0.1:10];
% tspan = [0:0.01:10];
n0 = [0 0];
[t,n] = ode23s(@(t,n) transient_model(t,n,kon1,kon2,koff1,koff2,mr,mx,ml,Ac), tspan, n0);
n_tot = sum(n,2);
Pa = 1-exp(-n_tot); 
Pa_at_exp = zeros(length(t_exp),1); 
for i=1:length(t_exp)
    diff = abs(t-t_exp(i));
    [tmin,imin] = min(diff);
    Pa_at_exp(i) = Pa(imin);   
end 
%calculate n at 5s (contact time for bfp experiment) so that we know bond distribution
diff_5s = abs(t-5);
[tmin_5,imin_5] = min(diff_5s);
n_at_5s = n(imin_5,:);
n = n';
end

function dndt = transient_model(t,n,kon1,kon2,koff1,koff2,mr,mx,ml,Ac)
nr = n(1);
n3 = n(2);
dndt = zeros(2,1);
dndt(1) = Ac.*mr*ml.*kon1 - koff1.*nr + koff2.*n3 - Ac.*mx.*kon2.*nr;
dndt(2) = Ac.*mx.*kon2.*nr - koff2.*n3;
end


