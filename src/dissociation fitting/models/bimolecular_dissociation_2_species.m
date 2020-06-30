function [p_at_exp,t,p_tot] = bimolecular_dissociation_2_species(t_exp,type,densities,indep_vars,k,p0)

mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

koff1 = k(1);
koff2 = k(2);
% koff1 = indep_vars(1);
% koff2 = indep_vars(2);

t = 0:0.01:max(t_exp);

p_tot = p0(1).*exp(-koff1.*t) + p0(2).*exp(-koff2.*t);
p_at_exp = p0(1).*exp(-koff1.*t_exp) + p0(2).*exp(-koff2.*t_exp);

end