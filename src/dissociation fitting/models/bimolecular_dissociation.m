function [p_at_exp,t,p_tot] = bimolecular_dissociation(t_exp,type,densities,indep_vars,k,p0)

mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

koff1 = k(1);

t = 0:0.01:max(t_exp);

p_tot = exp(-koff1.*t);
p_at_exp = exp(-koff1.*t_exp);

end