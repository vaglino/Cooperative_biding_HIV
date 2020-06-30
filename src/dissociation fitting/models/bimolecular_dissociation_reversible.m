function [p_at_exp,t,p_tot] = bimolecular_dissociation_reversible(t_exp,type,densities,indep_vars,k,p0)

mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

kon1 = k(1);
% kon1 = indep_vars(1);
koff1 = k(2);

t = 0:0.01:max(t_exp);

p_tot = (koff1.*exp(-koff1.*t) + Ac.*kon1.*mr.*ml...
        - Ac.*kon1.*mr.*ml.*exp(-koff1.*t))./koff1;
    
p_at_exp = (koff1.*exp(-koff1.*t_exp) + Ac.*kon1.*mr.*ml...
        - Ac.*kon1.*mr.*ml.*exp(-koff1.*t_exp))./koff1;

% p_tot = exp(-koff1.*t);
% p_at_exp = exp(-koff1.*t_exp);

end