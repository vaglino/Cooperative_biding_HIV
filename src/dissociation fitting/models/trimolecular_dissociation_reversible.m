function [p_at_exp,t,p_tot] = trimolecular_dissociation_reversible(t_exp,type,densities,indep_vars,k,p0)

mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

kon1 = k(1);
koff1 = k(2);
kon2 = k(3);
koff2 = k(4);

t = 0:0.01:max(t_exp);

%analytical implementation for speed

p = analytical_sol(kon1,kon2,koff1,koff2,mr,ml,mx,Ac,t,p0);
p_tot = sum(p,1)';

p_exp = analytical_sol(kon1,kon2,koff1,koff2,mr,ml,mx,Ac,t_exp,p0);
p_at_exp = sum(p_exp,1)';

end

function p = analytical_sol(kon1,kon2,koff1,koff2,mr,ml,mx,Ac,t,p0)
p = zeros(2,numel(t));
p(1,:) = (exp(-koff1.*t).*(koff1.*p0(1) - Ac.*kon1.*mr.*ml...
    + Ac.*kon1.*mr.*ml.*exp(koff1.*t)))./koff1;
p(2,:) = (exp(-koff2.*t).*(koff2.*p0(2) - Ac.*kon2.*ml.*mx...
    + Ac.*kon2.*ml.*mx.*exp(koff2.*t)))./koff2;
end