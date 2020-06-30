function [p_at_exp,t,p_tot] = CD4_X4_HIV_dissociation_fixed(t_exp,type,densities,indep_vars,k,p0)

mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

kon3=k(1);
koff3=k(2);
kon4=k(3);
koff4=k(4);

koff1=indep_vars(1); 
koff2=indep_vars(2);
kon3 = indep_vars(3);
kon4 = indep_vars(4);
% p0 = [0.3 0.3 0.4];  %starting distribution of species


tspan = 0:0.1:max(t_exp);

[t,p] = ode23s(@(t,p) transient_model(t,p,kon3,kon4,koff1,koff2,koff3,koff4,mr,ml,mx,Ac), tspan, p0);
p_tot = sum(p,2);
% Pa = 1-exp(-n_tot); 
p_at_exp = zeros(length(t_exp),1); 
for i=1:length(t_exp)
    diff = abs(t-t_exp(i));
    [tmin,imin] = min(diff);
    p_at_exp(i) = p_tot(imin);   
end 
end

function dpdt = transient_model(t,p,kon3,kon4,koff1,koff2,koff3,koff4,mr,ml,mx,Ac)
pr = p(1);
px = p(2);
p3 = p(3);
dpdt = zeros(3,1);
dpdt(1) = -Ac.*mx.*kon4.*pr - koff1.*pr + koff4.*p3;
dpdt(2) = -Ac.*mr.*kon3.*px - koff2.*px + koff3.*p3;
dpdt(3) = Ac.*mx.*kon4.*pr + Ac.*mr.*kon3.*px - koff3.*p3 - koff4.*p3;
end
