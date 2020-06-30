function [p_at_exp,t,p_tot] = trimolecular_dissociation_slow(t_exp,type,densities,indep_vars,k,p0)
% load('tri_diss_solution.mat')
tri_diss_solution = type;

mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

%rates for the bimolecular dissociation pathway and
%slow dissociation pathway are previously obtained from bimolecular
koff1=indep_vars(1);
kint=indep_vars(2);
ks = indep_vars(3);

kon_tri = k(1);
koff_tri = k(2);

tspan = 0:0.01:max(t_exp);

%% analytical sol

p1 = p0(1);
p2 = p0(3);
p_tot = tri_diss_solution(tspan,kint,koff1,koff_tri,kon_tri,ks,mx,p1,p2);
p_at_exp = tri_diss_solution(t_exp,kint,koff1,koff_tri,kon_tri,ks,mx,p1,p2);
t = tspan;

%% numerical sol
%p0 should be [pr ps p3] = []
%ode solver used to be 23s
% [t,p] = ode23s(@(t,p) transient_model(t,p,koff1,kint,ks,kon_tri,koff_tri,mr,ml,mx,Ac), tspan, p0);
% p_tot = sum(p,2);
% % Pa = 1-exp(-n_tot); 
% p_tot = p_tot';
% p_at_exp = zeros(length(t_exp),1); 
% for i=1:length(t_exp)
%     diff = abs(t-t_exp(i));
%     [tmin,imin] = min(diff);
%     p_at_exp(i) = p_tot(imin);   
% end 
end

function dpdt = transient_model(t,p,koff1,kint,ks,kon_tri,koff_tri,mr,ml,mx,Ac)
pr = p(1);
ps = p(2);
p3 = p(3);
dpdt = zeros(3,1);
dpdt(1) = koff_tri.*p3 - kon_tri.*mx.*pr - koff1.*pr - kint.*pr;
dpdt(2) = kint.*pr - ks.*ps;
dpdt(3) = kon_tri.*mx*pr - koff_tri.*p3;
end
