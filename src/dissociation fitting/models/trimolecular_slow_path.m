function [p_at_exp,t,p_tot] = trimolecular_slow_path(t_exp,type,densities,indep_vars,k,p0)

mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

kf1 = indep_vars(1);
ka1 = indep_vars(2);
ks1 = indep_vars(3);
kf2 = indep_vars(4);
ka2 = indep_vars(5);
ks2 = indep_vars(6);

p0_1 = log10(k(1));

p0 = [p0_1, 1-p0_1];

t = 0:0.01:max(t_exp);
%analitical solution
p_tot=model_slow_species(t,kf1,ka1,ks1,kf2,ka2,ks2,p0);
p_at_exp=model_slow_species(t_exp,kf1,ka1,ks1,kf2,ka2,ks2,p0);

% %numerical solution
% tspan = 0:0.1:max(t_exp);

% % p0 = [1,0];
% [t,p] = ode23s(@(t,p) transient_model(t,p,koff1,kint,ks,mr,ml,mx,Ac), tspan, p0);
% p_tot = sum(p,2);
% % Pa = 1-exp(-n_tot); 
% p_at_exp = zeros(length(t_exp),1); 
% for i=1:length(t_exp)
%     diff = abs(t-t_exp(i));
%     [tmin,imin] = min(diff);
%     p_at_exp(i) = p_tot(imin);   
% end 

end

function p = model_slow_species(t,kf1,ka1,ks1,kf2,ka2,ks2,p0)%analitical solution

w1 = ka1/(kf1+ka1-ks1);
b1 = kf1+ka1;
pr = w1.*exp(-ks1.*t)+(1-w1).*exp(-b1.*t);

w2 = ka2/(kf2+ka2-ks2);
b2 = kf2+ka2;
px = w2.*exp(-ks2.*t)+(1-w2).*exp(-b2.*t);

p = p0(1).*pr + p0(2).*px; 
end


function dpdt = transient_model(t,p,kf1,ka1,ks1,mr,ml,mx,Ac)%numerical solution
pr = p(1);
ps = p(2);

dpdt = zeros(2,1);
dpdt(1) = (-kf1-ka1).*pr;
dpdt(2) = ka1.*pr -ks1.*ps;

end
