function [p_at_exp,t,p_tot] = bimolecular_slow_path(t_exp,type,densities,indep_vars,k,p0)

mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

kf = k(1);
ka = k(2);
ks = k(3); 

t = 0:0.01:max(t_exp);
%analitical solution
p_tot=model_slow_species(t,kf,ka,ks);
p_at_exp=model_slow_species(t_exp,kf,ka,ks);

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

function p = model_slow_species(t,kf,ka,ks)%analitical solution

w = ka/(kf+ka-ks);
b = kf+ka;
p = w.*exp(-ks.*t)+(1-w).*exp(-b.*t);
end


function dpdt = transient_model(t,p,koff1,kint,ks,mr,ml,mx,Ac)%numerical solution
pr = p(1);
ps = p(2);

dpdt = zeros(2,1);
dpdt(1) = (-koff1-kint).*pr;
dpdt(2) = kint.*pr -ks.*ps;

end
