function [p_at_exp,t,p_tot] = trimolecular_dissociation_enhancement(t_exp,type,densities,indep_vars,k,p0)

mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

% alpha=log10(k(1));
% beta=log10(k(2));
alpha=k(1);
beta=k(2);
% alpha=0;
% beta=1;

kon1=indep_vars(1);
koff1=indep_vars(2);
kon2=indep_vars(3);
koff2=indep_vars(4);
% kon1=k(1);
% koff1=k(2);
% kon2=k(3);
% koff2=k(4);

% p0 = [0.5 0.5];

tspan = 0:0.01:max(t_exp);

[t,p] = ode23s(@(t,p) transient_model(t,p,kon1,kon2,koff1,koff2,alpha,beta,mr,ml,mx,Ac), tspan, p0);
p_tot = sum(p,2);
p_tot = p_tot';
% Pa = 1-exp(-n_tot); 
p_at_exp = zeros(length(t_exp),1); 
for i=1:length(t_exp)
    diff = abs(t-t_exp(i));
    [tmin,imin] = min(diff);
    p_at_exp(i) = p_tot(imin);   
end 


% t = 0:0.01:max(t_exp);
% p = analytical_sol(kon1,kon2,koff1,koff2,alpha,beta,mr,ml,mx,Ac,t,p0);
% p_tot = sum(p,1)';
% 
% p_exp = analytical_sol(kon1,kon2,koff1,koff2,alpha,beta,mr,ml,mx,Ac,t_exp,p0);
% p_at_exp = sum(p_exp,1)';



end

function dpdt = transient_model(t,p,kon1,kon2,koff1,koff2,alpha,beta,mr,ml,mx,Ac)
pr = p(1);
px = p(2);
% p3 = p(3);
dpdt = zeros(2,1);
dpdt(1) = mr.*ml.*Ac.*kon1*(1+alpha*px) - koff1.*pr;
dpdt(2) = mx.*ml.*Ac.*kon2.*(1+beta*pr) - koff2.*px;

end

function p = analytical_sol(kon1,kon2,koff1,koff2,alpha,beta,mr,ml,mx,Ac,t,p0)
p = zeros(2,numel(t));
p(1,:) = (exp(-koff1.*t).*(koff1.*p0(1) - alpha.*Ac.*kon1.*mr.*ml...
    + alpha.*Ac.*kon1.*mr.*ml.*exp(koff1.*t)))./koff1;
p(2,:) = (exp(-koff2.*t).*(koff2.*p0(2) - beta.*Ac.*kon2.*ml.*mx...
    + beta.*Ac.*kon2.*ml.*mx.*exp(koff2.*t)))./koff2;
end