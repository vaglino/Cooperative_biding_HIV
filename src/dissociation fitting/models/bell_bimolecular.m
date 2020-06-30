function [p_at_exp,t,p_tot] = bell_bimolecular(t_exp,type,densities,indep_vars,k,forces,p0)
load('bell_solution.mat');
n_bins = length(t_exp); %number of force bins to fit

%constant densities
mr=densities(1);
ml=densities(2);
mx=densities(3);
Ac=densities(4);

%parameters to fit
% f_s_0 = k(1);
% k_s_0 = k(2);
% f_f_0 = k(3);
% k_ona_0 = k(4);
% k_offa_0 = k(5);
% % f_ona_0 = k(5);
%parameters to fit
% f_s_0 = k(1);
k_s_0 = k(1);

k_ona_0 = k(2);
k_offa_0 = k(3); k_offa_0 = 15;
f_ona_0 = k(4);
f_offa_0 = k(5);
f_s_0 = k(6);

%fixed parameters
k_f_0 = indep_vars(1);

p_tot = cell(1,n_bins);
p_at_exp = cell(1,n_bins);
t = cell(1,n_bins);

p1 = 1;
p2 = 0;

for i=1:n_bins
    %experimental times for this bin
    t_exp_bin = t_exp{i};
    %time span for each bin
    t_span = 0:0.1:max(t_exp{i});
    t{i} = t_span;
    %force for this specific bin
    f_bin = forces(i);
    
    %calculate two bell curve parameters for force bin
    k_s = k_s_0*exp(f_bin/f_s_0);
%     k_f = k_f_0*exp(f_bin/f_f_0);
%     k_a = k(i+3);
    k_on_a = k_ona_0*exp(f_bin/f_ona_0);
    k_off_a = k_offa_0*exp(-f_bin/f_offa_0);
%     k_on_a = k_ona_0;
%     k_off_a = k_offa_0;
%     k_s = k_s_0;
    k_f = k_f_0;
    
    %analitical solution
%     p_tot{i}=model_slow_species(t_span,k_f,k_a,k_s);
%     p_at_exp{i}=model_slow_species(t_exp_bin,k_f,k_a,k_s);
    
    p_tot{i} = bell_solution(t_span,k_f,k_off_a,k_on_a,k_s,p1,p2);
    p_at_exp{i} = bell_solution(t_exp_bin,k_f,k_off_a,k_on_a,k_s,p1,p2);
end
% t=t_span;
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


% function dpdt = transient_model(t,p,koff1,kint,ks,mr,ml,mx,Ac)%numerical solution
% pr = p(1);
% ps = p(2);
% 
% dpdt = zeros(2,1);
% dpdt(1) = (-koff1-kint).*pr;
% dpdt(2) = kint.*pr -ks.*ps;
% 
% end
