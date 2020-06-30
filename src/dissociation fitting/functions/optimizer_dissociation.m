function [k] = optimizer_dissociation(model,k_best,t,Pa,p0,type,...
                                    densities,indep_vars,n_k,lb,ub,options)

% sets up optimization problem

low = ones(1,n_k).*lb;
up = ones(1,n_k).*ub;

fun = @(k)cost_fun(model,k,t,Pa,p0,type,densities,indep_vars,options);
% options = optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance',1e-8);
opt = optimoptions('fmincon','Display','off');
k = fmincon(fun,k_best,[],[],[],[],low,up,[],opt);

% disp('optimization end')
% toc
end