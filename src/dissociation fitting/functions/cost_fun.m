function [L] = cost_fun(model,k,t,y_exp,p0,type,densities,indep_vars,options)
% defines problem to be optimized by fmincon

n=length(y_exp);
% evaluate model
[y_fit] = model(t,type,densities,indep_vars,10.^k,p0);

% calculate loss
L = error_fn(n,y_fit,y_exp,t,options);

end