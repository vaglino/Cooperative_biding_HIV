function [L] = cost_fun_association(model,k,t,y_exp,type,densities,indep_vars,options)
% defines problem to be optimized by fmincon

n=length(y_exp);
%evaluate model
[y_fit] = model(t,type,densities,indep_vars,10.^k);
%calculate loss
L = error_fn_association(n,y_fit,y_exp,t,options);

end