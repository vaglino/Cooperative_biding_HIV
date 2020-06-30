function [Pa_fit,t_fit,ks,boot_stats,RSS,n_fit]...
    = adhesion_curve_fit(t,Pa,type,indep_vars,densities,options)

% ADHESION_CURVE_FIT
% Performs fitting of given model to experimental data provided.
% an initial set of parameters is determined via hypercube scanning of the
% parameter set. The parameters are then optimized to fit the data, 
% starting with the best hypercube parameters. Both scanning and optimization
% are done in log space to better explore the wide parameter space. Then
% standard deviation of the fitted parameters is calculated by bootstrapping.
% 
% INPUTS:
% t = contact times
% Pa = experimental Pa
% type = model type id (specify ids in pick_model_type function)
% indep_vars = rates constants that we keep fixed
% densities = coating densities from experiment
% options = hyperparameters for fitting
% OUTPUTS:
% Pa_fit = fitted adhesion frequency
% t_fit = contact time range
% ks = fitted rates
% boot_stats = structure with three fields:
%                                     std = standard deviation of ks
%                                     dist = bootstrap k distribution
%                                     n_ss = bootrstrap s.s <n> dist
%                                     n_boot = bootrstrap <n> dist


% RSS = residual sum of squares
% n_fit = fitted <n>tot

% Stefano Travaglino, Zhu Lab, 2020/05/16
% -------------------------------------------------------------------------


%determine model type and number of free parameters
[model,n_k,lb,ub] = pick_model_type(type,options);

%% adhesion curve fit

pbaspect([1 1 1])

% parameter space random scanning by hypercube method 
[k_best, L_hyper] = hypercube_association(model,t,Pa,type,densities,...
    indep_vars,n_k,lb,ub,options);

% optimize from best starting parameter set
k = optimizer(model,k_best,t,Pa,type,densities,indep_vars,n_k,lb,ub,options);

% find std of k by bootstrapping
% [k_stds, k_boot] = bootstrap(model,k,t,Pa,type,densities,indep_vars,n_k,lb,ub,options);
[k_stds, k_boot] = bootstrap(model,k_best,t,Pa,type,densities,indep_vars,n_k,lb,ub,options);


%Evaluate best fit. Add the 10 so that fit is extended up to 10
[Pa_fit_at_exp,t_fit,Pa_fit,n_fit] = model([t; 10],type,densities,indep_vars,10.^k); 

% calculate fitting error
[L,RSS] = error_fn_association(length(Pa),Pa_fit_at_exp(1:end-1),Pa,t,options);

% std for n at steady state
ns_steady_state = zeros(size(n_fit,1), size(k_boot,2));
n_boot = zeros(size(n_fit,2),size(n_fit,1),size(k_boot,2));
for i=1:size(k_boot,2)
    [~,~,~,n] = model(0:10,type,densities,indep_vars,10.^k_boot(:,i));
    n_ss = n(:,end);
    ns_steady_state(:,i) = n_ss;
    n_boot(:,:,i) = n';
end

% organize outputs
ks=k';
boot_stats(1).std = k_stds;
boot_stats(1).dist = k_boot;
boot_stats(1).n_ss = ns_steady_state;
boot_stats(1).n_boot = n_boot;

end



