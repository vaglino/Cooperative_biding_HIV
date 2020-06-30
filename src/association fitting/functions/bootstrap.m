function [k_sd, ks] = bootstrap(model,k_strt,t,Pa,type,densities,indep_vars,n_k,lb,ub,options)

% BOOTSTRAP - bootstrap algorithm 
% ------------------------------------------------------------------------
% INPUTS: 
% model - pass in model function
% k_strt - initial parameter guess (from opt)
% t - time pts
% Pa - adhesion freq at t
% type - empty, or analytical solution
% densities - site densities
% indep_vars - constant rates
% n_k = number of free rates
% lb, ub - lower and upper boundary for optimization
% options - struct with hyperparameter settings
% OUTPUTS:
% k_sd - standard deviation of rates
% ks - all rates obtained from each bootstrap sample

% Stefano Travaglino, ZhuLab, 2020
% ------------------------------------------------------------------------

B_size = options.B_size; % bootstrap resampling size 
% B_size = 300; 
n=length(t); % n number of observations

% resampling index matrix,
% random index between 1 and n, each column contains the randomize index
% for a single bootstrap resampling 
i_rsmpl = randi(n,n,B_size);
% resample
t_rsmpl = t(i_rsmpl);
Pa_rsmpl = Pa(i_rsmpl);

% calculate kinetic rates for each resample
ks = zeros(length(k_strt), B_size); % save all rates
for i=1:B_size
    
    k = optimizer(model,k_strt,t_rsmpl(:,i),Pa_rsmpl(:,i),type,densities,indep_vars,n_k,lb,ub,options);
    
    ks(:,i) = k'; % save new set of rates
end

k_sd = std(ks,[],2);

end