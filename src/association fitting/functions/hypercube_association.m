function [best_k, L] = hypercube_association(model,t,y_obs,type,densities,indep_vars,n_k,lb,ub,options)

% HYPERCUBE_ASSOCIATION - hypercube parameter space exploration 
% ------------------------------------------------------------------------
% INPUTS: 
% model - pass in model function
% k_strt - initial parameter guess (from opt)
% t - time pts
% y_obs - adhesion freq at experimental time points
% type - empty, or analytical solution
% densities - site densities
% indep_vars - constant rates
% n_k = number of free rates
% lb, ub - lower and upper boundary for optimization
% options - struct with hyperparameter settings
% OUTPUTS:
% best_k - best set of parameter ks
% L - loss

% Stefano Travaglino, ZhuLab, 2020
% ------------------------------------------------------------------------
rng('default');
rng(1);
%divide log parameter space in equally spaced intervals
m = length(y_obs);

sampling = options.hypercube_sampling;
% sampling = 10000;

%create parameter vectors
k = linspace(0,1,sampling)' * ones(1,n_k);
int_size = ub-lb;
k = k .* int_size + lb;

perm_mat = [];
for i=1:n_k
    perm_mat = [perm_mat, randperm(size(k,1))'];
end
%randomize each parameter vector
rand_k = k(perm_mat);
rand_k = 10.^rand_k;

y_list = zeros(sampling,m);%keep track of p of bonds and their associated errors
L_list = zeros(sampling,1);

for j=1:sampling
    %Find y for each hypercube coeff combination

    [p_at_exp,~,~] = model(t,type,densities,indep_vars,rand_k(j,:));

    y_list(j,:) = p_at_exp;

    %find overall error for parameter combination
    L_list(j) = error_fn_association(m,y_list(j,:)',y_obs,t,options);
end


%extract best parameters based on index of the lowest error
[L,i_min] = min(L_list);
best_k = rand_k(i_min,:);
best_k = log10(best_k);
% y = y_list(i_min,:);

plot_parameter_space(n_k,rand_k,L_list,best_k,L,model)

end

function plot_parameter_space(n_k,rand_k,L_list,best_k,L,model)


if n_k==1
    figure
    hold off
    scatter(log10(rand_k(:,1)),L_list,30,L_list,'filled')
    hold on
    scatter(best_k(1),L,50,'r','filled')
    title(['parameter space ' char(model)])
    colormap(jet)
elseif n_k==2
    figure
    hold off
    scatter3(log10(rand_k(:,1)),log10(rand_k(:,2)),L_list,30,L_list,'filled')
    hold on
    scatter3(best_k(1),best_k(2),L,50,'r','filled')
    title(['parameter space ' char(model)])
    colormap(jet)
elseif n_k>=3
    figure
    hold off
    scatter3(log10(rand_k(:,1)),log10(rand_k(:,2)),log10(rand_k(:,3)),25,L_list,'filled','markerfacealpha',0.7)
    hold on
    scatter3(best_k(1),best_k(2),best_k(3),50,'r','filled')
    title(['parameter space ' char(model)])
    xlabel('p1')
    ylabel('p2')
    zlabel('p3')
    colormap(jet)
end
end
