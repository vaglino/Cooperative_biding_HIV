function [best_k, L] = hypercube_dissociation(model,t,y_exp,y0,type,densities,indep_vars,n_k,lb,ub,options)
%
% HYPERCUBE_ASSOCIATION - hypercube parameter space exploration 
% model kinetic rates exploration by hypercube sampling the parameter space
% ------------------------------------------------------------------------
% INPUTS: 
% model - pass in model function
% k_strt - initial parameter guess (from opt)
% t - time pts
% y_exp - adhesion freq at experimental time points
% y0 - initial p
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


disp('hypercube start')
tic
        
rng('default');
rng(1);
%divide log parameter space in equally spaced intervals
m = length(y_exp);
sampling = options.hypercube_sampling;
% sampling = 30000;

%create parameter vectors
% k =linspace(low_limit,hi_limit,sampling)'*ones(1,n_k);
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

    [y_fit_at_exp,~,~] = model(t,type,densities,indep_vars,rand_k(j,:),y0);
%     [p_at_exp,~,~] = CD4_X4_HIV_dissociation(t,type,densities,indep_vars,rand_k(j,:));

    y_list(j,:) = y_fit_at_exp;

    %find overall error for parameter combination
    L_list(j) = error_fn(m,y_fit_at_exp,y_exp,t,options);
end

%extract best parameters based on index of the lowest error
[L,i_min] = min(L_list);
best_k = rand_k(i_min,:);
best_k = log10(best_k)';
% y = y_list(i_min,:);


% plot_parameter_space(n_k,rand_k,L_list,best_k,L,model)

disp('hypercube end')
toc
end

function plot_parameter_space(n_k,rand_k,L_list,best_k,L,model)


if n_k==1
    figure
    hold off
    scatter(log10(rand_k(:,1)),L_list,30,L_list,'filled')
    hold on
    scatter(best_k(1),L,25,'r','filled')
    title(['parameter space ' char(model)])
    colormap(jet)
elseif n_k==2
    figure
    hold off
    scatter3(log10(rand_k(:,1)),log10(rand_k(:,2)),L_list,30,L_list,'filled')
    hold on
    scatter3(best_k(1),best_k(2),L,25,'r','filled')
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
