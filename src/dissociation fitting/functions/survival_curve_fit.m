function [lnp_fits,ks,avg_forces,boot_stats]...
    = survival_curve_fit(bond_lifetimes,type,indep_vars,densities,p0,k_start,options)

% SURVIVAL_CURVE_FIT
% Performs fitting of given model to experimental data provided.
% an initial set of parameters is determined via hypercube scanning of the
% parameter set. The parameters are then optimized to fit the data, 
% starting with the best hypercube parameters. Both scanning and optimization
% are done in log space to better explore the wide parameter space. Then
% standard deviation of the fitted parameters is calculated by bootstrapping.
% 
% INPUTS:
% mat_bins = cell array made up of survival data for all bins for an
% experimental condition.
% type = model type id (specify ids in pick_model_type function)
% indep_vars = rates constants that we keep fixed
% densities = coating densities from experiment
% p0 = initial probability distribution of bond species
% OUTPUTS:
% lnp_fit = fitted ln of p(t >= t_b)
% ks = fitted rates (n_k by n_bins), each column is rates for one force bin
% avg_forces = mean F for each force bin
% boot_stats = structure with three fields:
%                                     std = standard deviation of ks
%                                     dist = bootstrap k distribution
% n_fit = fitted <n>tot

% Stefano Travaglino, Zhu Lab, 2020/05/16
% -------------------------------------------------------------------------

% densities = [mr,ml,mx,Ac];
n_bins = size(bond_lifetimes,2);
indep_k = indep_vars;

%determine model type and number of free parameters
[model,n_k,lb,ub] = pick_model_type_dissociation(type,options);
if type == 5 % load analytical model for trimolecular dissociation for speed
        type = load('tri_diss_solution.mat');
        type = type.tri_diss_solution;
end

disp('--------------------------------------------------------------------')
disp(char(model))

%initialize outputs
lnp_fits={};
ks = zeros(n_k,n_bins);
avg_forces = zeros(1,n_bins);
k_stds = zeros(n_k,n_bins);
k_stds_clean = zeros(n_k,n_bins);
k_sems = zeros(n_k,n_bins);
k_distributions = cell(1,n_bins);


cm = colormap(jet(n_bins));%define colormap for survival plots
if n_bins==1
    cm = [0 0 1];
end

lowest_lnp = [];
% loop though each bin and fit model
for i=1:n_bins
    % parse data for single force bin
    bin = bond_lifetimes{i};
    forces = bin(:,1);
    lifetimes = bin(:,2);
    lnp = bin(:,4);
    survival_p = exp(lnp);
    avg_force = mean(forces);


    %% survival curve for bin
    a = scatter(lifetimes,lnp,25,cm(i,:),'filled','MarkerFaceAlpha',0.7);
    uistack(a,'bottom')
    pbaspect([1 1 1])
    
    % if independent variables are force dependent, parse indep_vars for force bin 
    if  ~isempty(indep_k) && ~isvector(indep_k) 
        %koff1,kint,kslow
        indep_vars = indep_k(:,i);
%         k = [];
    end
    
    if ~isempty(k_start)
        % use adhesion frequency fit as starting point for optimization
        k_best = log10(k_start);
    else
        % parameter space random scanning by hypercube method 
        [k_best, L] = hypercube_dissociation(model,lifetimes,survival_p,p0,...
             type,densities,indep_vars,n_k,lb,ub,options);
    end

    % optimize from best starting parameter set
    k = optimizer_dissociation(model,k_best,lifetimes,survival_p,p0,...
        type,densities,indep_vars,n_k,lb,ub,options);

    % find std of k by bootstrapping
    [k_sd, k_boot] = bootstrap_dissociation(model,k,lifetimes,...
        survival_p,p0,type,densities,indep_vars,n_k,lb,ub,options);
    
    % Evaluate best fit. Add the 10 so that fit is extended up to 10
    [~,t,survival_p_fit] = model(0:60,type,densities,indep_vars,10.^k,p0);
    
    % plot best fit
    lnp_fit = log(survival_p_fit);
    hold on
    b = plot(t,lnp_fit,'color',cm(i,:),'linewidth',1.5);
    uistack(b,'bottom')

    lowest_lnp = [lowest_lnp min(lnp)];
    
    % save results
    lnp_fit = lnp_fit(1:10:end);
    lnp_fit = lnp_fit(lnp_fit>-7);
    lnp_fits=[lnp_fits {lnp_fit'}];
    
    ks(:,i) = k;
    avg_forces(i) = avg_force;
    k_stds(:,i) = k_sd;
    k_distributions{i} = k_boot;
%     k_boot_clean = rmoutliers(k_boot','gesd')';
    k_boot_clean = rmoutliers(k_boot','thresholdfactor', 6)'; % 6 sigma
    k_sd_clean = std(k_boot_clean,[],2);
    k_stds_clean(:,i) = k_sd_clean; 
    k_sem = std(k_boot,[],2) ./ sqrt(size(k_boot,2));
    k_sems(:,i) = k_sem;
    
   %     Ls(i) = L;
    
    
end

ylim([min(lowest_lnp) 0])
xlim([0 50])
if n_bins ==1
    xlim([0 3])
end



% boot_stats = struct('std',k_stds,'dist',k_distributions);
boot_stats(1).std = k_stds;
boot_stats(1).std_clean = k_stds_clean';
boot_stats(1).sem = k_sems';

boot_stats(1).dist = k_distributions;


pause(0.5)
end

