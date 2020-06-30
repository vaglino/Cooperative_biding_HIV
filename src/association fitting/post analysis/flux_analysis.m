function [normalized_rates, normalized_fluxes]...
    = flux_analysis(type,ks,indep_vars,dens, k_boot, options)

% FLUX_ANALYSIS
% given the rate constants obtained from optimization, calculates steady
% state rates and net rates (fluxes) and their standard deviations


model = pick_model_type(type,options);

dens = dens';

[~,~,~,n] = model(0:20,type,dens,indep_vars,ks);
n_ss = n(:,end); % <n>i at steady state

%calculate steady state fluxes for all densities combinations
fluxes = rate2flux(ks,indep_vars,dens,n_ss);
fluxes_log = log10(fluxes);
flux_boot = rate2flux(k_boot,[],dens,n_ss);

% std rates
std_fluxes = std(flux_boot, [], 2);
% flux_boot = flux_boot_clean;

% std log rates
flux_boot_log = log10(flux_boot);
std_log_fluxes = std(flux_boot_log, [], 2);

figure % plot
subplot(1,3,1)
errorbar(fluxes_log,std_log_fluxes); pbaspect([1 1 1])
subplot(1,3,2)
boxplot(flux_boot_log'); pbaspect([1 1 1])
subplot(1,3,3)
plot(flux_boot_log); pbaspect([1 1 1])

% remove outlier rates
% flux_boot_clean = rmoutliers(log10(flux_boot)','grubbs','ThresholdFactor',0.1);
% flux_boot_clean = log10(rmoutliers(flux_boot,2,'ThresholdFactor',2))';%, gp120 
flux_boot_clean = log10(rmoutliers(flux_boot,2,'ThresholdFactor',1.5))';%,'ThresholdFactor',0.01); 


flux_boot_log = flux_boot_clean';
flux_boot = 10.^flux_boot_log;


%% find CI for differeces between rates

% propagate_uncertainty_sum(std_fluxes([1,3,5,7]),std_fluxes([2,4,6,8]))



flux_diffs = [flux_boot(1,:)-flux_boot(2,:);...
            flux_boot(3,:)-flux_boot(4,:);...
            flux_boot(5,:)-flux_boot(6,:);...
            flux_boot(7,:)-flux_boot(8,:)];

net_fluxes = fluxes(1:2:7,:)-fluxes(2:2:8,:);
std_net_fluxes = std(flux_diffs,[],2);

figure


CI_error = confidence_interval(flux_diffs, 0.95);
% errorbar(net_fluxes,CI_error);
errorbar(net_fluxes,std_net_fluxes);
% CI_bounds = [net_fluxes-CI_error, net_fluxes+CI_error];
xlabel('fluxes')
ylabel('s^-1')

normalized_rates(1).rates = fluxes_log;
normalized_rates(1).error = std_log_fluxes;
normalized_fluxes(1).fluxes = net_fluxes;
normalized_fluxes(1).error = std_net_fluxes;

end




%% function that converts rates to fluxes, normalizing by densities and n
function [fluxes] = rate2flux(k,indep_vars,dens,ns_ss)
%INPUTS: 
% k = fitted rates
% indep_vars = constant rates
% dens = matrix of all densities combinations, or single vector of densities
%          [mr_1 - - - - mr_n
%           ml_1 - - - - ml_n
%           mx_1 - - - - mx_n];
% ns_ss = steady state <n> of all species i.e. 
%           [nr_1 - - - nr_n
%            nx_1 - - - nx_n
%            nt_1 - - - nt_n]
%OUTPUTS:
%normalized fluxes

% Ac = 1;
% n_combs = size(dens,2);
rates = [indep_vars; k];

%vectorized code for speed
mr = dens(1,:);
ml = dens(2,:);
mx = dens(3,:);
Ac = dens(4,:);
nr = ns_ss(1,:);
nx = ns_ss(2,:); 
nt = ns_ss(3,:); 

multiplier = [Ac.*mr.*ml;  % Ac*mcd4*menv*kon1
            nr;            %      <n>cd4*koff1
            Ac.*mx.*ml;    %  Ac*mx4*menv*kon2
            nx;            %       <n>x4*koff2
            mr.*nx;        %   mcd4*<n>x4*kon3
            nt;            %      <n>tri*koff3
            mx.*nr;        %   mx4*<n>cd4*kon4
            nt];           %      <n>tri*koff4
        
fluxes = rates.*multiplier;

end

% function outlier_detection()


