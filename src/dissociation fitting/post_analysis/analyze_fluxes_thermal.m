function [normalized_rates,normalized_fluxes]...
    = analyze_fluxes_thermal(type,ks,indep_vars,dens,p0,k_boot)

mr = dens(1);
ml = dens(2);
mx = dens(3);
Ac = dens(4);

dens_comb = combvec(mr,ml,mx); %creates a 3 by n1*n2*n3 matrix of all 
n_combs = size(dens_comb,2);
dens_comb = [dens_comb; Ac.*ones(1,n_combs)];

ns_ss_save = p0;

%calculate steady state fluxes for all densities combinations
fluxes = rate2flux(ks,indep_vars,dens_comb,ns_ss_save);
fluxes_log = log10(fluxes);
% std_fluxes = rate2flux(sd_ks,sd_indep,dens_comb,ns_ss_save);
flux_boot = rate2flux(k_boot,[],dens_comb,ns_ss_save);
flux_boot_log = log10(flux_boot);
std_fluxes = std(flux_boot_log, [], 2);


figure
subplot(1,3,1)
errorbar(fluxes_log,std_fluxes);
subplot(1,3,2)
boxplot(flux_boot_log')
subplot(1,3,3)
plot(flux_boot_log)


%% find error for differeces between rates

% % b_length = size(k_boot_1,2);
figure
flux_diffs = [-flux_boot(1,:);...
            -flux_boot(2,:);...
            flux_boot(3,:)-flux_boot(4,:);...
            flux_boot(5,:)-flux_boot(6,:)];

net_fluxes = [0;0;fluxes([3,5])] - fluxes([1,2,4,6]);
CI_error = confidence_interval(flux_diffs, 0.95);
errorbar(net_fluxes,CI_error);
% CI_bounds = [net_fluxes-CI_error, net_fluxes+CI_error];

normalized_rates(1).rates = fluxes_log;
normalized_rates(1).error = std_fluxes;
normalized_fluxes(1).fluxes = net_fluxes;
normalized_fluxes(1).error = CI_error;

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

multiplier = [nr;          %      <n>cd4*koff1
            nx;            %       <n>x4*koff2
            mr.*nx;        %   mcd4*<n>x4*kon3
            nt;            %      <n>tri*koff3
            mx.*nr;        %   mx4*<n>cd4*kon4
            nt];           %      <n>tri*koff4
        
fluxes = rates.*multiplier;


end


