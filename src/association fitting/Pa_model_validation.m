%% Master script to fit adhesion frequencies
% Stefano Travaglino - 2019/10/29

% Script fits data for "Cooperative Binding of HIV-1 Envelope Glycoprotein 
% to CD4 and Coreceptor" paper, and in general, provides examples on how to
% use model fitting library

% Fitting of 2-step reversible reaction for HIV:
% (CD4+VLP -> CD4-VLP+CXCR4 <-> CD4-VLP+CXCR4),
% raw data is provided by Dr. Bai Ke, Zhu lab

close all
% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
% hold off
rng('default');
rng(1);


% choose system
% ligand_type = 'Env';
ligand_type = 'gp120';


%% load data
switch ligand_type
    case 'Env'
        num1 = xlsread('Pa_VLP_cd4');
        num2 = xlsread('Pa_VLP_x4');
        num3 = xlsread('Pa_VLP_cd4x4');
        load('options_Env')
        options.B_size = 3;
    case 'gp120'
        num1 = xlsread('Pa_gp120_cd4');
        num2 = xlsread('Pa_gp120_x4');
        num3 = xlsread('Pa_gp120_cd4x4');
        load('options_gp120')
        options.B_size = 3;
        options.hypercube_size = 1000;
%         options.lb = -7;
end

t1_mean = num1(:,1);
Pa1_mean = num1(:,2);
sd1_mean = num1(:,3);
sd1_n = num1(:,4);
t2_mean = num2(:,1);
Pa2_mean = num2(:,2);
sd2_mean = num2(:,3);
sd2_n = num2(:,4);
t3_mean = num3(:,1);
Pa3_mean = num3(:,2);
sd3_mean = num3(:,3);
sd3_n = num3(:,4);

%% fixed parameters
% Ac = 6;%um^2
Ac = 1;%um^2
% switch ligand_type
%     case 'Env'     
%         %             [CD4, Env,  X4, Ac]
%         densities_1 = [279, 7.6,  0,  Ac]; %env
%         densities_2 = [0,   7.6,  18, Ac];
%         densities_3 = [242, 7.6,  18, Ac];
%     case 'gp120'
%         densities_1 = [279, 10.1, 0,  Ac]; %gp120
%         densities_2 = [0,   6.17, 18, Ac];
%         densities_3 = [242, 10.1, 18, Ac];
% end
switch ligand_type
    case 'Env'     
        %             [CD4, Env,  X4, Ac]
        load('densities_env');
        densities_1 = densities_env.dens_1;
        densities_2 = densities_env.dens_2;
        densities_3 = densities_env.dens_3;
    case 'gp120'
        load('densities_gp120');
        densities_1 = densities_gp120.dens_1;
        densities_2 = densities_gp120.dens_2;
        densities_3 = densities_gp120.dens_3;
        options.B_size = 500;
end

%% format data
%format as columns, where the first column is t, second column is Pa.
%and discard NaNs

[t1, Pa1] = clean_association_data(num1);
[t2, Pa2] = clean_association_data(num2);
[t3, Pa3] = clean_association_data(num3);


%% CD4-HIV bimolecular bond fit 

type=0;
[Pa_fit_1,t,k_log_1,boot_1,L_1,n_fit_1] = adhesion_curve_fit(t1,Pa1,...
    type,[],densities_1,options);

figure(1)
hold off
plot(t,Pa_fit_1,'k','linewidth',1.5)
hold on
e = errorbar(t1_mean,Pa1_mean,sd1_mean,'vertical','s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
e.Color = 'k';
e.LineWidth = 1;

%% X4-HIV bimolecular bond fit 

type=0;
[Pa_fit_2,t,k_log_2,boot_2,L_2,n_fit_2] = adhesion_curve_fit(t2,Pa2,...
    type,[],densities_2([3,2,1,4]),options);

figure(1)
hold on
plot(t,Pa_fit_2,'b','linewidth',1.5)
% scatter(t2,Pa2,20,'b','filled')
e = errorbar(t2_mean,Pa2_mean,sd2_mean,'vertical','s','MarkerSize',5,...
    'MarkerEdgeColor','b','MarkerFaceColor','b');
e.Color = 'b';
e.LineWidth = 1;


%% extract biomolecular rates

kon1 = 10.^k_log_1(1);
koff1 = 10.^k_log_1(2);
kon2 = 10.^k_log_2(1);
koff2 = 10.^k_log_2(2);

% for trimolecular fit, keep bimolecular rates fixed
% indep_vars = [kon1; koff1; kon2; koff2]; 

%% CD4-HIV-X4 trimoleclular bond fit (sequential 1-path, 2-step)

indep_vars = [kon1; koff1]; % for CD4-Env + X4

% indep_vars = [kon2; koff2]; % for X4-Env + CD4
% densities_3 = densities_3([3,2,1,4]); % organize as [X4, Env, CD4, Ac]

type=1; %1 for 2 step 1 path model
[Pa_fit_3,t,k_log_3,boot_3,L_3,n_fit_3] = adhesion_curve_fit(t3,Pa3,...
    type,indep_vars,densities_3,options);

figure(1)
plot(t,Pa_fit_3,'r','linewidth',1.5)
% scatter(t3,Pa3,20,'g','filled')
e = errorbar(t3_mean,Pa3_mean,sd3_mean,'vertical','s','MarkerSize',5,...
    'MarkerEdgeColor','r','MarkerFaceColor','r');
e.Color = 'r';
e.LineWidth = 1;


%% format Pa plot
xlabel('t')
ylabel('Pa')
% xlim([0 5
ylim([0 0.6])
legend('cd4','cd4','x4','x4','cd4.x4 (2-step)','cd4.x4',...
    'x4.cd4 (2-step)','cd4.x4 (4-step)','Location','northeastoutside')
legend boxoff
title(ligand_type)
box off
pbaspect([1 1 1])

%% Calculate and plot <n> curves
figure(2)
hold off
% convert Pa to <n> bonds
n1_mean = log(1./(1-Pa1_mean));
n2_mean = log(1./(1-Pa2_mean));
n3_mean = log(1./(1-Pa3_mean));
n1 = sum(n_fit_1,1);
n2 = sum(n_fit_2,1);
n3 = sum(n_fit_3,1);

plot(t,n1,'k','linewidth',1.5)
hold on
e = errorbar(t1_mean,n1_mean,sd1_n,'vertical','s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
e.Color = 'k';
e.LineWidth = 1;

plot(t,n2,'b','linewidth',1.5)
e = errorbar(t2_mean,n2_mean,sd2_n,'vertical','s','MarkerSize',5,...
    'MarkerEdgeColor','b','MarkerFaceColor','b');
e.Color = 'b';
e.LineWidth = 1;

plot(t,n3,'r','linewidth',1.5)
e = errorbar(t3_mean,n3_mean,sd3_n,'vertical','s','MarkerSize',5,...
    'MarkerEdgeColor','r','MarkerFaceColor','r');
e.Color = 'r';
e.LineWidth = 1;

summed = n1 + n2;
plot(t,summed,'g','linewidth',1.5)

xlabel('t')
ylabel('<n>')
% xlim([0 5])
ylim([0 0.8])
legend('cd4','cd4','x4','x4','cd4.x4 (4-step)','cd4.x4',...
    'predicted sum','Location','northeastoutside')
legend boxoff
title(ligand_type)
box off
pbaspect([1 1 1])

predicted = 1-exp(-summed);
predicted = predicted';
figure(1)
plot(t,predicted,'g','linewidth',1.5)


%% clean up rates and output as structures
k_fit_1 = 10.^k_log_1;
k_fit_2 = 10.^k_log_2;
k_fit_3 = 10.^k_log_3;
k_boot_1 = 10.^boot_1.dist;
k_boot_2 = 10.^boot_2.dist;
k_boot_3 = 10.^boot_3.dist;
sd_1 = std(k_boot_1,[],2);
sd_2 = std(k_boot_2,[],2);
sd_3 = std(k_boot_3,[],2);


% bimolecular_C = struct('kon1',k_fit_1(1),'koff1',k_fit_1(2))
% bimolecular_X = struct('kon1',k_fit_2(1),'koff1',k_fit_2(2))
% % trimolecular_2step = struct('kon1',k_fit_1(1),'koff1',k_fit_1(2),...
% %     'kon2',k_fit_3(1),'koff2',k_fit_3(2))
% trimolecular_2path = struct('kon1',k_fit_1(1),'koff1',k_fit_1(2),...
%     'kon2',k_fit_2(1),'koff2',k_fit_2(2),'kon3',k_fit_3(1),...
%     'koff3',k_fit_3(2),'kon4',k_fit_3(3),'koff4',k_fit_3(4))

%% Flux analysis
% k_boot = [k_boot_1;k_boot_2;k_boot_3];
% sd_indep = [sd_1; sd_2];
% sd_ks = sd_3;
% 
% % calculate normalized rates and fluxes at steady state
% type = 2;
% [normalized_rates, fluxes] = analyze_fluxes(type,k_fit_3,indep_vars,...
%                                         densities_3, k_boot, options);
% 
% n_ss_log_1 = log10(n_fit_1(:,end)); % get log of <n>_i at steady state
% n_ss_log_2 = log10(n_fit_2(:,end));
% n_ss_log_3 = log10(n_fit_3(:,end));
% sd_n_ss_1 = std(log10(boot_1.n_ss),[],2); % calculate std of log <n>_i,s.s
% sd_n_ss_2 = std(log10(boot_2.n_ss),[],2);
% sd_n_ss_3 = std(log10(boot_3.n_ss),[],2);
% 
% k_save = [k_fit_1;k_fit_2;k_fit_3];
% % save('rates_env','k_save')
% % save('rates_gp120','k_save')


%% cooperativity analysis

% calculate <n> std from bootstrap distributions of <n> fitted curves
n1_std = std(boot_1.n_boot,[],3);
n2_std = std(boot_2.n_boot,[],3);
n3_boot = sum(boot_3.n_boot,2); % sum of all three bound species
n3_std = std(n3_boot,[],3);

% figure % cooperativity index calculation
% [coop_indx, sd_coop_indx] = cooperativity_index(t,n1',n2',n3',...
%                                             n1_std,n2_std,n3_std);
% 
%                                                                         
% sd_n_fit_3 =  std(boot_3.n_boot,[],3);
% sd_n_cd4 = sd_n_fit_3(:,1);
% sd_n_x4 = sd_n_fit_3(:,2);
% sd_n_tri = sd_n_fit_3(:,3);

% calculate Pa std from bootstrap distributions of <n> fitted curves
% converted from <n> to Pa
Pa1_fit_std = std(n2Pa(boot_1.n_boot),[],3);
Pa2_fit_std = std(n2Pa(boot_2.n_boot),[],3);
n3_boot = sum(boot_3.n_boot,2); % sum of all three bound species
Pa3_fit_std = std(n2Pa(n3_boot),[],3);



function [Pa] = n2Pa(n)
Pa = 1 - exp(-n);
end

function [n] = Pa2n(Pa)
n = log(1 ./ (1-Pa));
end
