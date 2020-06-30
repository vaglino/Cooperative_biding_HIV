%% Master script to fit dissociations
% Stefano Travaglino - 2019/10/29

% Script fits data for "Cooperative Binding of HIV-1 Envelope Glycoprotein 
% to CD4 and Coreceptor" paper, and in general, provides examples on how to
% use model fitting library to fit survival probability curves
% raw data is provided by Dr. Bai Ke, Zhu lab

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
hold off
close all
rng('default');
rng(1);

% Define site densities:
%              [CD4,      Env,     X4,   Ac]
densities_v = [35.114045, 3.17, 58.76262, 1]; % Env
densities_g = [35.114045, 5.01, 58.76262, 1]; % rgp120

% load adhesion frequency results
ks_v = load('rates_env.mat');   ks_v = ks_v.k_save;
ks_g = load('rates_gp120.mat'); ks_g = ks_g.k_save;
opt_Env = load('options_Env');      opt_Env = opt_Env.options;
opt_gp120 = load('options_gp120');  opt_gp120 = opt_gp120.options;

type = 2; % two-path, two-step model
n_v_norm = initial_distribution(type,ks_v,densities_v,'trimolecular',opt_Env);
n_g_norm = initial_distribution(type,ks_g,densities_g,'trimolecular',opt_gp120);

p0_v_norm = n_v_norm; % [ncd4, nx4, ntri]
p0_g_norm = n_g_norm;


%% load data
%VLP
sheet_name = 'bond_lifetime_distributions_v2.xlsx';
data_vlp_cd4_x4 = readmatrix(sheet_name,'Sheet','VLP vs CD4.X4 Jurkat');
data_vlp_cd4 = readmatrix(sheet_name,'Sheet','VLP vs CD4 Jurkat');
data_vlp_x4 = readmatrix(sheet_name,'Sheet','VLP vs CD4.X4 Jurkat aCD4');
data_vlp_cd4_x4_c52l = readmatrix(sheet_name,'Sheet','VLP vs CD4.X4 Jurkat w C52L');
%gp120
data_gp120_cd4_x4 = readmatrix(sheet_name,'Sheet','gp120 vs CD4.X4 Jurkat');
data_gp120_cd4 = readmatrix(sheet_name,'Sheet','gp120 vs CD4 Jurkat');
data_gp120_x4 = readmatrix(sheet_name,'Sheet','gp120 vs CD4.X4 Jurkat aCD4');
%VLP thermal fluctuation
sheet_name_survival = 'thermal_survival_HIV.xlsx';
data_vcx_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct VLP vs CD4.X4.J');
data_vc_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct VLP vs CD4.J');
data_vx_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct VLP vs X4.J');
%gp120 thermal fluctuation
data_gcx_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct gp120 vs CD4.X4.J');
data_gc_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct gp120 vs CD4.J');
data_gx_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct gp120 vs X4.J');

%for each experimental condition clean data (of zeros) and organize bins in
%cell matrix
vcx = clean_data(data_vlp_cd4_x4);
vc = clean_data(data_vlp_cd4);
vx = clean_data(data_vlp_x4);
vcx_c52l = clean_data(data_vlp_cd4_x4_c52l);
gcx = clean_data(data_gp120_cd4_x4);
gc = clean_data(data_gp120_cd4);
gx = clean_data(data_gp120_x4);
vcx_thermal = {data_vcx_thermal};
vc_thermal = {data_vc_thermal};
vx_thermal = {data_vx_thermal};
gcx_thermal = {data_gcx_thermal};
gc_thermal = {data_gc_thermal};
gx_thermal = {data_gx_thermal};


%% fit thermal fluctuation data
opt_thrml = load('options_thrml.mat'); opt_thrml = opt_thrml.options;

% Env
figure(1)
subplot(2,3,1)
type=0;%normally 0
[lnp_vc_thermal,k_vc_thermal,~,vc_thrml_boot]... 
    = survival_curve_fit(vc_thermal,type,[],densities_v,[],[],opt_thrml);
title('vlp-cd4')

subplot(2,3,2)
type=0;
[lnp_vx_thermal,k_vx_thermal,~,vx_thrml_boot]...
    = survival_curve_fit(vx_thermal,type,[],densities_v,[],[],opt_thrml);
title('vlp-x4')

indep_vars_v_thrml = 10.^[k_vc_thermal; k_vx_thermal];

subplot(2,3,3)
type = 4;
[lnp_vcx_thermal,k_vcx_thermal,~,vcx_thrml_boot] = survival_curve_fit(vcx_thermal,...
    type,indep_vars_v_thrml,densities_v,p0_v_norm,ks_v(5:end),opt_thrml);
title('vlp-cd4-x4')

%--------------------------------------------------------------------------
% rgp120
subplot(2,3,4)
type=0;
[lnp_gc_thermal,k_gc_thermal,~,gc_thrml_boot]...
    = survival_curve_fit(gc_thermal,type,[],densities_g,[],[],opt_thrml);
title('gp120-cd4')

subplot(2,3,5)
type=0;
[lnp_gx_thermal,k_gx_thermal,~,gx_thrml_boot] ...
    = survival_curve_fit(gx_thermal,type,[],densities_g,[],[],opt_thrml);
title('gp120-x4')

indep_vars_g_thrml = 10.^[k_gc_thermal; k_gx_thermal];
subplot(2,3,6)
type=4;
[lnp_gcx_thermal,k_gcx_thermal,~,gcx_thrml_boot] = survival_curve_fit(gcx_thermal,...
    type,indep_vars_g_thrml,densities_g,p0_g_norm,ks_g(5:end),opt_thrml);
title('gp120-cd4-x4')

%% clean up thermal fluctuation rates

k_vc_thrml = 10.^k_vc_thermal;
k_vx_thrml = 10.^k_vx_thermal;
k_vcx_thrml = 10.^k_vcx_thermal;
k_gc_thrml = 10.^k_gc_thermal;
k_gx_thrml = 10.^k_gx_thermal;
k_gcx_thrml = 10.^k_gcx_thermal;

k_boot_v = 10.^[vc_thrml_boot.dist{1};
                vx_thrml_boot.dist{1};
                vcx_thrml_boot.dist{1}];
k_boot_g = 10.^[gc_thrml_boot.dist{1};
                gx_thrml_boot.dist{1};
                gcx_thrml_boot.dist{1}];
            
[norm_rates_v, fluxes_v] = analyze_fluxes_thermal([],k_vcx_thrml,...
                indep_vars_v_thrml,densities_v,p0_v_norm,k_boot_v);
[norm_rates_g, fluxes_g] = analyze_fluxes_thermal([],k_gcx_thrml,...
                indep_vars_g_thrml,densities_g,p0_g_norm,k_boot_g);

% output fitted ks
trimolecular_v_thermal = struct('koff1',k_vc_thrml,'koff2',k_vx_thrml,...
    'kon3',k_vcx_thrml(1),'koff3',k_vcx_thrml(2),...
    'kon4',k_vcx_thrml(3),'koff4',k_vcx_thrml(4))
trimolecular_g_thermal = struct('koff1',k_gc_thrml,'koff2',k_gx_thrml,...
    'kon3',k_gcx_thrml(1),'koff3',k_gcx_thrml(2),...
    'kon4',k_gcx_thrml(3),'koff4',k_gcx_thrml(4))

close(2:5) % close rate figures
%% fit lifetime distribution data (force dependent) Env

opt_force = load('options_force.mat'); opt_force = opt_force.options;
% opt_force.B_size = 1000;
% opt_force.hypercube_sampling = 50000;

% remember to change error function
figure(2)
%for each experimental condition fit each bin and plot 
subplot(2,3,1) %bimolecular VLP-CD4
type=1;
[lnp_vc,k_vc,fs_vc,vc_boot] = survival_curve_fit(vc,type,[],densities_v,[],[],opt_force);
% type = 4; %uncomment for bell model fitting
% indep_vars_vc = 8.82;
% [lnp_vc,k_vc,fs_vc,std_vc,k_dist_vc] = survival_curve_fit_v3(vc,type,indep_vars_vc,densities_v,[]);
title('Env-CD4')

figure(3)
fs_vcx = extract_avg_forces(vcx);
indep_vars_vc = 10.^interpolate_rates(fs_vc,k_vc,fs_vcx,'poly','Env-CD4',true);

figure(2)
subplot(2,3,2) %bimolecular VLP-X4
type=1;
[lnp_vx,k_vx,fs_vx,vx_boot] = survival_curve_fit(vx,type,[],densities_v,[],[],opt_force);
title('Env-X4')

figure(4)
% subplot(1,2,1)
indep_vars_vx = 10.^interpolate_rates(fs_vx,k_vx,fs_vcx,'poly','Env-X4',true);

 %trimolecular VLP-CD4-X4

p0_vlp = [p0_v_norm(2) 0 p0_v_norm(3)];
p0_vlp_norm = p0_vlp./sum(p0_vlp);

figure(2)
subplot(2,3,3)
% indep_vars_v = [indep_vars_vc; indep_vars_vx];
indep_vars_v = [indep_vars_vx];

type = 5;%5
tic
[lnp_vcx,k_vcx,fs_vcx,vcx_boot] = survival_curve_fit(vcx,type,indep_vars_v,...
                        densities_v([3 2 1 4]),p0_vlp_norm,[],opt_force);

toc
title('Env-CD4-X4')
figure(5)
interpolate_rates(fs_vcx,k_vcx,[],'poly','Env-CD4.X4',true);

%% trimolelcular with C52L treatment

figure(6)
fs_vcx_c52l = extract_avg_forces(vcx_c52l);

type = 5;%5

indep_vars_v_c52l = 10.^interpolate_rates(fs_vx,k_vx,fs_vcx_c52l,'poly','Env-X4',false);
[lnp_vcx_c52l,k_vcx_c52l,fs_vcx_c52l,vcx_c52l_boot] = survival_curve_fit(vcx_c52l,...
    type,indep_vars_v_c52l,densities_v([3 2 1 4]),p0_vlp_norm,[],opt_force);

% indep_vars_v_c52l = 10.^interpolate_rates(fs_vc,k_vc,fs_vcx_c52l,'poly','Env-cd4',false);
% p0_vlp_c52l = [p0_v_norm(1) 0 p0_v_norm(3)];
% p0_vlp_c52l_norm = p0_vlp_c52l./sum(p0_vlp_c52l);
% [lnp_vcx_c52l,k_vcx_c52l,fs_vcx_c52l,vcx_c52l_boot] = survival_curve_fit(vcx_c52l,...
%     type,indep_vars_v_c52l,densities_v,p0_vlp_c52l_norm,[],opt_force);



%% fit lifetime distribution data (force dependent) rgp120
figure(2)
subplot(2,3,4) %bimolecular rgp120-CD4
type=1;
[lnp_gc,k_gc,fs_gc,gc_boot] = survival_curve_fit(gc,type,[],densities_g,[],[],opt_force);
title('rgp120-CD4')

fs_gcx = extract_avg_forces(gcx);
figure(7)
indep_vars_gc = 10.^interpolate_rates(fs_gc,k_gc,fs_gcx,'poly','rgp120-CD4',true);

figure(2)
subplot(2,3,5) %bimolecular rgp120-X4
type=1;
[lnp_gx,k_gx,fs_gx,gx_boot] = survival_curve_fit(gx,type,[],densities_g,[],[],opt_force);
title('rgp120-X4')

figure(8)
indep_vars_gx = 10.^interpolate_rates(fs_gx,k_gx,fs_gcx,'poly','rgp120-X4',true);

%trimolecular rgp120-CD4-X4
p0_gp120 = [p0_g_norm(1) 0 p0_g_norm(3)];
p0_gp120_norm = p0_gp120./sum(p0_gp120);


figure(2)
subplot(2,3,6)
% indep_vars_g = [indep_vars_gc; indep_vars_gx];
indep_vars_g = [indep_vars_gc];

type=5;
[lnp_gcx,k_gcx,fs_gcx,gcx_boot] = survival_curve_fit(gcx,type,indep_vars_g,...
                                densities_g,p0_gp120_norm,[],opt_force);
title('rgp120-CD4-X4')
figure(9)
interpolate_rates(fs_gcx,k_gcx,[],'poly','rgp120-CD4.X4',true);


%% Catch bonds analysis

% CD4 with Env or gp120 force lifetime curve
figure(10)
subplot(2,2,1); pbaspect([1 1 1]);
[curve_vc, means_exp_vc] = catch_bond_analysis(vc,1,k_vc,[],[],densities_v,[],opt_force);
[curve_gc, means_exp_gc] = catch_bond_analysis(gc,1,k_gc,[],[],densities_g,[],opt_force);

[means_fit_vc, ~] = force_lifetime_curve(lnp_vc,fs_vc,vc);
[means_fit_gc, ~] = force_lifetime_curve(lnp_gc,fs_gc,gc);

% X4 with Env or gp120 force lifetime curve
subplot(2,2,2); pbaspect([1 1 1]);
[curve_vx, means_exp_vx] = catch_bond_analysis(vx,1,k_vx,[],[],densities_v,[],opt_force);
[curve_gx, means_exp_gx] = catch_bond_analysis(gx,1,k_gx,[],[],densities_g,[],opt_force);

[means_fit_vx, ~] = force_lifetime_curve(lnp_vx,fs_vx,vx);
[means_fit_gx, ~] = force_lifetime_curve(lnp_gx,fs_gx,gx);

% CD4.X4 with Env or gp120 force lifetime curve
subplot(2,2,3); pbaspect([1 1 1]);
[curve_vcx, means_exp_vcx] = catch_bond_analysis(vcx,5,k_vcx,k_vx,fs_vx,densities_v([3 2 1 4]),p0_vlp_norm,opt_force);
[curve_gcx, means_exp_gcx] = catch_bond_analysis(gcx,5,k_gcx,k_gc,fs_gc,densities_g,p0_gp120_norm,opt_force);

[means_fit_vcx, ~] = force_lifetime_curve(lnp_vcx,fs_vcx,vcx);
[means_fit_gcx, curve_gcx] = force_lifetime_curve(lnp_gcx,fs_gcx,gcx);

% CD4.X4 with Env and C52L treatment
subplot(2,2,4); pbaspect([1 1 1]);
[means_fit_vcx_c52l, ~] = force_lifetime_curve(lnp_vcx_c52l,fs_vcx_c52l,vcx_c52l);
[curve_vcx_c52l, means_exp_vcx_c52l] = catch_bond_analysis(vcx_c52l,5,k_vcx_c52l,k_vx,fs_vx,densities_v([3 2 1 4]),p0_vlp_norm,opt_force);
% [curve_vcx_c52l, means_exp_vcx_c52l] = catch_bond_analysis(vcx_c52l,5,k_vcx_c52l,k_vc,fs_vc,densities_v,p0_vlp_c52l_norm,opt_force);

% means of thermal fluctuations
subplot(2,2,1)
[means_fit_vc_thrml, ~] = force_lifetime_curve(lnp_vc_thermal,0,vc_thermal);
[means_fit_gc_thrml, ~] = force_lifetime_curve(lnp_gc_thermal,0,gc_thermal);
subplot(2,2,2)
[means_fit_vx_thrml, ~] = force_lifetime_curve(lnp_vx_thermal,0,vx_thermal);
[means_fit_gx_thrml, ~] = force_lifetime_curve(lnp_gx_thermal,0,gx_thermal);
subplot(2,2,3)
[means_fit_vcx_thrml, ~] = force_lifetime_curve(lnp_vcx_thermal,0,vcx_thermal);
[means_fit_gcx_thrml, ~] = force_lifetime_curve(lnp_gcx_thermal,0,gcx_thermal);


%% time-cooperativity analysis

% calculate initial number of bonds of bimolecular interactions
type = 0; % bimolecular association
[~, n_ss_vc] = initial_distribution(type,ks_v(1:2),densities_v,'bimolecular',opt_Env);
[~, n_ss_vx] = initial_distribution(type,ks_v(3:4),densities_v([3,2,1,4]),'bimolecular',opt_Env);
[~, n_ss_gc] = initial_distribution(type,ks_g(1:2),densities_g,'bimolecular',opt_gp120);
[~, n_ss_gx] = initial_distribution(type,ks_g(3:4),densities_g([3,2,1,4]),'bimolecular',opt_gp120);


figure(11)
coop_env = time_cooperativity(curve_vc,curve_vx,curve_vcx,[n_ss_vc,n_ss_vx]);
hold on
coop_gp120 = time_cooperativity(curve_gc,curve_gx,curve_gcx,[n_ss_gc,n_ss_gx]);



