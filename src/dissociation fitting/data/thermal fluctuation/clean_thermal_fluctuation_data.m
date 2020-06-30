%% cleane thermal fluctuation data and calculate lnp for each condition

sheet_name_survival = 'thermal_survival_HIV.xlsx'
data_vcx_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct VLP vs CD4.X4.J');
data_vc_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct VLP vs CD4.J');
data_vx_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct VLP vs X4.J');
data_gcx_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct gp120 vs CD4.X4.J');
data_gc_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct gp120 vs CD4.J');
data_gx_thermal = readmatrix(sheet_name_survival,'Sheet','Therm fluct gp120 vs X4.J');

lnp_vcx = calculate_lnp(data_vcx_thermal(:,1));
lnp_vc = calculate_lnp(data_vc_thermal(:,1));
lnp_vx = calculate_lnp(data_vx_thermal(:,1));
lnp_gcx = calculate_lnp(data_gcx_thermal(:,1));
lnp_gc = calculate_lnp(data_gc_thermal(:,1));
lnp_gx = calculate_lnp(data_gx_thermal(:,1));

function lnp = calculate_lnp(t_bin)
bin_size = length(t_bin);
survival_bin = (bin_size.*ones(bin_size,1)-[0:bin_size-1]')./bin_size;
lnp = log(survival_bin);
scatter(t_bin,lnp,20,'filled','MarkerFaceAlpha',0.5)
end
