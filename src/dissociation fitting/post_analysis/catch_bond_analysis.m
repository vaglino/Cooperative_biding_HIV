function [curve, mean_pts] = catch_bond_analysis(raw_data,type,k,indep_k,f_indep,dens,p0,options)

% determine proper model based on type id

model = pick_model_type_dissociation(type,options);
if type == 5 % analytical solution of trimolecular dissociation for speed
        type = load('tri_diss_solution.mat');
        type = type.tri_diss_solution;
end

n_bins = numel(raw_data); % number of force bins

%initialize variables

mean_t = zeros(n_bins,1); % mean lifetime of each bin
sem_t = zeros(n_bins,1);  % sem of lifetimes
mean_f = zeros(n_bins,1); % mean force of each bin
sem_f = zeros(n_bins,1);  % sem of forces

F = [];
lifetime = [];
max_t = zeros(1,n_bins);

% evaluate and plot experimental force-lifetime curve

for i=1:n_bins %loop trhough each force bin to calculate force lifetime curves

    % extract data for single force bin
    bin = raw_data{i};
    bin_f = bin(:,1);
    bin_t = bin(:,2);
    
    % find mean and sem for exp. t vs f
    mean_f(i) = mean(bin_f);
    mean_t(i) = mean(bin_t);
    sem_f(i) = std(bin_f) / sqrt(length(bin_f));
    sem_t(i) = std(bin_t) / sqrt(length(bin_t)); 
    
    F = [F; bin_f];
    lifetime = [lifetime; bin_t];
    max_t(i) = max(bin_t);
    
end

% plot experimental curve

% figure(1)
e = errorbar(mean_f,mean_t,sem_t,sem_t,sem_f,sem_f,...
    's', 'MarkerSize',5, 'MarkerEdgeColor','k',...
    'MarkerFaceColor','k', 'Color','k', 'linewidth',1);
hold on

% for experimental force range, plot force-lifetime curve obtained from
% model fitting

% first we need to estimate model parameters by interpolating and
% extrapolating from previously fitted parameters for each force level
f_range = linspace(min(mean_f), max(mean_f), 300);

k_interp = zeros(1,length(f_range)); % initialize interpolated vectors
indep_interp = zeros(1,length(f_range));

if isempty(indep_k) % if bimoelcular
%     figure(2)
    k_interp = 10.^interpolate_rates(mean_f,k,f_range,'spline','',false);
    
elseif ~isempty(indep_k) % if trimolecular (there are independent vars)
%     figure(2)
    k_interp = 10.^interpolate_rates(mean_f,k,f_range,'spline','',false);
%     figure(3)
    indep_interp = 10.^interpolate_rates(f_indep,indep_k,f_range,'poly','',false);
end

% then we need to evaluate the model for the entire force range and compute
% the mean lifetime from fit

% figure(4)

% % calculate boundary of lifetimes based on convex hull (does not work as well)
% 
% % clean forces and lifetimes based on lowest and highest force bin
% % F_mask = F >= min(mean_f) & F <= max(mean_f);
% % F = F(F_mask);
% % lifetime = lifetime(F_mask);
% % [boundary,x_bound,y_bound] = upper_conv_hull(F,lifetime);
% % t_limit = interp1(x_bound,y_bound,f_range);

% calculate boundary of maximum lifetime, based on interpolation of
% experimental maxima

t_limit = interp1(mean_f,max_t,f_range);
t_limit = smooth(t_limit, 40); 
% plot(f_range,t_limit)

% t_limit_spline = ppval(boundary,f_range);%spline boundary evaluation

t_fit = zeros(1,length(f_range));
for j=1:length(f_range)
%     if abs(f_range(j) - 7.66) <0.1
%         j = j
%     end
    % evaluate model
    [~,t,p_fit] = model(0:50,type,dens,indep_interp(:,j),k_interp(:,j),p0);
    
    % calculate <t> for force level from fit
    t_limit_at_f = t_limit(j);
    t_fit(j) = avg_lifetime(t,p_fit,t_limit_at_f);
end

% plot fitted force-lifetime curve
% figure(1)

% uncomment for plotting

% plot(f_range,t_fit,'r','linewidth',1.5)
pbaspect([1 1 1])
xlim([0 30])
ylim([0 7])


% organize output
curve = [f_range' t_fit'];
mean_pts = table(mean_f,mean_t,sem_f,sem_t);

end


%% -------------CONVEX HULL FUNCTIONS------------------------------------
% %% function to find lifetime boundary
% function [boundary,x_split,y_split] = upper_conv_hull(F,t)
% 
% % get convex hull
% k = convhull(F,t);
% 
% x_hull = F(k);
% y_hull = t(k);
% scatter(F,t); hold on
% % plot(x_hull,y_hull,'linewidth',1.5)
% % get upper portion of convex hull
% [x_split, y_split] = split_hull(x_hull,y_hull);
% % interpolate upper convex hull
% 
% boundary = spline(x_split,y_split); %may change this to smoothing
% xx = linspace(0,25,101);
% % plot(x_split,y_split,'o',xx,ppval(boundary,xx),'-');
% end
% 
% 
% %% function to get upper part of hull
% function [x_split, y_split] = split_hull(x,y)
% 
% [x_left, left_ind] = min(x);
% y_left = y(left_ind); 
% [x_right, right_ind] = max(x);
% y_right = y(right_ind); 
% % find line equation between two x extreme points
% coeff = polyfit([x_left, x_right], [y_left, y_right], 1);
% % m = coeff(1);
% % b = coeff(2);
% 
% above_line = zeros(1,length(x));
% %determine which hull are above the 
% for i = 1:length(x)
%     y_line_at_x = polyval(coeff,x(i));
%     if y(i) >= y_line_at_x
%         above_line(i) = 1;
%     end
% end
% above_line([left_ind, right_ind]) = 1;
% above_line = logical(above_line);
% x_split = x(above_line);
% y_split = y(above_line);
% plot(x_split,y_split,'--r')
% end