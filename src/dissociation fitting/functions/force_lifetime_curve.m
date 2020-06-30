function [mean_fit_pts, mean_fit_curve] = force_lifetime_curve(lnp_fit,forces,raw_data)

%loop trhough each force bin to calculate force lifetime curves

n_bins = size(raw_data,2);

%initialize variables
avg_lifetimes_exp = zeros(1,n_bins);
sem_lifetimes_exp = zeros(1,n_bins);
avg_lifetimes_fit = zeros(1,n_bins);
avg_lifetimes_exp_w = zeros(1,n_bins);
F = [];
lifetime = [];


for i=1:n_bins
    %extract data for single force bin
    bin = raw_data{i};
    lifetimes_exp = bin(:,2);
    max_lifetime = max(lifetimes_exp);
    
%     lnp = bin(:,4);
%     survival_p = exp(lnp);
    p_fit = exp(lnp_fit{i});
    
    %calculate mean lifetime for force bin
    avg_lifetimes_exp(i) = mean(lifetimes_exp); %experimental
    sem_lifetimes_exp(i) = std(lifetimes_exp)/sqrt(length(lifetimes_exp));
%     avg_lifetimes_exp_w(i) = mean(lifetimes_exp.*survival_p); %experimental weighted
    
    avg_lifetimes_fit(i) = avg_lifetime(p_fit,max_lifetime); %fit
    
    F = [F; bin(:,1)];
    lifetime = [lifetime; bin(:,2)];
    
end
% figure(1)
e = errorbar(forces,avg_lifetimes_exp,sem_lifetimes_exp,'vertical','s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
e.Color = 'k';
e.LineWidth = 1;
hold on
plot(forces,avg_lifetimes_fit,'linewidth', 1.5)
% e = errorbar(forces,avg_lifetimes_exp_w,sem_lifetimes_exp,'vertical','s','MarkerSize',5,...
%     'MarkerEdgeColor','b','MarkerFaceColor','b');
% e.Color = 'b';
% e.LineWidth = 1;
xlabel('Force (pN)')
ylabel('Lifetime (s)')
% ylim([0 6])
xlim([0 30])

%% interpolate lifetime fit

f_range = linspace(min(forces), max(forces), 300);
if numel(forces) >= 3 % need at least three points to calculate spline
    avg_lifetime_curve = spline(forces,avg_lifetimes_fit,f_range);
    plot(f_range, avg_lifetime_curve,'linewidth', 1.5)
end

% % figure(2)
[F_sort, sort_i] = sort(F);
lifetime = lifetime(sort_i);
% % max_lifetime = movmax(lifetime,100);
% k = convhull(F_sort,lifetime);
% plot(F_sort(k),lifetime(k),'linewidth',1.5)
% scatter(F_sort(k),lifetime(k))
% scatter(F_sort,lifetime)
% 
% lifetime_smooth = smooth(lifetime,0.3,'sgolay');
% lifetime_smooth2 = smooth(lifetime_smooth,20);
lifetime_smooth1 = smooth(lifetime,100);
lifetime_smooth2 = smooth(lifetime_smooth1,25);
% plot(F_sort,lifetime_smooth2,'linewidth',1.5)

% organize output
mean_fit_pts = [forces',avg_lifetimes_fit'];
if numel(forces) >= 3
    mean_fit_curve = [f_range', avg_lifetime_curve'];
else
    mean_fit_curve = [];
end

end

function [avg_lifetime] = avg_lifetime(p_fit,max_lifetime)

%avg lifetime is integral of t*p(t)dt (lifetime weighted by probability)
t = [0:0.1:max_lifetime]';
if length(t) <= length(p_fit)
    p = p_fit(1:length(t));
else
    p = p_fit;
    t = [0:0.1:0.1*(length(p_fit)-1)]';
end
% plot(t,p)

avg_lifetime = sum(t.*p);
% avg_lifetime = sum(p.*0.1);
avg_lifetime = trapz(t,p);
% avg_lifetime = sum(t.*p.*0.1)/sum(p.*0.1);
end

    
    
    
    