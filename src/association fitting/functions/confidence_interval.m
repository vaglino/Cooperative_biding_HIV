function CI = confidence_interval(dist, C)
% dist is distribution, as row vector, or matrix where rows are
% distributions
% C is confidence level specified

n = size(dist,2);
dof = n-1;
alpha = (1 - C/2);

% evaluate t dist at alpha level
t_star = cdf('T', alpha, dof); 
% evaluate standard deviation
s = std(dist,[],2); 

% calculate confidence interval
CI = t_star .* (s./sqrt(n));

end


