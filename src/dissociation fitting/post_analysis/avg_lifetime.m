function [avg_lifetime] = avg_lifetime(t,p_fit,max_t)

%avg lifetime is integral of t*p(t)dt (lifetime weighted by probability)

% find lifetimes below maximum t
t_mask = t <= max_t;
t = t(t_mask);
p_fit = p_fit(t_mask);

% p_i = p_fit >= 1e-2;
% t = t(p_i);
% p_fit = p_fit(p_i);
avg_lifetime = trapz(t,p_fit); % numerical integral

end