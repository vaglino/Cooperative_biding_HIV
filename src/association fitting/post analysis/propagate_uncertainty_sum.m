function [u_sum] = propagate_uncertainty_sum(u1,u2)

% Uncertainity propagation for independent samples 1 and 2 (sum or difference)

u_sum = sqrt(u1.^2 + u2.^2);

end