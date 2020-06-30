function [u_Q] = propagate_uncertainty_ratio(Q,A,B,uA,uB)

% Uncertainity propagation for independent samples products or ratios

% Finds uncertainty of ratio Q = A/B
%
% Q = A./B;
% uQ/Q = sqrt((uA/A)^2 + (uB/B)^2)
u_Q_over_Q = sqrt( (uA./A).^2 + (uB./B).^2 ); 

% uQ
u_Q = u_Q_over_Q .* Q;

end