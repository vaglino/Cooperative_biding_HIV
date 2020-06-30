%       Env     rgp120
RSS = [1.3801	0.084123809   % path I
        1.1395	0.071863134   % 2 path
        2.6986	0.21466       % path II
        1.1395	0.071863134]; % 2 path

% Env
% pathway I
[AIC_e1, w_e1, rel_l_e1] = AIC_analysis([RSS(1,1);RSS(2,1)],[2;4]);

% pathway II
[AIC_e2, w_e2, rel_l_e2] = AIC_analysis([RSS(3,1);RSS(4,1)],[2;4]);


% rgp120
% pathway I
[AIC_g1, w_g1, rel_l_g1] = AIC_analysis([RSS(1,2);RSS(2,2)],[2;4]);

% pathway II
[AIC_g2, w_g2, rel_l_g2] = AIC_analysis([RSS(3,2);RSS(4,2)],[2;4]);


function [AICc, w, rel_l_mat] = AIC_analysis(RSS,k)
    
% raw = readmatrix('models_comparison_matrix.xlsx');
% RSS = raw(:,3);
% k = raw(:,2);
n = 43;

AICc = AIC_corrected(RSS, k, n);

% [AIC_min, AIC_min_i] = min(AICc);

delta_i = AICc - min(AICc);
likelihoods = exp(-delta_i./2);
sum_l = sum(likelihoods);
w = likelihoods ./ sum_l; %this is w_i

rel_l_mat = (w' * 1./w)';
rel_l_mat = triu(rel_l_mat);
    
end

function AIC = AIC_corrected(RSS, k, n)
    AIC = n .* log(RSS./n) + 2.*k; % AIC, where MLE of variance is RSS/n
    correction = 2.*k.*(k+1) ./ (n-k-1); % correction for small sample size
    AIC = AIC + correction;
    % AICc = AIC;
end