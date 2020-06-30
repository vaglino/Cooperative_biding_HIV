function [L] = error_fn(n,y_fit,y_exp,t,options)
%calculate loss

log_space = options.log_space;

if log_space
    % better thermal fluctuation fit with log space error
    L = sqrt((1/n)*sum(((log(y_fit)-log(y_exp)).^2).*t)); 
else
    % better force dependent dissociation fit with linear space error
    L = sqrt((1/n)*sum(((y_fit-y_exp).^2).*(t.^2)));
end

% w = weighting(t);
% L = sqrt((1/n)*sum(((log(y_fit)-log(y_exp)).^2).*w));
end

% function w = weighting(t)
% diff_w = diff([0; t]);
% previous_diff = diff_w(1);
% for i=1:length(diff_w)
%     if diff_w(i)==0
%         diff_w(i)=previous_diff;
%     else
%         previous_diff=diff_w(i);
%     end
% end
% w = diff_w;
% end