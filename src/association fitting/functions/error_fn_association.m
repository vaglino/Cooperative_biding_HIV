function [L,sse] = error_fn_association(n,y_fit,y_exp,t,options)
%calculate loss
% w = weigthting(t);
time_weighted = options.time_weighted;

if time_weighted
    L = sqrt((1/(n)*sum(((y_fit-y_exp).^2).*t))); % time weighted (Env)
    sse = sum(((y_fit-y_exp).^2).*t); % time weighted
else
    L = sqrt((1/(n)*sum(((y_fit-y_exp).^2)))); % (rgp120)
    sse = sum(((y_fit-y_exp).^2));
end
    
% L = sqrt((1/(n)*sum(((y_fit-y_exp).^2).*t))); % time weighted (Env)
%     sse = sum(((y_fit-y_exp).^2).*t); % time weighted    
end

% function w = weigthting(t)
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