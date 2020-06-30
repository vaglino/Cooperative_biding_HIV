function cooperativity = time_cooperativity(R,X,RX,n)
% R = receptor only curve
% X = coreceptor only curve
% RX = receptor + coreceptor curve
% all inputs are n by 2 matrix, where first columns is all force levels and
% second column is fitted mean lifetime
% n = [nr, nx, ntri] average number of bonds at steady state

Rf = R(:,1);
Rt = R(:,2);
Xf = X(:,1);
Xt = X(:,2);
RXf = RX(:,1);
RXt = RX(:,2);

% low and up boundary are determine by shortest fit

lb = max([min(Rf), min(Xf), min(RXf)]); % maximum of minima
ub = min([max(Rf), max(Xf), max(RXf)]); % minumum of maxima

% resample each curve (each curve has different force levels), so that we
% can compare them

f_resample = linspace(lb, ub, 100)';

Rt_res = spline(Rf,Rt,f_resample);
Xt_res = spline(Xf,Xt,f_resample);
RXt_res = spline(RXf,RXt,f_resample);

% calculate cooperativity index

coop = coop_index(Rt_res,Xt_res,RXt_res, n, f_resample);

% plot cooperativity index

plot(f_resample, coop, 'linewidth', 1.5)

% compile output
cooperativity = table(f_resample,coop);
end

function coop = coop_index(Rt,Xt,RXt, n, F) 

nR = n(1);
nX = n(2);

predicted = (Rt.*nR + Xt.*nX) ./ (nR + nX);

coop = (RXt ./ predicted) - 1;

end