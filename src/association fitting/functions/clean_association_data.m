function [t, Pa] = clean_association_data(num)
%num is raw data matrix provided by Ke, columns 5 to end are raw data

raw = num(:,5:end); 
n_col = size(raw,2); % n of raw data columns, each is one experiment

raw_vert = raw(1:end)'; % rearrange raw data vertically, stacking columns
t_vert = repmat(num(:,1), n_col, 1); %create vertical vector of repeated timepoints

% remove missing data points (nans)
id_missing = isnan(raw_vert);
Pa = raw_vert(~id_missing); 
t = t_vert(~id_missing);
end


