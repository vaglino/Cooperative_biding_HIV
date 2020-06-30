function [mat] = clean_data(raw_data)

%data is organized by columns [Force, lifetime, ln(n), ln(p)] ... xN
mat = {};
n = size(raw_data,2)./4;  %4 columns per force bin
for i=0:n-1
    ind_bin = i*4+1:i*4+4;
    data_bin = raw_data(:,ind_bin);
    mask_zeros = data_bin(:,1) == 0; %data is padded with zeros at the end
    data_bin(mask_zeros,:) = []; %delete zeros at the end
    if size(data_bin,1)>1
        mat{i+1} = data_bin; %organize data for each force bin in a cell
    end
end
