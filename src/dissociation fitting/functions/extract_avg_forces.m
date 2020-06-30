function [avg_forces] = extract_avg_forces(mat_bins)

n_bins = size(mat_bins,2);
avg_forces=[];
for i=1:n_bins
    %extract data for single force bin
    bin = mat_bins{i};
    forces = bin(:,1);
  
    avg_force = mean(forces);
    avg_forces = [avg_forces avg_force]; 
end
end