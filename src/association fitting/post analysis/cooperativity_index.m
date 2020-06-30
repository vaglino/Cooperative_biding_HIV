function [coop_indx, sd_coop_indx] = cooperativity_index(varargin)

% COOPERATIVITY_INDEX
% calculates cooperativity index based on the adhesion data for bimolecular
% bonds alone and that of a mix of both species

% parse inputs
t = varargin{1};   % time pts
n1 = varargin{2};  % bimolecular 1
n2 = varargin{3};  % bimolecular 2
n3 = varargin{4};  % trimolecular
sd1 = varargin{5}; % respective sds ...
sd2 = varargin{6};
sd3 = varargin{7};


sum_bi = n1+n2; % sum of bimolecular bonds
sum_tri = n3; % sum of bonds when trimolecular allowed to form
% sum_tri = sum_tri-sum_tri(2)+sum_bi(2);


% propagate uncertainity for difference (sum_tri - sum_bi)
sd_sum_bi = propagate_uncertainty_sum(sd1,sd2);

% delta_n
delta_n = (sum_tri-sum_bi); 

% propagate uncertainity for difference (sum_tri - sum_bi)
sd_delta_n = propagate_uncertainty_sum(sd3, sd_sum_bi);

% cooperativity index is delta_n 
coop_indx = delta_n ./ sum_bi;
% coop_indx = coop_indx - coop_indx(2);
% propagate uncertainity for coop ratio delta_n ./ sum_bi
sd_coop_indx = propagate_uncertainty_ratio(coop_indx,delta_n,sum_bi, ...
                                        sd_delta_n,sd_sum_bi);

                                    
coop_indx = coop_indx';
% coop_indx = coop_indx - coop_indx(2);
coop_indx = coop_indx*100; % convert to percent
sd_coop_indx = abs(sd_coop_indx)*100;
plot(t,coop_indx,'b','linewidth',1.5)
errorbar(t,coop_indx,sd_coop_indx)
end



