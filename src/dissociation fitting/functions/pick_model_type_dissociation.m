%% helper function to switch between different models
function [model,n_k,lb,ub] = pick_model_type_dissociation(type,options)

lb = -5;
ub = 5;
switch type
    case 0
        model = @bimolecular_dissociation;%simple bimolecular
        n_k = 1;
%         lb = options.lb;
%         ub = options.lb;
    case 1
        model = @bimolecular_slow_path;%bimolecular with slow intermediate and activation
        n_k = 3;
    case 2
        model = @bimolecular_dissociation_2_species;%concurrent bimolecular dissociation of two species
        n_k = 2;
    case 3
        model = @trimolecular_slow_path;%concurrent bimolecular dissociation with slow paths
        n_k = 1;
        lb = ones(1,n_k).*0;
        ub = ones(1,n_k).*1;
    case 4
        model = @CD4_X4_HIV_dissociation;%trimolecular_dissociation
        n_k = 4;
    case 5
        model = @trimolecular_dissociation_slow;%trimolecular dissociation with slow intermediate
        n_k = 2;
%         lb = ones(1,n_k).*[-6 -3];
%         ub = ones(1,n_k).*[-3 2];
    case 6
        model = @bimolecular_dissociation_reversible;
        n_k = 2;
    case 7
        model = @trimolecular_dissociation_reversible;
        n_k = 4;
    case 8
        model = @trimolecular_dissociation_enhancement;
        n_k = 2;
        lb = ones(1,n_k).*-10;
        ub = ones(1,n_k).*5;
end
%     case 6
%         model = @CD4_X4_HIV_gma;
%         n_k = 2;
%         lb = ones(1,n_k).*-10;
%         ub = ones(1,n_k).*10;
  
end