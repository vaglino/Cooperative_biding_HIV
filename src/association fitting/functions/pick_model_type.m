%% helper function to switch between different models
function [model,n_k,lb,ub] = pick_model_type(type,options)

% PICK_MODEL_TYPE
% given model ID number, returns model function handle as well as number of
% free parameters and lower and upper boundary for optimization

% Stefano Travaglino, ZhuLab, 2020

lb = -10; % defaults
ub = 4;


switch type
    case 0
        model = @bimolecular_model;
        n_k = 2;
    case 1
        model = @CD4_X4_HIV_association_2step;
        n_k = 2;
        lb = options.lb;
        ub = options.ub;
    case 2
        model = @CD4_X4_HIV_association;
        n_k = 4;
        lb = options.lb;
        ub = options.ub;
%         lb = ones(1,n_k).*-10;
%         ub = ones(1,n_k).*3;
    case 3
        model = @CD4_X4_HIV_association_no_off3;
        n_k = 3;
    case 4
        model = @CD4_X4_HIV_association_no_3;
        n_k = 2;
%         lb = ones(1,n_k).*-10;
%         ub = ones(1,n_k).*2;
    case 5
        model = @CD4_X4_HIV_association_no_on3;
        n_k = 3;
    case 6
        model = @CD4_X4_HIV_gma;
        n_k = 2;
%         lb = ones(1,n_k).*-10;
%         ub = ones(1,n_k).*10;
    case 7
        model = @bimolecular_model_coop_zhu;
        n_k = 2;   
        lb = ones(1,n_k).*0;
        ub = ones(1,n_k).*5;
    case 8
        model = @CD4_X4_HIV_gma_alpha_only;
        n_k = 1;
%         lb = ones(1,n_k).*0;
%         ub = ones(1,n_k).*10;
        lb = ones(1,n_k).*-10;
        ub = ones(1,n_k).*10;
    case 9
        model = @bimolecular_model_coop_zhu_a_only;
        n_k = 1;   
        lb = ones(1,n_k).*0;
        ub = ones(1,n_k).*5;
        
end
end