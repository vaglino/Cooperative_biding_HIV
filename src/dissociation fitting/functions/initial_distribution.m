function [n_norm, n_ss] = initial_distribution(type,ks,densities,interaction,options)


[model, n_k] = pick_model_type(type,options);

switch interaction
    case 'bimolecular'
        indep_vars = [];
        k = ks;
    case 'trimolecular'
        indep_vars = ks(1:end-n_k);
        k = ks(end-n_k+1 : end);
end

[~,~,~,n] = model([0; 10],type,densities,indep_vars,k);

n_ss = n(:,end);
n_norm = n_ss ./ sum(n_ss);

end