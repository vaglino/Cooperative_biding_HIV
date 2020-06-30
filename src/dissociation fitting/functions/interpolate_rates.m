function [k_new] = interpolate_rates(f,k,f_new,algorithm,species,plotting)

% INTERPOLATE_RATES
% function to interpolate force dependent kinetic rates at new force levels
% so that rate constants can be obtained for trimolecular models
% INPUTS
% f - original force levels
% k - rate constants at f
% f_new - new force levels we want to interpolate at
% algorithm - either spline or poly
% species - string for figure title
% plotting - either true or false
% OUTPUTS
% k_new - rates at new_f

% Stefano Travaglino, Zhu Lab, 2020
% -------------------------------------------------------------------------

%            f1 - - fn
% k = [koff1 - - - - -;
%      kint  - - - - -;
%      kslow - - - - -]

n_k = size(k,1);
cm = colormap(hsv(n_k));
k_new = zeros(n_k,length(f_new));
interp_save = zeros(100,n_k);
for i=1:n_k
%     
    %uncomment for plotting
    if isempty(f_new)
        k_max = max(f);
        k_min = min(f);
    else
        k_max = max(max(f_new),max(f));
        k_min = min(min(f_new),min(f));
    end
%     x1 = linspace(0,k_max);
    if plotting
        scatter(f,k(i,:),30,cm(i,:),'filled')%scatter ks from bimolecular
        hold on
    end
    
    xx = linspace(k_min, k_max,100);
    switch algorithm
        case 'spline'
            % spline fitting
            p = csape(f,k(i,:),'second',[0 0]); %interpolate with spline and bondary conditions
            y1 = ppval(p,xx);
            if plotting
                plot(xx,y1,'color',cm(i,:),'linewidth',1.5);
            end
            k_new(i,:) = fnval(p,f_new); %evaluate interpolated line at new force values
        case 'poly'
            % polynomial fit
            p3 = polyfit(f,k(i,:),2);
            y1 = polyval(p3,xx);
            if plotting
                plot(xx,y1,'color',cm(i,:),'linewidth',1.5)   
            end
            k_new(i,:) = polyval(p3,f_new);
    end
    if plotting
        scatter(f_new, k_new(i,:),30,'k');
        ylabel('log(k_i)')
        xlabel('F (pN)')
        title(species)
        xlim([0, 30])
        ylim([-5, 3])
    end
    
    interp_save(:,i) = y1';
    
end
xx = xx';
end
% hold off