%MLE_cost_autoreg_art_data.m written 2-12-18 by JTN to compute cost funciton
% while performing autoregression model. In this case, phi1 and phi2 have
% already been computed

function [J,res,umodel] = MLE_cost_autoreg_art_data(cell_data,q_est,dx,xn,x_int,xbd_0,xbd_1,dt,...
    tn,IC,A,Abd,x,xc,num_meth,t,tc,phi1,phi2)
    
    %determine rate of advection from params
    g = advection_rate('root',q_est(1),q_est(2));

    %run simulation using forward euler
    umodel = advection_computation(q_est,g,dx,xn,x_int,xbd_0,xbd_1,...
            dt,tn,IC,A,Abd,x,xc,num_meth,t,tc);
        
    umodel = umodel';
        
    %total data points considered
    N = numel(cell_data);
    res = zeros(length(tc),length(xc));
    
    for i = 1:length(tc)
        
              
        res_tmp = umodel(i,:) - cell_data(i,:);
        
        if i == 1
            
            res(i,:) = res_tmp;
        
        else
            
%             if any(i == [3,5,6])
%                 disp('hey')
%             end
            
            res_max = max(abs(res_tmp));
            res_max_loc = xc(abs(res_tmp)==res_max);

            
            if res_tmp(abs(res_tmp)==res_max) >= 0 %max res value positive

                ind_past_shock = find(xc>=res_max_loc);
                ind_before_shock = find(xc<res_max_loc);

            elseif res_tmp(abs(res_tmp)==res_max) < 0 %max res value negative


                ind_past_shock = find(xc>res_max_loc);
                ind_before_shock = find(xc<=res_max_loc);


            end

            

            [B1,B2] = autoreg_mat(phi1(i),phi2(i),ind_past_shock,ind_before_shock,length(xc));

            res(i,:) = (B1+B2)*(res_tmp');
        
        end
        
    end
    
%     %calculate residuals
%     res = umodel(:) - cell_data(:);

    res = res(:);

    %%% for fmincon
    J = 1/N*sum(res.^2); %slightly biased sample variance
    
    %%% for lsqnonlin
%     J = 1/sqrt(N)*res;

end