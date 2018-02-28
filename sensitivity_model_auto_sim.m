%MLE_cost_D0 written 6-13-17 by JTN to do 29-point WLS scheme

function [umodel,s1_auto_model,s2_auto_model,J,res] = sensitivity_model_auto_sim(cell_data,q_est,dx,xn,x_int,xbd_0,xbd_1,dt,...
    tn,IC,A,Abd,x,xc,num_meth,t,tc,phi1,phi2)
    
    %determine rate of advection from params
    [g,~,~,galpha,gbeta] = advection_rate('root',q_est(1),q_est(2));

    %run simulation using forward euler
    [umodel,s1model,s2model] = advection_sens_computation(g,galpha,gbeta...
        ,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xc,num_meth,t,tc);
        
       
    %total data points considered
    N = numel(cell_data);
    res = zeros(length(tc),length(xc));
    s1_auto_model = zeros(length(tc),length(xc));
    s2_auto_model = zeros(length(tc),length(xc));
    
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
                s1_auto_model(i,:) = (B1+B2)*s1model(i,:)';
                s2_auto_model(i,:) = (B1+B2)*s2model(i,:)';

            end

    end

    %     %calculate residuals
    %     res = umodel(:) - cell_data(:);

        res = res(:);

        %%% for fmincon
        J = 1/(N-2)*sum(res.^2); %slightly biased sample variance

            
end