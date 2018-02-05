%MLE_cost_D0 written 6-13-17 by JTN to do 29-point WLS scheme

function [J,res,umodel] = MLE_cost_art_data(cell_data,q_est,dx,xn,x_int,xbd_0,xbd_1,dt,...
    tn,IC,A,Abd,x,xc,num_meth,t,tc)
    
    %determine rate of advection from params
    g = advection_rate('root',q_est(1),q_est(2));

    %run simulation using forward euler
    umodel = advection_computation(q_est,g,dx,xn,x_int,xbd_0,xbd_1,...
            dt,tn,IC,A,Abd,x,xc,num_meth,t,tc);
        
    umodel = umodel';
        
    %total data points considered
    N = numel(cell_data);
    
    
    %calculate residuals
    res = umodel(:) - cell_data(:);
    
    %%% for fmincon
    J = 1/N*sum(res.^2); %slightly biased sample variance
    
    %%% for lsqnonlin
%     J = 1/sqrt(N)*res;

end