%art_advec_fitting_f.m written 2-2-18 by JTN to fit numerical model to
%artifical data from u_t+(g(x)u)_x=0.

function [q_ols,J_ols,q_autoreg,J_autoreg,phi1,phi2] =  autoreg_art_advec_fitting_f(xni,m,num_meth,IC_str)

    simnum = 1;

    lambda = 1/2;
    

    %initial condition
    phi = IC_spec(IC_str(2:end));


    %create grids for computaiton
    xnsize = [21,41,81,161,321,641,2*640+1];


    %load data
    load(['advection_art_data' IC_str '.mat'])


    xdi = ceil(m/2);
    sigmaj = mod(m,2);
    
    if sigmaj == 0
        sigmaj = 2;
    end
    
    
    %select exp. data
    cell_data = data{xdi,sigmaj};
    tdata = td;
    xdata = xd{xdi};



    %space
    x = linspace(0,1,xnsize(xni));
    dx = x(2) - x(1);
    xn = xnsize(xni);

    %time
    dt = lambda*dx;
    tfin = 10;
    t = 0:dt:tfin;
    tn = length(t);

    %points for computaiton
    [x_int,xbd_0,xbd_1] = int_bd_def(xn);

    %initial condition
    IC = phi(x);



    %load matrices for computation
    [A,Abd] = aMatrixupwind(xn,num_meth);



    %%%% Now fit to migration data. First initialize q and cost vectors
    q_ols = zeros(simnum,length(q0));
    q0_all = cell(simnum,1);
    J_ols = zeros(simnum,1);

    q_autoreg = zeros(simnum,length(q0));
    J_autoreg = zeros(simnum,1);
    
    phi1 = zeros(1,6);
    phi2 = zeros(1,6);
    
    options = optimset();%'display','iter');

    for i = 1%:simnum


             q0_all{i} = q0;
             LB = zeros(2,1);
             UB = inf(2,1);


            % obtain OLS estimate
            [q_ols(i,:),J_ols(i)] = fmincon(@(q) MLE_cost_art_data(cell_data,...
                q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,num_meth,...
                t,tdata),q0_all{i},[],[],[],[],LB,UB,[],options);
            
            %now estimate phihat values
            [~,~,umodel] = MLE_cost_art_data(cell_data,q_ols(i,:),dx,xn,x_int,xbd_0,...
                xbd_1,dt,tn,IC,A,Abd,x,xdata,num_meth,t,tdata);
            
            res_OLS = umodel - cell_data;
            
            for j = 2:6
                 res_tmp = res_OLS(j,:);

                 res_max = max(abs(res_tmp));
                 res_max_loc = xdata(abs(res_tmp)==res_max);

                 res_past_shock = res_tmp(xdata>=res_max_loc);
                 res_before_shock = res_tmp(xdata<res_max_loc);
                 
                 phi1(j) = abs(phihat_estimate(res_past_shock));
                 phi2(j) = abs(phihat_estimate(res_before_shock));       
                
            end
            
                        
            % obtain OLS estimate
            [q_autoreg(i,:),J_autoreg(i)] = fmincon(@(q) MLE_cost_autoreg_art_data(cell_data,...
                q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,num_meth,...
                t,tdata,phi1,phi2),q0_all{i},[],[],[],[],LB,UB,[],options);
           

    end
    
%     save(['/scratch/summit/jona8898/chem_fitting/chem_fitting_art_data_xdn_'...
%         num2str(xdi) '_xmn_' num2str(xni) '_sigma_' num2str(sigmaj)...
%         '_fl_' num2str(flims) '.mat' ],'q_all','q0_all','J_all')


end

