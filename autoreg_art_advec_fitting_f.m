%art_advec_fitting_f.m written 2-2-18 by JTN to fit numerical model to
%artifical data from u_t+(g(x)u)_x=0.

function [q_ols,J_ols,q_autoreg,J_autoreg,phi1,phi2,num_it,num_it_auto] =  autoreg_art_advec_fitting_f(xni,m,num_meth,IC_str,model_str,data_str)

    simnum = 5;

    lambda = 1/2;
    

    %initial condition
    phi = IC_spec(IC_str(2:end));


    %create grids for computaiton
    xnsize = [21,41,81,161,321,641,2*640+1];


    %load data
    load(data_str)
    q0 = q0;


    
    xdi = ceil(m/length(eta));
    sigmaj = mod(m,length(eta));
    
    if sigmaj == 0
        sigmaj = length(eta);
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



%     %%%% Now fit to migration data. First initialize q and cost vectors
    q_ols = cell(simnum,1);
    q0_all = cell(simnum,1);
    J_ols = zeros(simnum,1);

    output_ols_all = cell(simnum,1);
    output_auto_all = cell(simnum,1);
    
    phi1 = zeros(6,1);
    phi2 = zeros(6,1);
    
    options = optimset();%'display','iter');

     LB = zeros(2,1);
     UB = inf(2,1);

    
    parfor i = 1:simnum


             q0_all{i} = q0 + .1*randn(size(q0));
            

            % obtain OLS estimate
            [q_ols{i},J_ols(i),~,output_ols_all{i}] = fmincon(@(q) MLE_cost_art_data(cell_data,...
                q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,num_meth,...
                t,tdata,model_str),q0_all{i},[],[],[],[],LB,UB,[],options);
            
    end
%     
    q_ols = q_ols{J_ols == min(J_ols)};
    J_ols = J_ols(J_ols==min(J_ols));
    num_it = output_ols_all{J_ols==min(J_ols)}.funcCount;
            
        %now estimate phihat values
        [~,~,umodel] = MLE_cost_art_data(cell_data,q_ols,dx,xn,x_int,xbd_0,...
            xbd_1,dt,tn,IC,A,Abd,x,xdata,num_meth,t,tdata,model_str);

        res_OLS = umodel - cell_data;

        for j = 2:6
             res_tmp = res_OLS(j,:);

             res_max = max(abs(res_tmp));
             res_max_loc = xdata(abs(res_tmp)==res_max);

                if res_tmp(abs(res_tmp)==res_max) >= 0 %max res value positive

                    res_past_shock = res_tmp(xdata>=res_max_loc);
                    res_before_shock = res_tmp(xdata<res_max_loc);

%                         ind_past_shock = find(xdata>=res_max_loc);
%                         ind_before_shock = find(xdata<res_max_loc);

                elseif res_tmp(abs(res_tmp)==res_max) < 0 %max res value negative


                    res_past_shock = res_tmp(xdata>res_max_loc);
                    res_before_shock = res_tmp(xdata<=res_max_loc);

%                         ind_past_shock = find(xdata>res_max_loc);
%                         ind_before_shock = find(xdata<=res_max_loc);


                end


             phi1(j) = abs(phihat_estimate(res_past_shock));
             phi2(j) = abs(phihat_estimate(res_before_shock));       

        end

        q_autoregc = cell(simnum,1);
        J_autoregc = zeros(simnum,1);
        

        for i = 1:simnum
            if i == 1
                q_guess = q0;
            elseif i == 2
                q_guess = q_ols;
            else
                q_guess = q_ols + .1*randn(size(q_ols));
            end
            % obtain autoreg estimate
            [q_autoregc{i},J_autoregc(i),~,output_auto_all{i}] = fmincon(@(q) MLE_cost_autoreg_art_data(cell_data,...
                q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,num_meth,...
                t,tdata,phi1,phi2),q_guess,[],[],[],[],LB,UB,[],options);
        end

        q_autoreg = q_autoregc{J_autoregc==min(J_autoregc)};
        J_autoreg = J_autoregc(J_autoregc==min(J_autoregc));
        num_it_auto = output_auto_all{J_autoregc==min(J_autoregc)}.funcCount;
    
%     save(['/scratch/summit/jona8898/chem_fitting/chem_fitting_art_data_xdn_'...
%         num2str(xdi) '_xmn_' num2str(xni) '_sigma_' num2str(sigmaj)...
%         '_fl_' num2str(flims) '.mat' ],'q_all','q0_all','J_all')


end

