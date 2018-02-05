%art_advec_fitting_f.m written 2-2-18 by JTN to fit numerical model to
%artifical data from u_t+(g(x)u)_x=0.

function [q_f,J_f] =  art_advec_fitting_f(xni,m,num_meth)

    simnum = 5;

    lambda = 1/2;
    

    %initial condition
    phi = IC_spec('step');


    %create grids for computaiton
    xnsize = [21,41,81,161,321,641,2*640+1];


    %load data
    load('advection_art_data.mat')


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
    q_all = cell(simnum,1);
    q0_all = cell(simnum,1);
    J_all = zeros(simnum,1);

    options = optimset();%'display','iter');

    for i = 1%:simnum

             %q = [chi,alpha,beta,D_v]^T;
             q0_all{i} = [1 1] + .75*randn(1,2);
             LB = zeros(2,1);
             UB = inf(2,1);

             while any(q0_all{i})<0
                 q0_all{i} = [1 1] + .75*randn(1,2);
             end

%             tic
            [q_all{i},J_all(i)] = fmincon(@(q) MLE_cost_art_data(cell_data,...
                q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,num_meth,...
                t,tdata),q0_all{i},[],[],[],[],LB,UB,[],options);

%             toc
    end

    q_f = q_all{1};
    J_f = J_all(1);
    
%     save(['/scratch/summit/jona8898/chem_fitting/chem_fitting_art_data_xdn_'...
%         num2str(xdi) '_xmn_' num2str(xni) '_sigma_' num2str(sigmaj)...
%         '_fl_' num2str(flims) '.mat' ],'q_all','q0_all','J_all')


end

