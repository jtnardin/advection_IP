%simulation_convergence.m written 2-1-18 by JTN to compute 
%numerical simulation order convergence based on true-solution comparison

clear all; clc

IC_str = '_front';

load('rel_range_sims.mat')
load(['advection_art_data' IC_str '.mat'])
phi = IC_spec(IC_str(2:end));

%generate true solution
alpha   = q0(1);
beta    = q0(2);

q = [alpha,beta];

%rate of advection
[g,sigma,sigma_inv] = advection_rate('root',alpha,beta);



%final solution form.
soln = @(t,x) (x>sigma_inv(t,0)).*g(sigma_inv(-t,x))./g(x).*phi(sigma_inv(-t,x));

tfin = 10;
xfin = 1;



%create grids for computation
xnsize = [21,41,81,161,321,641,2*640+1,4*640+1];%,8*640+1,16*640+1];%,4*640+1];%,8*640+1];



tdata = 0:2:tfin;
xdata = linspace(0,xfin,xnsize(1));
[Td,Xd] = meshgrid(tdata,xdata);

lambda = 1/2;


sim_order = zeros(4,1);

for j = 1:4

    if j == 1
        num_method = 'upwind';
    elseif j == 2
        num_method = 'laxfried';
    elseif j == 3
        num_method = 'laxwend';
    elseif j == 4
        num_method = 'beamwarm';
    end


    %matrix whose columns are simulations with increasing precision.
    model_sims = zeros(length(xdata),6,length(xnsize));

    %vector whose entries are the true error between numerical sim and true
    %soln (E = ||u(h)-\hat{u}||_1)
    E = zeros(1,length(xnsize));

    tic

    %compute model simulations
    for i = 1:length(xnsize)

            %compute model

            %space
            x = linspace(0,1,xnsize(i));
            dx = x(2) - x(1);
            xn = xnsize(i);

            %time
            dt = lambda*dx;
            t = 0:dt:tfin;
            tn = length(t);


            %points for computaiton
            [x_int,xbd_0,xbd_1] = int_bd_def(xn);

            %create initial condition
            IC = phi(x);

            %load matrices for computation
            [A,Abd] = aMatrixupwind(xn,num_method);

            umodel = advection_computation(q,g,dx,xn,x_int,xbd_0,xbd_1,...
                dt,tn,IC,A,Abd,x,xdata,num_method,t,tdata);


            model_sims(:,:,i) = umodel;


    end

    toc

    %now compute errors and ratios

    udata = soln(Td,Xd);
    udata(isnan(udata))=phi(Xd(isnan(udata)));

    %estimate true order
    for i = 1:length(xnsize)

        E(i) = sum(sum(abs(model_sims(:,:,i)-udata)));

    end
    % 
    % R = log2(E(1:end-1)./E(2:end))

    figure
    loglog(1./xnsize,E,'.-')


    if strcmp(IC_str,'_gauss')
        rel_range = rel_range_gauss{j};
    elseif strcmp(IC_str,'_front')
        rel_range = rel_range_front{j};
    end

    
    fit_line = polyfit(log(1./xnsize(rel_range)),log(E(rel_range)),1);
    sim_order(j) = fit_line(1);
    
end

sim_order