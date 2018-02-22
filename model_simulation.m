%Used to check local minima from optimizer
%tell which scenario you want to sim, specify parameter, get resulting cost
%function

clear all;clc

%%%%specify which initial condition, method,
%%%%and data set we're compating
IC_str = '_front';
num_meth = 4;
m = 7;
xni = 7;

%load best-fit params, data, and initial condition

if strcmp(IC_str,'_front')
    load(['advection_rates_autoreg' IC_str '_IC.mat'])
elseif strcmp(IC_str,'_gauss')
    load(['advection_rates' IC_str '_IC.mat'])
end

load(['advection_art_data' IC_str '.mat'])
phi = IC_spec(IC_str(2:end));


for i = 1:2

    if i == 1
        q = q_ols{xni,m,num_meth};
    elseif i == 2
        q = q_ols{xni-1,m,num_meth};
    end

    %grid sizes
    xnsize = [21,41,81,161,321,641,2*640+1];
    lambda = 1/2;

    xndata = [length(xd{1}), length(xd{2})];


    xdi = ceil(m/length(eta));
    sigmaj = mod(m,length(eta));

    if sigmaj == 0
        sigmaj = length(eta);
    end

    xnstr(m) = xndata(xdi);
    eta_str(m) = eta(sigmaj);

    %select exp. data
    cell_data = data{xdi,sigmaj};
    tdata = td;
    xdata = xd{xdi};

    num_meth_cell = cell(4,1);
    num_meth_cell{1} = 'upwind';
    num_meth_cell{2} = 'laxfried';
    num_meth_cell{3} = 'laxwend';
    num_meth_cell{4} = 'beamwarm';






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
    [A,Abd] = aMatrixupwind(xn,num_meth_cell{num_meth});

    %get model sim
    [J,res,model_sims] = MLE_cost_art_data(cell_data,...
        q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,...
        num_meth_cell{num_meth},t,tdata);

    J

end

