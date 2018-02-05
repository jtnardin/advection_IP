%simulation_convergence.m written 2-1-18 by JTN to compute 
%numerical simulation order convergence based on true-solution comparison

clear all; clc
%generate true solution
alpha   = .3;
beta    = 0.4;

num_method = 'beamwarm';

q = [alpha,beta];

%rate of advection
[g,sigma,sigma_inv] = advection_rate('root',alpha,beta);


%initial condition
phi = IC_spec('step');



%final solution form.
soln = @(t,x) (x>sigma_inv(t,0)).*g(sigma_inv(-t,x))./g(x).*phi(sigma_inv(-t,x));

tfin = 10;
xfin = 1;



%create grids for computation
xnsize = [21,41,81,161,321,641,2*640+1];%,4*640+1];%,8*640+1];

tdata = 0:2:tfin;
xdata = linspace(0,xfin,xnsize(1));
[Td,Xd] = meshgrid(tdata,xdata);

lambda = 1/2;


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
udata(isnan(udata))=0;

%estimate true order
for i = 1:length(xnsize)
    
    E(i) = sum(sum(abs(model_sims(:,:,i)-udata)));
    
end

R = log2(E(1:end-1)./E(2:end))