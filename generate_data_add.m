%generate_Data.m written 2-2-18 by JTN to generate data from the advection
%equation u_t + (g(x)u)_x = 0.

clear all; clc

td = linspace(0,10,6);

data = cell(2,2);

IC_str = 'front';

if strcmp(IC_str,'gauss')
    %%%%% gauss
    load('advection_art_data_gauss.mat')
    xndata = [10, 50];
    eta = [0 0.05 .01 .02];
    %params, terms to determine data
    alpha   = .3;
    beta    = 0.4;
    
    data_old = data;

elseif strcmp(IC_str,'front')
    %%%%% front
    load('advection_art_data_front.mat')
    xndata = [11, 51];
    eta = [0 0.1 .2 .3];
    %params, terms to determine data
    alpha   = .3;
    beta    = 0.5;
    
    data_old = data;
end

xd = cell(length(xndata),1);

q0 = [alpha,beta];

%rate of advection
[g,sigma,sigma_inv] = advection_rate('root',alpha,beta);
%initial condition
phi = IC_spec(IC_str);

%final solution form.
soln = @(t,x) (g(x)~=0).*(x>=sigma_inv(t,0)).*g(sigma_inv(-t,x))./g(x).*phi(sigma_inv(-t,x));



for i = 1:length(xndata)
    xd{i} = linspace(0,1,xndata(i));
end

for i = 1:length(xndata)
    for j = 1:length(eta)
        
       %re-use old data
        if j <= 2
        
            data{i,j} = data_old{i,j};
            
        else
            
            [Xd,Td] = meshgrid(xd{i},td);

            data_0 = soln(Td,Xd);

            data_0(isnan(data_0))=phi(Xd(isnan(data_0)));

            data_1 = data_0+eta(j)*randn(size(data_0));



            data{i,j} = data_1;
        end
        
    end
end
        

save(['advection_art_data_' IC_str '.mat'],'data','q0','xd','td','eta')