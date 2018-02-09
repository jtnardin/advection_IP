%generate_Data.m written 2-2-18 by JTN to generate data from the advection
%equation u_t + (g(x)u)_x = 0.

clear all; clc

td = linspace(0,10,6);

data = cell(2,2);

xndata = [11, 51];
eta = [0 0.1];

xd = cell(3,1);

%params, terms to determine data
alpha   = .3;
beta    = 0.5;

q0 = [alpha,beta];

%rate of advection
[g,sigma,sigma_inv] = advection_rate('root',alpha,beta);
%initial condition
phi = IC_spec('front');

%final solution form.
soln = @(t,x) (g(x)~=0).*(x>=sigma_inv(t,0)).*g(sigma_inv(-t,x))./g(x).*phi(sigma_inv(-t,x));



for i = 1:2
    xd{i} = linspace(0,1,xndata(i));
end

for i = 1:length(xndata)
    for j = 1:length(eta)
        
        
        [Xd,Td] = meshgrid(xd{i},td);
        
        data_0 = soln(Td,Xd);
        
        data_0(isnan(data_0))=phi(Xd(isnan(data_0)));
        
        data_1 = data_0+eta(j)*randn(size(data_0));
        
        
        
        data{i,j} = data_1;
        
    end
end
        

save('advection_art_data_front.mat','data','q0','xd','td','eta')