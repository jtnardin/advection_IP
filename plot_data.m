clear all; clc

load('advection_art_data_gauss.mat')

%params, terms to determine data
alpha   = .3;
beta    = 0.4;

%rate of advection
[g,sigma,sigma_inv] = advection_rate('root',alpha,beta);
%initial condition
phi = IC_spec('gauss');



%final solution form.
soln = @(t,x) (x>sigma_inv(t,0)).*g(sigma_inv(-t,x))./g(x).*phi(sigma_inv(-t,x));

td = linspace(0,10,6);
xf = linspace(0,1,100);

colors = 'rgbckm';

figure
hold on

for i = 1:length(td)
    plot(xf,soln(td(i),xf),colors(i))
    plot(xd{1},data{1,2}(i,:),[colors(i) '*'])
    
end

xlabel('x')
ylabel('t')

title('artificial data, y')
 

% exportfig(gcf,'advec_art_data.eps','color','rgb','fontsize',1.5)
% saveas(gcf,'advec_art_data.fig')

