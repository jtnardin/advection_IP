clear all; clc

IC_str = '_gauss';

load(['advection_art_data' IC_str '_all.mat'])

%params, terms to determine data
alpha   = q0(1);
beta    = q0(2);

%rate of advection
[g,sigma,sigma_inv] = advection_rate('root',alpha,beta);
%initial condition
phi = IC_spec(IC_str(2:end));



%final solution form.
soln = @(t,x) (x>sigma_inv(t,0)).*g(sigma_inv(-t,x))./g(x).*phi(sigma_inv(-t,x));

td = linspace(0,10,6);
xf = linspace(0,1,100);

[Td,X] = meshgrid(td,xd{1});

colors = 'rgbckm';
markers = '*sx^vo';

figure
hold on

for i = 1:length(td)
%     plot(xd{1},data{1,2}(i,:),[colors(i) markers(i)])
    plot(xd{1},data{1,6}(i,:),[colors(i) markers(i)])
    
end
for i = 1:length(td)
    plot(xf,soln(td(i),xf),colors(i))
end

h=legend('$t = 0$','$t = 2$','$t = 4$','$t = 6$','$t = 8$','$t = 10$'...
    ,'$u_0(t,x)$','location','northeast');
set(h,'interpreter','latex')



xlabel('Location ($x$)','interpreter','latex')
ylabel('Density ($u$)','interpreter','latex')

title('Artificial data, $y$','interpreter','latex')
 

exportfig(gcf,['advec_art_data' IC_str '.eps'],'color','rgb','fontsize',1.5)
saveas(gcf,['advec_art_data' IC_str '.fig'])

