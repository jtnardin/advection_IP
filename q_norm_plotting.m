IC_str = '_front';


%load best-fit params, data, and initial condition
load(['advection_rates_autoreg' IC_str '_IC.mat'])
load(['advection_art_data' IC_str '.mat'])
phi = IC_spec(IC_str(2:end));



xndata = [10, 50];
eta = [0 0.05];

xnstr = zeros(1,4);
eta_str = zeros(1,4);


for m = 1:4
    
    xdi = ceil(m/2);
    sigmaj = mod(m,2);
    
    if sigmaj == 0
        sigmaj = 2;
    end
    
    xnstr(m) = xndata(xdi);
    eta_str(m) = eta(sigmaj);
    
end


qnorm1 = zeros(7,4);
qnorm2 = zeros(7,4);
for i = 1:7
    for j = 1:4
        qnorm1(i,j) = norm(q_ols{i,j,1}-q0);
        qnorm2(i,j) = norm(q_autoreg{i,j,1}-q0);
    end
end

figure('units','normalized','outerposition',[0 0 1 1])

for i = 1:4
    subplot(2,2,i)
    loglog(1./xnsize,qnorm1(:,i))
    hold on
    loglog(1./xnsize,qnorm2(:,i))
    
    if i == 1
        legend('OLS','autoreg','location','northwest')
    end
    
     title(strcat('$\| \theta_0 - \hat{\theta} \|_2, N$ = ',num2str(xnstr(i)),', $\eta^2 = $',...
            num2str(eta_str(i))),'interpreter','latex')
    
    xlabel('h')
    ylabel('$\|\cdot\|_2$','interpreter','latex')

end


exportfig(gcf,['q_auto_plot' IC_str '.eps'],'fontsize',1.5,'color','rgb')
saveas(gcf,['q_auto_plot' IC_str '.fig'])