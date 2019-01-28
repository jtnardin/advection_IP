IC_str = '_front';


%load best-fit params, data, and initial condition
% load(['advection_rates_autoreg' IC_str '_IC_2_28_BW.mat'])


load('advection_rates_autoreg_front_IC_all.mat')


load(['advection_art_data' IC_str '_all.mat'])
phi = IC_spec(IC_str(2:end));

num_meth = 1;

num_meth_cell = cell(4,1);
num_meth_cell{1} = 'upwind';
num_meth_cell{2} = 'laxfried';
num_meth_cell{3} = 'laxwend';
num_meth_cell{4} = 'beamwarm';

for i = 1:length(xd)
    xndata(i) = length(xd{i});
end
for i = 1:length(eta)
    eta_vec(i) = eta(i);
end

xnsize = [21,41,81,161,321,641,2*640+1];

xnstr = zeros(1,numel(data));
eta_str = zeros(1,numel(data));


for m = 1:numel(data)
    
    xdi = ceil(m/length(eta_vec));
    sigmaj = mod(m,length(eta_vec));
    
    if sigmaj == 0
        sigmaj = length(eta_vec);
    end
   
    
    if sigmaj == 0
        sigmaj = 2;
    end
    
    xnstr(m) = xndata(xdi);
    eta_str(m) = eta_vec(sigmaj);
    
end


qnorm1 = zeros(7,numel(data));
qnorm2 = zeros(7,numel(data));
for i = 1:7
    for j = 1:numel(data)
        qnorm1(i,j) = norm(q_ols{i,j,num_meth}-q0);
        qnorm2(i,j) = norm(q_autoreg{i,j,num_meth}-q0);
    end
end

figure('units','normalized','outerposition',[0 0 1 1])

for i = 1:numel(data)
    subplot(size(data,1),size(data,2),i)
    loglog(1./xnsize,qnorm1(:,i),'linewidth',2)
    hold on
    loglog(1./xnsize,qnorm2(:,i),'linewidth',2)
    
    if i == 1
        legend('OLS','autoreg','location','northwest')
    end
    
     title(strcat('$N = ',num2str(xnstr(i)),', \eta^2 = $',...
            num2str(eta_str(i))),'interpreter','latex')
    
    xlabel('h')
    ylabel('$\| \theta_0 - \hat{\theta} \|_2,$','interpreter','latex')
    

    if num_meth == 1
        axis([10^-3.5 10^-1.15 10^-3 1])
    elseif num_meth == 2
        axis([10^-3.5 10^-1.15 10^-1 10])
    end

end


exportfig(gcf,['q_auto_compare' IC_str '_' num_meth_cell{num_meth} '.eps'],'fontsize',1.5,'color','rgb')
saveas(gcf,['q_auto_compare' IC_str '_' num_meth_cell{num_meth} '.fig'])