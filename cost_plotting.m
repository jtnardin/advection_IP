clear all; clc


IC_str = '_step';


%load data and best-fit rates
if strcmp(IC_str,'_step')
    
    load('advection_rates_step_IC_2.mat')
    load('advection_art_data.mat')
    
    phi = IC_spec('step');
    
elseif strcmp(IC_str,'_gauss')

    load('advection_rates_gauss_IC.mat')
    load('advection_art_data_gauss.mat')
    
    phi = IC_spec('gauss');
    
end




xnsize = [21,41,81,161,321,641,2*640+1];

num_meth_cell = cell(4,1);
num_meth_cell{1} = 'upwind';
num_meth_cell{2} = 'Lax-Friedrich';
num_meth_cell{3} = 'Lax-Wendroff';
num_meth_cell{4} = 'Beam warming';


xndata = [10, 50];
eta = [0 0.05];

xnstr = zeros(1,4);
eta_str = zeros(1,4);

order_table = zeros(4,4);
poly_table = cell(4,4);

for i = 1:4
    for j = 1:4
        
        poly_table{i,j} = polyfit(log(1./xnsize)',log(squeeze(J(:,i,j))),1);
        order_table(i,j) = poly_table{i,j}(1);

        
        
    end
end



for m = 1:4
    
    xdi = ceil(m/2);
    sigmaj = mod(m,2);
    
    if sigmaj == 0
        sigmaj = 2;
    end
    
    xnstr(m) = xndata(xdi);
    eta_str(m) = eta(sigmaj);
    
end


figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:4
    
    subplot(2,2,i)
    loglog(1./xnsize,J(:,:,i))
    
    if i == 1
        
        h=legend(strcat('$N$ = ',num2str(xnstr'),', $\eta^2 = $',...
            num2str(eta_str')),'location','northeast');
        
        set(h,'interpreter','latex');
        
    end
    
    title(['J, ' num_meth_cell{i}])

    xlabel('h')
    ylabel('$J(h,\hat{\theta})$','interpreter','latex')

end

exportfig(gcf,['J_h_plot' IC_str '.eps'],'fontsize',1.5,'color','rgb')
saveas(gcf,['J_h_plot' IC_str '.fig'])


q_norm = zeros(7,4,4);

for i = 1:7
    for j = 1:4
        for k = 1:4
            
            if J(i,j,k)>0
                q_norm(i,j,k) = norm(q{i,j,k}-q0);
            end
            
        end
    end
end



figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:4
    
    subplot(2,2,i)
    loglog(1./xnsize,q_norm(:,:,i))
    
    if i == 1
        
        h=legend(strcat('$N$ = ',num2str(xnstr'),', $\eta^2 = $',...
            num2str(eta_str')),'location','northeast');
        
        set(h,'interpreter','latex');
        
    end
    
    axis([10^-4 10^-1 10^-5 10^1])
    
    title(['$\| \hat{\theta} - \theta_0 \|_2$, ' num_meth_cell{i}],'interpreter','latex')
    
    
    xlabel('h')
    ylabel('$J(h,\hat{\theta})$','interpreter','latex')

    
end


exportfig(gcf,['q_h_plot' IC_str '.eps'],'fontsize',1.5,'color','rgb')
saveas(gcf,['q_h_plot' IC_str '.fig'])

