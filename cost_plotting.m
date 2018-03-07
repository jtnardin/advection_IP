clear all; clc


IC_str = '_gauss';


%load best-fit params, data, and initial condition
if strcmp(IC_str,'_gauss')
     load(['advection_rates' IC_str '_IC_3_6.mat'])
    load(['advection_art_data' IC_str '_3_6.mat'])
elseif strcmp(IC_str,'_front')
    load(['advection_rates_autoreg' IC_str '_IC.mat'])
    load(['advection_art_data' IC_str '.mat'])

end
    
    phi = IC_spec(IC_str(2:end));


xnsize = [21,41,81,161,321,641,2*640+1];

num_meth_cell = cell(4,1);
num_meth_cell{1} = 'upwind';
num_meth_cell{2} = 'Lax-Friedrich';
num_meth_cell{3} = 'Lax-Wendroff';
num_meth_cell{4} = 'Beam warming';
num_meth_cell{4} = 'upwind fl';



for i = 1:length(xd)
    xndata(i) = length(xd{i});
end
for i = 1:length(eta)
    eta_vec(i) = eta(i);
end


xnstr = zeros(1,length(xd));
eta_str = zeros(1,length(eta));





for m = 1:numel(data)
    

    xdi = ceil(m/length(eta_vec));
    sigmaj = mod(m,length(eta_vec));
    
    if sigmaj == 0
        sigmaj = length(eta_vec);
    end
    
    xnstr(m) = xndata(xdi);
    eta_str(m) = eta_vec(sigmaj);
    
end


figure('units','normalized','outerposition',[0 0 1 1])

for i = 1:numel(data)
    
    subplot(size(data,1),size(data,2),i)
%     loglog(1./xnsize,J_ols(:,:,i))
    loglog(1./xnsize,squeeze(J_ols(:,i,:)),'.-','linewidth',2,'markersize',20)
    
    hold on
    
    plot([1e-16 1],repmat(eta_str(i)^2,1,2),'k--') 
    
    if i == 1
        
%         h=legend(strcat('$N$ = ',num2str(xnstr'),', $\eta^2 = $',...
%             num2str(eta_str')),'location','northeast');
%         
        if strcmp(IC_str,'_front')
            h=legend('Upwind','Lax-Friedrichs','Lax-wendroff','beam warming',...
                'Upwind FL','$\eta^2$','location','southeast');
        elseif strcmp(IC_str,'_gauss')
            h=legend('Upwind','Lax-Friedrichs','Lax-wendroff','beam warming',...
                '$\eta^2$','location','northwest');
        end
        set(h,'interpreter','latex');
        
    end
    
%     title(['J, ' num_meth_cell{i}])

       
    
    if strcmp(IC_str,'_front')
%         axis([1e-4 1e-1 1e-4 1e0])
          axis([1e-4 1e-1 1e-4 10^.5])
    elseif strcmp(IC_str,'_gauss')
        axis([10^-3.5 10^-1.15 10^-7 1])
    end


    title(strcat('J, $N$ = ',num2str(xnstr(i)),', $\eta^2 = $',...
            num2str(eta_str(i)^2)),'interpreter','latex')

    xlabel('h')
    ylabel('$J(h,\hat{\theta})$','interpreter','latex')

end

% exportfig(gcf,['J_h_plot' IC_str '_3_6.eps'],'fontsize',1.5,'color','rgb')
% saveas(gcf,['J_h_plot' IC_str '_3_6.fig'])


q_norm = zeros(size(q_ols));

for i = 1:size(q_ols,1)
    for j = 1:size(q_ols,2)
        for k = 1:size(q_ols,3)
            
            if ~isempty(q_ols{i,j,k})
                q_norm(i,j,k) = norm(q_ols{i,j,k}-q0);
            end
            
        end
    end
end



figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:numel(data)
    
    subplot(size(data,1),size(data,2),i)
    loglog(1./xnsize,squeeze(q_norm(:,i,:)),'.-','linewidth',2,'markersize',20)
    
    if i == 1
        
%         h=legend(strcat('$N$ = ',num2str(xnstr'),', $\eta^2 = $',...
%             num2str(eta_str')),'location','northeast');

          h=legend('Upwind','Lax-Friedrichs','Lax-wendroff','beam warming',...
            'Upwind FL','location','southeast');
      

        set(h,'interpreter','latex');
        
    end
    
    if strcmp(IC_str,'_front')
        axis([10^-4 10^-1 10^-3.5 10^1])
    elseif strcmp(IC_str,'_gauss')
        axis([10^-3.5 10^-1.15 10^-5 10])
    end
    
%     title(['$\| \hat{\theta} - \theta_0 \|_2$, ' num_meth_cell{i}],'interpreter','latex')
    


    title(strcat('$\| \theta_0 - \hat{\theta} \|_2, N$ = ',num2str(xnstr(i)),', $\eta^2 = $',...
            num2str(eta_str(i)^2)),'interpreter','latex')
    
    xlabel('h')
    ylabel('$\|\cdot\|_2$','interpreter','latex')

    
end


% exportfig(gcf,['q_h_plot' IC_str '_3_6.eps'],'fontsize',1.5,'color','rgb')
% saveas(gcf,['q_h_plot' IC_str '_3_6.fig'])

