clear all; clc


IC_str = '_front';


%load best-fit params, data, and initial condition
if strcmp(IC_str,'_gauss')
    load(['advection_rates' IC_str '_IC_all_3_26.mat'])
    load(['advection_art_data' IC_str '_all_3_26.mat'])
elseif strcmp(IC_str,'_front')
    load(['advection_art_data' IC_str '_all_rev.mat'])
    load('advection_rates_autoreg_front_IC_rev.mat')
%     
%     load(['advection_art_data' IC_str '_all_missing.mat'])
%     load('advection_rates_front_IC_missing.mat')
    

end
    
    phi = IC_spec(IC_str(2:end));


xnsize = [21,41,81,161,321,641,2*640+1];

if strcmp(IC_str,'_front')

    num_meth_cell = cell(5,1);
    num_meth_cell{1} = 'upwind';
    %num_meth_cell{2} = 'Lax-Friedrich';
    num_meth_cell{2} = 'Lax-Wendroff';
    num_meth_cell{3} = 'Beam-Warming';
    num_meth_cell{4} = 'upwind fl';

elseif strcmp(IC_str,'_gauss')

    num_meth_cell = cell(5,1);
    num_meth_cell{1} = 'upwind';
    num_meth_cell{2} = 'Lax-Wendroff';
    num_meth_cell{3} = 'Beam-Warming';
    num_meth_cell{4} = 'upwind fl';
end


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


markers = '.^sxv';
markersize = [15 8 8 8 8];


for i = 1:numel(data)
    
    subplot(size(data,1),size(data,2),i)
%     loglog(1./xnsize,J_ols(:,:,i))
    
    if strcmp(IC_str,'_gauss')
%        loglog(1./xnsize,squeeze(J_ols(:,i,:)),'.-','linewidth',2,'markersize',17)
       
       markerct = 1;

        for k = [1 2 3 4] 
            loglog(1./xnsize,squeeze(J_ols(:,i,k)),[markers(markerct) '-'],'linewidth',1,'markersize',markersize(markerct))

            if k == 1
                hold on
            end

            markerct = markerct+1;
        end
    
       
    elseif strcmp(IC_str,'_front')
%        loglog(1./xnsize,squeeze(J_ols(:,i,[1 3:5])),'.-','linewidth',2,'markersize',17)

        markerct = 1;
    
        for k = [1 2 3 4] 
            loglog(1./xnsize,squeeze(J_ols(:,i,k)),[markers(markerct) '-'],'linewidth',1,'markersize',markersize(markerct))

            if k == 1
                hold on
            end

            markerct = markerct+1;
        end
    
    end
    
       
      
    plot([1e-16 1],repmat(eta_str(i)^2,1,2),'k--') 
    
    if i == numel(data)
        
        h=legend(strcat('$N$ = ',num2str(xnstr'),', $\eta^2 = $',...
            num2str(eta_str')),'location','northeast');
        
        if strcmp(IC_str,'_front')
            h=legend('Upwind','Lax-wendroff','Beam-Warming',...
                'Upwind FL','$\eta^2$','location','southeast');
        elseif strcmp(IC_str,'_gauss')
            h=legend('Upwind','Lax-wendroff','beam warming','upwind FL',...
                '$\eta^2$','location','northwest');
        end
        set(h,'interpreter','latex','units','normalized','position',[.9 .7 .05 .15]);
        
    end
    

       
    
    if strcmp(IC_str,'_front')
%         axis([1e-4 1e-1 1e-4 1e0])
          axis([10^-3.5 1e-1 10^-(4) 10^.5])
    elseif strcmp(IC_str,'_gauss')
        axis([10^-3.5 10^-1 10^-8 10^-1])
    end


    if i <= length(eta)
          title(['$\mathbf{\eta^2 = ',sprintf('%1.e',eta_str(i)^2) '}$'],...
              'interpreter','latex','fontsize',15)
    end
            
    if mod(i,length(eta)) == 1
        ylabel(['$\mathbf{N = ' num2str(xnstr(i)) '}$'],'interpreter',...
            'latex','fontsize',15)
    end
    
    xlabel('$h$','interpreter','latex')
    xticks([10^-3 10^-2])

end


exportfig(gcf,['J_h_plot' IC_str '.eps'],'fontsize',1.5,'color','rgb')
saveas(gcf,['J_h_plot' IC_str '.fig'])


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

if strcmp(IC_str,'_front')
    qnorma = zeros(size(q_autoreg));
    
    
    for i = 1:size(q_autoreg,1)
        for j = 1:size(q_autoreg,2)
            for k = 1

                if ~isempty(q_autoreg{i,j,k})
                    q_norma(i,j,k) = norm(q_autoreg{i,j,k}-q0);
                end
            end
        end
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:numel(data)
    
    subplot(size(data,1),size(data,2),i)
    if strcmp(IC_str,'_gauss')
        markerct = 1;
        for k = [1 2 3 4] 
            loglog(1./xnsize,squeeze(q_norm(:,i,k)),[markers(markerct) '-'],'linewidth',1,'markersize',markersize(markerct))

            if k == 1
                hold on
            end

            markerct = markerct+1;
        end
    elseif strcmp(IC_str,'_front')
        markerct = 1;
        for k = [1 2 3 4] 
            loglog(1./xnsize,squeeze(q_norm(:,i,k)),[markers(markerct) '-'],'linewidth',1,'markersize',markersize(markerct))

            if k == 1
                hold on
            end

            markerct = markerct+1;
        end
    end
    hold on
    if strcmp(IC_str,'_front')
        loglog(1./xnsize,squeeze(q_norma(:,i,1)),[markers(markerct) '-'],'linewidth',1,'markersize',markersize(markerct))
    end
    
    if i == 1
        

    if strcmp(IC_str,'_front')
        L=legend('Upwind','Lax-wendroff','Beam-Warming',...
            'Upwind FL','Upwind auto');
    elseif strcmp(IC_str,'_gauss')
        L=legend('Upwind','Lax-wendroff','beam warming','upwind FL');
    end

        set(L,'interpreter','latex','units','normalized','position',[.92 .6 .05 .15]);
        
    end
    
    if strcmp(IC_str,'_front')
        axis([10^-3.5 10^-1 10^-3.1 10^1])
    elseif strcmp(IC_str,'_gauss')
        axis([10^-3.5 10^-1.15 10^-5.5 10])
    end
    
    if i <= length(eta)
              title(['$\mathbf{\eta^2 = ',sprintf('%1.e',eta_str(i)^2) '}$'],...
                  'interpreter','latex','fontsize',15)
    end
            
    if mod(i,length(eta)) == 1
        ylabel(['$\mathbf{N = ' num2str(num2str(xnstr(i))) '}$'],'interpreter',...
            'latex','fontsize',15)
    end
            
    
    
    xlabel('$h$','interpreter','latex')
    %ylabel('$\| \theta_0 - \hat{\theta} \|_2$','interpreter','latex')
    xticks([10^-3 10^-2])
    
    
end


exportfig(gcf,['q_h_plot' IC_str '.eps'],'fontsize',1.5,'color','rgb')
saveas(gcf,['q_h_plot' IC_str '.fig'])

