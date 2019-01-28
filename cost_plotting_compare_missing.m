clear all; clc


IC_str = '_front';


%load best-fit params, data, and initial condition
if strcmp(IC_str,'_gauss')
    load(['advection_rates' IC_str '_IC_all_3_26.mat'])
    load(['advection_art_data' IC_str '_all_3_26.mat'])
elseif strcmp(IC_str,'_front')
%     load(['advection_art_data' IC_str '_all_rev.mat'])
%     load('advection_rates_autoreg_front_IC_rev.mat')
%     
    load(['advection_art_data' IC_str '_all_missing.mat'])
    load('advection_rates_front_IC_missing.mat')
    
    J_ols_missing = J_ols;
    q_ols_missing = q_ols;
    
    load(['advection_art_data' IC_str '_all_rev.mat'])
    load('advection_rates_autoreg_front_IC_rev.mat')
    
end
    
phi = IC_spec(IC_str(2:end));


xnsize = [21,41,81,161,321,641,2*640+1];

num_meth_cell = cell(5,1);
num_meth_cell{1} = 'upwind';
%num_meth_cell{2} = 'Lax-Friedrich';
num_meth_cell{2} = 'Lax-Wendroff';
num_meth_cell{3} = 'Beam-Warming';
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


q_norm = zeros(size(q_ols));
q_norm_missing = zeros(size(q_ols_missing));
for i = 1:size(q_ols,1)
    for j = 1:size(q_ols,2)
        for k = 1:size(q_ols,3)
            
            if ~isempty(q_ols{i,j,k}) && norm(q_ols_missing{i,j,k})~=0
                q_norm(i,j,k) = norm(q_ols{i,j,k}-q0);
                q_norm_missing(i,j,k) = norm(q_ols_missing{i,j,k}-q0);
            end
            
        end
    end
end


figure('units','normalized','outerposition',[0 0 1 1])


markers = '.^sxv';
markersize = [20 10 10 15 8];

eta_sel = 5;


subplot(2,2,1)
loglog(1./xnsize,squeeze(J_ols(:,eta_sel,1)),[markers(1) '-'],'linewidth',2,'markersize',markersize(1))
hold on
loglog(1./xnsize,squeeze(J_ols(:,eta_sel+length(eta),1)),[markers(2) '-'],'linewidth',2,'markersize',markersize(2))
loglog(1./xnsize,squeeze(J_ols_missing(:,eta_sel+length(eta),1)),[markers(3) '-'],'linewidth',2,'markersize',markersize(3))
loglog(1./xnsize,squeeze(J_ols_missing(:,eta_sel+2*length(eta),1)),[markers(4) '-'],'linewidth',2,'markersize',markersize(4))


plot([1e-16 1],repmat(eta_str(eta_sel)^2,1,2),'k--') 
    

subplot(2,2,3)

loglog(1./xnsize,squeeze(J_ols(:,eta_sel,3)),[markers(1) '-'],'linewidth',2,'markersize',markersize(1))
hold on
loglog(1./xnsize,squeeze(J_ols(:,eta_sel+length(eta),3)),[markers(2) '-'],'linewidth',2,'markersize',markersize(2))
loglog(1./xnsize,squeeze(J_ols_missing(:,eta_sel+length(eta),3)),[markers(3) '-'],'linewidth',2,'markersize',markersize(3))
loglog(1./xnsize,squeeze(J_ols_missing(:,eta_sel+2*length(eta),3)),[markers(4) '-'],'linewidth',2,'markersize',markersize(4))

plot([1e-16 1],repmat(eta_str(eta_sel)^2,1,2),'k--') 

subplot(2,2,2)

loglog(1./xnsize,squeeze(q_norm(:,eta_sel,1)),[markers(1) '-'],'linewidth',2,'markersize',markersize(1))
hold on
loglog(1./xnsize,squeeze(q_norm(:,eta_sel+length(eta),1)),[markers(2) '-'],'linewidth',2,'markersize',markersize(2))
loglog(1./xnsize,squeeze(q_norm_missing(:,eta_sel+length(eta),1)),[markers(3) '-'],'linewidth',2,'markersize',markersize(3))
loglog(1./xnsize,squeeze(q_norm_missing(:,eta_sel+2*length(eta),1)),[markers(4) '-'],'linewidth',2,'markersize',markersize(4))

subplot(2,2,4)

loglog(1./xnsize,squeeze(q_norm(:,eta_sel,3)),[markers(1) '-'],'linewidth',2,'markersize',markersize(1))
hold on
loglog(1./xnsize,squeeze(q_norm(:,eta_sel+length(eta),3)),[markers(2) '-'],'linewidth',2,'markersize',markersize(2))
loglog(1./xnsize,squeeze(q_norm_missing(:,eta_sel+length(eta),3)),[markers(3) '-'],'linewidth',2,'markersize',markersize(3))
loglog(1./xnsize,squeeze(q_norm_missing(:,eta_sel+2*length(eta),3)),[markers(4) '-'],'linewidth',2,'markersize',markersize(4))


for i = 1:2:3
    subplot(2,2,i) 
    axis([10^-3.5 1e-1 10^-(1.5) 10^-.5])
    if i == 1
        title(['Upwind Method, $\mathbf{\eta^2 = ',sprintf('%1.e',eta_str(eta_sel)^2) '}$'],...
                  'interpreter','latex','fontsize',15)
    elseif i == 3
        title(['Beam--Warming Method, $\mathbf{\eta^2 = ',sprintf('%1.e',eta_str(eta_sel)^2) '}$'],...
                  'interpreter','latex','fontsize',15)
    end
    ylabel('$J_{OLS}(h,\hat{\theta})$','interpreter',...
            'latex','fontsize',15)
    xlabel('$h$','interpreter','latex')
    xticks([10^-3 10^-2])
end

for i = 2:2:4
    subplot(2,2,i) 
    axis([10^-3.5 10^-1 10^-3 10^1])
    
    if i == 2
        title(['Upwind Method, $\mathbf{\eta^2 = ',sprintf('%1.e',eta_str(eta_sel)^2) '}$'],...
                  'interpreter','latex','fontsize',15)
    elseif i == 4
        title(['Beam-Warming Method, $\mathbf{\eta^2 = ',sprintf('%1.e',eta_str(eta_sel)^2) '}$'],...
                  'interpreter','latex','fontsize',15)
    end
    
    
    ylabel('$\|\hat{\theta}_{OLS}-\theta_0\|$',...
              'interpreter','latex','fontsize',15)
    
    
    xlabel('$h$','interpreter','latex')
    %ylabel('$\| \theta_0 - \hat{\theta} \|_2$','interpreter','latex')
    xticks([10^-3 10^-2])
end

subplot(2,2,1)
h=legend('N = 11','N = 31','N = 11, missing points','N = 31, missing points','location','northwest');
set(h,'interpreter','latex','fontsize',15);


exportfig(gcf,['J_h_missing_plot' IC_str '.eps'],'fontsize',2.5,'color','rgb')
saveas(gcf,['J_h_missing_plot' IC_str '.fig'])

