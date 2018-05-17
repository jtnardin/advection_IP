clear all; clc


IC_str = '_front';

numel_data =8;

%load best-fit params, data, and initial condition
if strcmp(IC_str,'_gauss')
     load(['advection_rates' IC_str '_IC_all.mat'])
    load(['advection_art_data' IC_str '_all.mat'])
elseif strcmp(IC_str,'_front')
    load(['advection_rates_autoreg' IC_str '_IC_all.mat'])
    load(['advection_art_data' IC_str '_all.mat'])

end
   
    phi = IC_spec(IC_str(2:end));


xnsize = [21,41,81,161,321,641,2*640+1];

num_meth_cell = cell(4,1);
num_meth_cell{1} = 'upwind';
num_meth_cell{2} = 'Lax-Friedrich';
num_meth_cell{3} = 'Lax-Wendroff';
num_meth_cell{4} = 'Beam-Warming';
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


figh = loglog(1./xnsize,squeeze(J_ols(:,numel_data,[1 3 4])),'.-','linewidth',4,'markersize',30);

hold on

plot([1e-16 1],repmat(eta_str(numel_data)^2,1,2),'k--') 


if strcmp(IC_str,'_front')
    h=legend('Upwind','Lax-wendroff','Beam-Warming',...
        '$\eta^2$','location','southeast');
elseif strcmp(IC_str,'_gauss')
    h=legend('Upwind','Lax-wendroff','beam warming',...
        '$\eta^2$','location','northwest');
end
set(h,'interpreter','latex','units','normalized','position',[.81 .19 .05 .15]);


if strcmp(IC_str,'_gauss')
    vloc = [10^-7.5, 10^-8.5, 10^-9.5];
elseif strcmp(IC_str,'_front')
    vloc = [10^-1.25, 10^-1.35, 10^-1.45];
end
for i = 1:3
   
    if strcmp(IC_str,'_gauss')
        if i == 1 
            p=polyfit(log(1./xnsize(1:5))',log(J_ols(1:5,numel_data,i)),1);
        elseif i == 2
            p=polyfit(log(1./xnsize(1:3))',log(J_ols(1:3,numel_data,i+1)),1);
        elseif i == 3
            p=polyfit(log(1./xnsize(1:4))',log(J_ols(1:4,numel_data,i+1)),1);
        end
    elseif strcmp(IC_str,'_front')
       if i == 1 
            p=polyfit(log(1./xnsize(1:7))',log(J_ols(1:7,numel_data,i)),1);
        elseif i == 2
            p=polyfit(log(1./xnsize(1:7))',log(J_ols(1:7,numel_data,i+1)),1);
        elseif i == 3
            p=polyfit(log(1./xnsize(1:7))',log(J_ols(1:7,numel_data,i+1)),1);
        end
 
    end
   text(10^-2,vloc(i),sprintf('%.3f',p(1)),'fontsize',18,'color',figh(i).Color)
end
    
      
    
if strcmp(IC_str,'_front')
%         axis([1e-4 1e-1 1e-4 1e0])
      axis([10^-3.5 10^-1.2 10^-1.6 10^-.4])
elseif strcmp(IC_str,'_gauss')
    axis([10^-3.3 10^-1.2 10^-10 1])
end


title(strcat('Numerical Cost Function, $N$ =  ',num2str(xnstr(numel_data)),', $\eta^2 = $ ',...
        num2str(eta_str(numel_data)^2)),'interpreter','latex','fontsize',15)

xlabel('$h$','interpreter','latex','fontsize',15)
ylabel('$J(h,\hat{\theta})$','interpreter','latex','fontsize',15)



% exportfig(gcf,['J_h_plot' IC_str '_ind.eps'],'fontsize',3.5,'color','rgb')
% saveas(gcf,['J_h_plot' IC_str '_ind.fig'])


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

pu = polyfit(log(1./(xnsize(1:6)-1)),log(q_norm(1:6,numel_data,1))',1);
pu = pu(1)


pLW = polyfit(log(1./(xnsize(1:6)-1)),log(q_norm(1:6,numel_data,3))',1);
pLW = pLW(1)

pBW = polyfit(log(1./(xnsize(1:4)-1)),log(q_norm(1:4,numel_data,4))',1);
pBW= pBW(1)

loglog(1./xnsize,squeeze(q_norm(:,numel_data,[1 3 4])),'.-','linewidth',4,'markersize',40)
hold on
if strcmp(IC_str,'_front')
    loglog(1./xnsize,squeeze(q_norma(:,numel_data,1)),'.-','linewidth',4,'markersize',40)
end


if strcmp(IC_str,'_front')
    h=legend('Upwind','Lax-wendroff','Beam-Warming','Upwind Autocorrelative');
elseif strcmp(IC_str,'_gauss')
    h=legend('Upwind','Lax-wendroff','beam warming');
end

set(h,'interpreter','latex','units','normalized','position',[.81 .19 .05 .15]);



if strcmp(IC_str,'_front')
    axis([10^-3.3 10^-1 10^-3 10^0])
elseif strcmp(IC_str,'_gauss')
    axis([10^-3.3 10^-1.15 10^-5 10])
end

%     title(['$\| \hat{\theta} - \theta_0 \|_2$, ' num_meth_cell{i}],'interpreter','latex')



title(strcat('$\| \theta_0 - \hat{\theta} \|_2,$ ',' $N$ =  ',num2str(xnstr(numel_data)),', $\eta^2 = $',...
        num2str(eta_str(numel_data)^2)),'interpreter','latex','fontsize',15)

xlabel('$h$','fontsize',15,'interpreter','latex')
ylabel('$\| \cdot \|_2$','interpreter','latex','fontsize',15)





exportfig(gcf,['q_h_plota' IC_str '_ind.eps'],'fontsize',3.5,'color','rgb')
saveas(gcf,['q_h_plota' IC_str '_ind.fig'])

