clear all; clc

%%%%specify which initial condition, method,
%%%%and data set we're compating
IC_str = '_front';
num_meth = 1;
m = 16;

model_str = 'root';


%load best-fit params, data, and initial condition
load(['advection_rates_autoreg' IC_str '_IC_all.mat'])
load(['advection_art_data' IC_str '_all.mat'])
phi = IC_spec(IC_str(2:end));


%grid sizes
xnsize = [21,41,81,161,321,641,2*640+1];
lambda = 1/2;



for i = 1:length(xd)
    xndata(i) = length(xd{i});
end
for i = 1:length(eta)
    eta_vec(i) = eta(i);
end

xnstr = zeros(1,4);
eta_str = zeros(1,4);

xdi = ceil(m/length(eta_vec));
sigmaj = mod(m,length(eta_vec));

if sigmaj == 0
    sigmaj = length(eta_vec);
end

xnstr(m) = xndata(xdi);
eta_str(m) = eta_vec(sigmaj);

%select exp. data
cell_data = data{xdi,sigmaj};
tdata = td;
xdata = xd{xdi};

num_meth_cell = cell(4,1);
num_meth_cell{1} = 'upwind';
num_meth_cell{2} = 'laxfried';
num_meth_cell{3} = 'laxwend';
num_meth_cell{4} = 'beamwarm';


model_sims_ols = cell(length(xnsize),1);
model_sims_auto = cell(length(xnsize),1);
mod_res = cell(length(xnsize),1);
res = cell(length(xnsize),1);

for i = 1:length(xnsize)
    
    
    %space
    x = linspace(0,1,xnsize(i));
    dx = x(2) - x(1);
    xn = xnsize(i);

    %time
    dt = lambda*dx;
    tfin = 10;
    t = 0:dt:tfin;
    tn = length(t);

    %points for computaiton
    [x_int,xbd_0,xbd_1] = int_bd_def(xn);

    %initial condition
    IC = phi(x);

    
    %load matrices for computation
    [A,Abd] = aMatrixupwind(xn,num_meth_cell{num_meth});
    
    
    %get model sim
    [J,res{i},model_sims_ols{i}] = MLE_cost_art_data(cell_data,...
        q_ols{i,m,num_meth},dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,...
        num_meth_cell{num_meth},t,tdata,model_str);
    
    %fix res structure
    res{i} = reshape(res{i},length(tdata),length(xdata));
    
%     [phi1{i,m,num_meth},phi2{i,m,num_meth}] = autoregression_terms_compute(res,xdata);
    
    
    %get model sim
    [J,mod_res{i},model_sims_auto{i}] = MLE_cost_autoreg_art_data(cell_data,...
        q_autoreg{i,m,num_meth},dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,...
        num_meth_cell{num_meth},t,tdata,phi1{i,m,num_meth},phi2{i,m,num_meth});
    
    %fix res structure
    mod_res{i} = reshape(mod_res{i},length(tdata),length(xdata));
    
end

% 
% figure
% for i = 2:length(xnsize)
%     subplot(3,2,i-1)
%     
%     contourf(xdata,tdata,res{i-1},'edgecolor','none')
%     colorbar
%     
% end
    
colors = 'rgbckm';


for i = 4
    
    for l = 1:2
        figure('unit','normalized','outerposition',[0 0 1 1])
        for j = 1:length(tdata)
            subplot(3,2,j)

            yyaxis left
                    hold on
                    if l == 1
                        plot(xdata,model_sims_ols{i}(j,:),[colors(j) '-'])
                    elseif l == 2
                        plot(xdata,model_sims_auto{i}(j,:),[colors(j) '-'])
                    end
                    
            ylabel('model')
            set(gca,'ycolor','k')
            
            yyaxis right
                    hold on
                    if l == 1
                        plot(xdata,res{i}(j,:),[colors(j) '*'])
                    elseif l == 2
                        plot(xdata,mod_res{i}(j,:),[colors(j) '*'])
                    end
                    plot([xdata(1) xdata(end)],zeros(1,2),'-')
                    axis([0 1 -1 1])

            xlabel('x')
            ylabel('model - data')
            if l == 1
                title(['OLS Residual, t = ' num2str(tdata(j))])
            elseif l == 2
                title(['Modified Residual, t = ' num2str(tdata(j))])
            end
            set(gca,'ycolor','k')


            if j == 1
                h=legend('model','residuals');
                set(h,'position',h.Position+[-0.0425 .01 0 0]);
            end


        end
    

        if l == 1
            exportfig(gcf,['residual_' IC_str '_' num2str(i) '_noise_' num2str(eta_str(m)) '.eps'],'fontsize',2,'color','rgb')
            saveas(gcf,['residual_' IC_str '_' num2str(i) '_noise_' num2str(eta_str(m)) '.fig'])
        elseif l == 2
            exportfig(gcf,['mod_residual_' IC_str '_' num2str(i) '_noise_' num2str(eta_str(m)) '.eps'],'fontsize',2,'color','rgb')
            saveas(gcf,['mod_residual_' IC_str '_' num2str(i) '_noise_' num2str(eta_str(m)) '.fig'])
        end
        
    end
    

end


% figure('unit','normalized','outerposition',[0 0 1 1])
% 
% i = 4;
% j = 4;
% 
% yyaxis left
% hold on
% 
% plot(xdata,model_sims_auto{i}(j,:),[colors(j) '-'])
% ylabel('model','fontsize',15)
% axis([0 1 0 5.1])
% 
% yyaxis right
% hold on
% plot(xdata,mod_res{i}(j,:),[colors(j) '.'],'markersize',30)
% plot([xdata(1) xdata(end)],zeros(1,2),'-','linewidth',3)
% axis([0 1 -1.2 1.2])
% 
% xlabel('$x$','interpreter','latex','fontsize',15)
% ylabel('model - data','fontsize',15)
% title(['Modified Residual, t = ' num2str(tdata(j))],'fontsize',15,'interpreter','latex')
% 
% legend('model','residuals')
%     
% exportfig(gcf,['mod_residual_' IC_str '_' num2str(i) '_noise_ind.eps'],'fontsize',3.5,'color','rgb')
% saveas(gcf,['mod_residual_' IC_str '_' num2str(i) '_noise_ind.fig'])
       