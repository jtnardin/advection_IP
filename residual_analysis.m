clear all; clc

%%%%specify which initial condition, method,
%%%%and data set we're compating
IC_str = '_front';
num_meth = 3;
m = 15;%16;

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

xnstr = zeros(1,numel(data));
eta_str = zeros(1,numel(data));

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
num_meth_cell{5} = 'upwindfl';


model_sims = cell(length(xnsize)-1,1);
res = cell(length(xnsize),1);

for i = 1:7
    
    
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
    
    %_ols{i,m,num_meth}
    %get model sim
    [J,res{i},model_sims{i}] = MLE_cost_art_data(cell_data,...
        q0,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,...
        num_meth_cell{num_meth},t,tdata,model_str);
    
    %fix res structure
    res{i} = reshape(res{i},size(model_sims{1}));
    
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

% figure
% for i = 2:length(xnsize)
%     subplot(3,2,i-1)
%     
%     for j = 1:length(tdata)
%         hold on
%         plot(xdata,model_sims{i-1}(j,:),[colors(j) '-'])
%         plot(xdata,cell_data(j,:),[colors(j) '*'])
% 
% %         plot(xdata,res{i-1}(j,:),[colors(j) '*'])
%         
% %       
%     end
%     
% %     axis([0 1 -.6 .6])
%     
%     if i == 2
%         legend('1','2','3','4','5','6','location','northeast')
%     end
% end


for i = 6
    figure('unit','normalized','outerposition',[0 0 1 1])
    for j = 1:length(tdata)
        subplot(3,2,j)
        
        yyaxis left
                ax1 = gca;
                hold on
                plot(xdata,model_sims{i-1}(j,:),[colors(j) '-'])
    %             hold on
    %             plot(xdata,cell_data(j,:),[colors(j) '*'])
        ylabel('model')
        set(ax1,'ycolor','k')
        axis([0 1 -.1 5.1])
        
        yyaxis right
                ax2 = gca;
                hold on
                plot(xdata,res{i-1}(j,:),[colors(j) '*'])
                plot([xdata(1) xdata(end)],zeros(1,2),'-')
                axis([0 1 -1 1])
                
        
        set(ax2,'ycolor','k')
        xlabel('x')
        ylabel('model - data')
        title(['Residual, t = ' num2str(tdata(j))])

        
        if j == 1
            h=legend('model','residuals');
            if num_meth == 3
                set(h,'position',h.Position+[-0.05 0.01 0 0]);
            else
                set(h,'position',h.Position+[-0.03 .01 0 0]);
            end
        end
        
        
    end
end

exportfig(gcf,['residual' IC_str '_' num2str(i) '_' num2str(num_meth) '_noise_' num2str(eta_str(m)) '.eps'],'fontsize',2,'color','rgb')
saveas(gcf,['residual' IC_str '_' num2str(i) '_' num2str(num_meth) '_noise_' num2str(eta_str(m)) '.fig'])