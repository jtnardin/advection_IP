%Used to check local minima from optimizer
%tell which scenario you want to sim, specify parameter, get resulting cost
%function

clear all;clc

%%%%specify which initial condition, method,
%%%%and data set we're compating
IC_str = '_front';
model_str = 'root';
num_meth = 1;
m = 3;
% xni =5;

%load best-fit params, data, and initial condition

if strcmp(IC_str,'_front')
    load(['advection_rates_autoreg' IC_str '_IC_all.mat'])
elseif strcmp(IC_str,'_gauss')
    load(['advection_rates' IC_str '_IC_all.mat'])
end

load(['advection_art_data' IC_str '.mat'])
phi = IC_spec(IC_str(2:end));



%grid sizes
xnsize = [21,41,81,161,321,641,2*640+1];
lambda = 1/2;

for i = 1:length(xd)
    xndata(i) = length(xd{i});
end

xdi = ceil(m/length(eta));
sigmaj = mod(m,length(eta));

if sigmaj == 0
    sigmaj = length(eta);
end

xnstr(m) = xndata(xdi);
eta_str(m) = eta(sigmaj);

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


num_meth_cell_title = cell(4,1);
num_meth_cell_title{1} = 'Upwind';
num_meth_cell_title{2} = 'LaxFried';
num_meth_cell_title{3} = 'LaxWend';
num_meth_cell_title{4} = 'BeamWarm';
num_meth_cell_title{5} = 'UpwindFL';



% alpha = linspace(q_ols{xni,m,num_meth}(1)-cwidth(k),q_ols{xni,m,num_meth}(1)+cwidth(k),10);
% beta = linspace(q_ols{xni,m,num_meth}(2)-cwidth(k),q_ols{xni,m,num_meth}(2)+cwidth(k),10);
% [Alpha,Beta] = meshgrid(alpha,beta);

count = 1;

model_sim = cell(2,5);

for xni = 3:4:7
    for num_meth = [1 3 4 5]

%     model_sim{count} = zeros(size(cell_data));
    
%     for i = 1:length(alpha)
%         for j = 1:length(beta)


            q = q_ols{xni,m,num_meth};



            %space
            x = linspace(0,1,xnsize(xni));
            dx = x(2) - x(1);
            xn = xnsize(xni);

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
            [J,res,~,model_sim{count,num_meth}] = MLE_cost_art_data(cell_data,...
                q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,...
                num_meth_cell{num_meth},t,tdata,model_str);

                                   
%         end
    end
    count = count + 1;
end

    %     figure
        % if mod(m,4)
        %     C = contour(Alpha,Beta,J_cont,[eta(m)^2*[1 2 4] .1:.1:max(max(J_cont))]);
        % else
        %     C = contour(Alpha,Beta,J_cont,[1e-3 1e-2 .1:.1:max(max(J_cont))]);
        % end

colors = 'rgbckm';
markers = '*sx^vo';
figure('units','normalized','outerposition',[0 0 1 1])
count = 1;
for i = 1:2
    for j = [1 3 4 5]
    subplot(2,4,count)
    
        hold on
        for k = 1:6
            plot(xdata,cell_data(k,:),[colors(k) markers(k)])%,'markersize',15)
            plot(linspace(0,1,xnsize(4*(i-1)+3)),model_sim{i,j}(:,k),[colors(k) '-'],'markersize',15)
        end
        
        %title({[num_meth_cell{j} ' method'],['$h = ' num2str(1./(xnsize(2*i+1)-1)) '$']}...
        %     ,'interpreter','latex')
        line1 = [num_meth_cell_title{j} ' method'];
        line2 = ['$h = ' num2str(1./(xnsize(4*(i-1)+3)-1)) '$'];
        title(sprintf('\\begin{tabular}{c} %s %s %s \\end{tabular}',line1,'\\',line2),...
        'interpreter','latex') 
         
         
        xlabel('$x$','interpreter','latex')
        ylabel('$u$','interpreter','latex')
        
        count = count+1;
        
        axis([0 .8 -.5 5.5])
        
    end
end

exportfig(gcf,['plots_v_data' IC_str '.eps'],'color','rgb','fontsize',2)
saveas(gcf,['plots_v_data' IC_str '.fig'])
