clear all; clc

%%%%specify which initial condition, method,
%%%%and data set we're compating
IC_str = '_gauss';
num_meth = 1;
m = 2;

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

%grid sizes
xnsize = [21,41,81,161,321,641,2*640+1];
lambda = 1/2;

xndata = [10, 50];
eta = [0 0.05];

xnstr = zeros(1,4);
eta_str = zeros(1,4);


%indices
xdi = ceil(m/2);
sigmaj = mod(m,2);

if sigmaj == 0
    sigmaj = 2;
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


model_sims = cell(length(xnsize)-1,1);
res = cell(length(xnsize)-1,1);

for i = 2:length(xnsize)
    
    
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
    [J,res{i-1},model_sims{i-1}] = MLE_cost_art_data(cell_data,...
        q{i,m,num_meth},dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xdata,...
        num_meth_cell{num_meth},t,tdata);
    
    %fix res structure
    res{i-1} = reshape(res{i-1},size(model_sims{1}));
    
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
    
colors = 'bgrkmc';

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


i = 4;
figure
for j = 1:length(tdata)
    subplot(3,2,j)
    yyaxis left
            hold on
            plot(xdata,model_sims{i-1}(j,:),[colors(j) '-'])
%             hold on
%             plot(xdata,cell_data(j,:),[colors(j) '*'])
        
        yyaxis right
            hold on
            plot(xdata,res{i-1}(j,:),[colors(j) '*'])
       
        
    
end
