%Used to check local minima from optimizer
%tell which scenario you want to sim, specify parameter, get resulting cost
%function

clear all;clc

%%%%specify which initial condition, method,
%%%%and data set we're compating
IC_str = '_gauss';
num_meth = 4;
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

alpha   = q0(1);
beta    = q0(2);


%rate of advection
[g,sigma,sigma_inv] = advection_rate('root',alpha,beta);
%initial condition
%final solution form.
soln = @(t,x) (g(x)~=0).*(x>=sigma_inv(t,0)).*g(sigma_inv(-t,x))./g(x).*phi(sigma_inv(-t,x));



%grid sizes
xnsize = [21,41,81,161,321,641,2*640+1];
lambda = 1/2;

xndata = [length(xd{1}), length(xd{2})];


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

xf = linspace(0,1,500);
tf = linspace(0,10,40);

count = 1;

model_sim = cell(1,3);

for xni = [4 5 6]
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
            [J,res,~,model_sim{count}] = MLE_cost_art_data(ones(length(xf),length(tf)),...
                q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xf,...
                num_meth_cell{num_meth},t,tf);

                                   
%         end
    count = count + 1;
end

xf = linspace(0,1,500);

colors = [0 0 1 ; .5 0 .5 ; 1 0 0];        
figure('units','normalized','outerposition',[0 0 1 1])


% %%%make video?
% title_m = ['advec_sim' num_meth_cell{num_meth} '.avi'];
% f1 = figure();
% vid = VideoWriter(title_m); %%title here
% vid.Quality = 100;
% vid.FrameRate = 7;
% open(vid);

set(gcf,'color',[1 1 1])

for k=1:length(tf)
    
    hold off
    g=plot(xf,soln(tf(k),xf),'k-','linewidth',1.5);
    hold on

    
    for i = 1:3

        hold on
%         plot(linspace(0,1,xnsize(2*i+1)),model_sim{i}(:,k),'color',colors(i,:),'linewidth',1.5)
plot(linspace(0,1,xnsize(i+3)),model_sim{i}(:,k),'color',colors(i,:),'linewidth',1.5)

    end
    
    title([num_meth_cell{num_meth} ' method' ],'interpreter','latex','fontsize',15)
    xlabel('$x$','interpreter','latex')
    ylabel('$u$','interpreter','latex')

    
    h=legend('$u_0(t,x)$','$u(t,x|h)$','$u(t,x|h/2)$','$u(t,x|h/4)$');
    set(h,'interpreter','latex','fontsize',15)
    
    uistack(g,'top')
    
    axis([0 .6 -.5 1.25])
    
    pause(.125)
%     writeVideo(vid, getframe(f1));
        
end

% close(vid)
% exportfig(gcf,['plots_v_data' IC_str '.eps'],'color','rgb','fontsize',1.5)
% saveas(gcf,['plots_v_data' IC_str '.fig'])
