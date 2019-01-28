clear all; clc

for k = 2

    for num_meth = [1 3] %[1 3 ]

        
        if k == 1
            IC_str = '_gauss';
        elseif k == 2
            IC_str = '_front';
        end

            model_str = 'root';


        clear eta eta_vec
        %load best-fit params, data, and initial condition
        if strcmp(IC_str,'_gauss')
            load(['advection_rates' IC_str '_IC_all.mat'])
            load(['advection_art_data' IC_str '_all.mat'])
        elseif strcmp(IC_str,'_front')
            load(['advection_rates_autoreg' IC_str '_IC_all.mat'])
            load(['advection_art_data' IC_str '_all.mat'])

        end

        %rate of advection
        [g,sigma,sigma_inv] = advection_rate(model_str,q0(1),q0(2));
        %initial condition
        phi = IC_spec(IC_str(2:end));

        %final solution form.
        soln = @(t,x) (g(x)~=0).*(x>=sigma_inv(t,0)).*g(sigma_inv(-t,x))./g(x).*phi(sigma_inv(-t,x));


        xnsize = [21,41,81,161,321,641,2*640+1];


        num_meth_cell = cell(5,1);
        num_meth_cell{1} = 'upwind';
        num_meth_cell{2} = 'laxfried';
        num_meth_cell{3} = 'laxwend';
        num_meth_cell{4} = 'beamwarm';
        num_meth_cell{5} = 'upwindfl';


        num_meth_short_cell = cell(5,1);
        num_meth_short_cell{1} = 'Upwind';
        num_meth_short_cell{2} = 'Lax-Friedrichs';
        num_meth_short_cell{3} = 'Lax-Wendroff';
        num_meth_short_cell{4} = 'Beam-Warming';
        num_meth_short_cell{5} = 'Upwind w/ flux limiters';



        for i = 1:length(xd)
            xndata(i) = length(xd{i});
        end
        for i = 1:length(eta)
            eta_vec(i) = eta(i);
        end

        %computation constants
        lambda = 1/2;

        %create grids for computaiton
        xnsize = [21,41,81,161,321,641,2*640+1];

        figure('units','normalized','outerposition',[0 0 1 1])
 
        
        markers = '.^sxv*d';
        markersize = [25 10 10 15 10 10 10];

        
        count = 1;
        
        for j = [1 2 4 8 9 11]%[1 4 7 15 18 21]

            xdi = ceil(j/length(eta_vec));
            sigmaj = mod(j,length(eta_vec));

            if sigmaj == 0
                sigmaj = length(eta_vec);
            end

            xnstr = xndata(xdi);
            eta_str = eta_vec(sigmaj);

            [Td,Xd] = meshgrid(td(2:end),xd{xdi});

            u0 = soln(Td,Xd);

            if strcmp(IC_str,'_front')
                u0(isnan(u0))=phi(Xd(isnan(u0)));
            elseif strcmp(IC_str,'_gauss')
                u0(isnan(u0))=0;
            end

            u0 = u0';

            epsilon = data{xdi,sigmaj}(2:end,:)-u0;

            A_J = 1/(numel(Td))*sum(epsilon(:).^2);
            B_J = zeros(7,1);
            C_J = zeros(7,1);
            D_J = zeros(7,1);
            E_J = zeros(7,1);
            F_J = zeros(7,1);
            
            G_J = zeros(7,1);

            for i = 1:7
                q = q_ols{i,j,num_meth};

                %rate of advection for qhat
                [gh,sigmah,sigma_invh] = advection_rate(model_str,q(1),q(2));
                %final solution form for qhat
                solnh = @(t,x) (g(x)~=0).*(x>=sigma_invh(t,0)).*gh(sigma_invh(-t,x))./gh(x).*phi(sigma_invh(-t,x));
                %analytical solution at qhat
                u0theta_hat = solnh(Td,Xd);
                u0theta_hat(isnan(u0theta_hat))=phi(Xd(isnan(u0theta_hat)));
                u0theta_hat = u0theta_hat';

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
                [J,res,uh,~] = MLE_cost_art_data(data{xdi,sigmaj},...
                    q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A,Abd,x,xd{xdi},...
                    num_meth_cell{num_meth},t,td,model_str);
                
                uh = uh(2:end,:);

                B_J(i) = 1/numel(Td)*sum((u0(:) - u0theta_hat(:)).^2);
                C_J(i) = 1/numel(Td)*sum((u0theta_hat(:) - uh(:)).^2);
                D_J(i) = 2/numel(Td)*sum(epsilon(:).*(u0(:)-u0theta_hat(:)));
                E_J(i) = 2/numel(Td)*sum(epsilon(:).*(u0theta_hat(:)-uh(:)));
                F_J(i) = 2/numel(Td)*sum((u0theta_hat(:)-uh(:)).*(u0(:) - u0theta_hat(:)));
            end

            J_final_abs = abs(A_J) + abs(B_J) + abs(C_J) + abs(D_J) + abs(E_J) + abs(F_J);
            J_final = A_J + B_J + C_J + D_J + E_J + F_J;


          

%             figure(fig_log)
            subplot(2,3,count)
            h = 1./(xnsize-1);
            p1 = loglog(h,J_final,[markers(1) '-'],'markersize',markersize(1),'linewidth',3);
            hold on
            p2 = loglog(h,abs(A_J*ones(7,1)),[markers(2) '-'],'markersize',markersize(2),'linewidth',2.5);
            p3 = loglog(h,abs(B_J),[markers(3) '-'],'markersize',markersize(3),'linewidth',2.5);
            p4 = loglog(h,abs(C_J),[markers(4) '-'],'markersize',markersize(4),'linewidth',2.5);
            p5 = loglog(h,abs(D_J),[markers(5) '-'],'markersize',markersize(5),'linewidth',2.5);
            p6 = loglog(h,abs(E_J),[markers(6) '-'],'markersize',markersize(6),'linewidth',2.5);
            p7 = loglog(h,abs(F_J),[markers(7) '-'],'markersize',markersize(7),'linewidth',2.5);
                        
            
            
            [J_shifth,A_shifth,B_shifth,C_shifth,D_shifth,E_shifth,F_shifth,...
                J_shiftv,A_shiftv,B_shiftv,C_shiftv,D_shiftv,E_shiftv,F_shiftv] = ...
                shift_det_f(num_meth,count);
            
            
            text(h(end)-2.5e-4-J_shifth, J_final(end)-J_shiftv, 'J', 'HorizontalAlignment','center', 'VerticalAlignment','middle','fontsize',10,'color',p1.Color)
            text(h(end)-2.5e-4-A_shifth, A_J-A_shiftv, 'A', 'HorizontalAlignment','center', 'VerticalAlignment','middle','fontsize',10,'color',p2.Color)
            text(h(end)-2.5e-4-B_shifth, B_J(end)-B_shiftv, 'B', 'HorizontalAlignment','center', 'VerticalAlignment','middle','fontsize',10,'color',p3.Color)
            text(h(end)-2.5e-4-C_shifth, C_J(end)-C_shiftv, 'C', 'HorizontalAlignment','center', 'VerticalAlignment','middle','fontsize',10,'color',p4.Color)
            text(h(end)-2.5e-4-D_shifth, abs(D_J(end))-D_shiftv, 'D', 'HorizontalAlignment','center', 'VerticalAlignment','middle','fontsize',10,'color',p5.Color)
            text(h(end)-2.5e-4-E_shifth, abs(E_J(end))-E_shiftv, 'E', 'HorizontalAlignment','center', 'VerticalAlignment','middle','fontsize',10,'color',p6.Color)
            text(h(end)-2.5e-4-F_shifth, abs(F_J(end))-F_shiftv, 'F', 'HorizontalAlignment','center', 'VerticalAlignment','middle','fontsize',10,'color',p7.Color)
                        
            if count == 3
                    h=legend('J','A','B','C','D','E','F');
                    set(h,'units','normalized','position',[.93,.45,.03,.1])
            end

            uistack(p1,'top')
            
              title(['$N$ = ',num2str(xndata(xdi)) ', $\eta^2 =$ ' ...
                num2str(eta_vec(sigmaj)^2)],'interpreter','latex')

            xlabel('$h$','interpreter','latex')
            
            if strcmp(IC_str,'_gauss')
                if num_meth == 1
                    axis([10^-3.5 10^-1 10^-9 1 ])
                elseif num_meth == 3 || num_meth == 4
                    axis([10^-3.5 10^-1 10^-12 1 ])
                end
                    
            elseif strcmp(IC_str,'_front')
                if num_meth == 1
                    axis([10^-3.75 10^-1 10^-5 10^.5 ])
                elseif num_meth == 3 || num_meth == 4
                    axis([10^-3.75 10^-1 10^-6 10^.5 ])
                end
            end

            xticks(10.^(-3:-1))
            yticks(10.^(-5:0))
            
            count = count + 1;
        end

        annotation('textbox',[.4 .99 .4 .03],'string',['J components, '...
            num_meth_short_cell{num_meth},' method'],'units','normalized','edgecolor','none','interpreter','latex')
        
        
        
        exportfig(gcf,['J_comp' IC_str '_' num2str(num_meth) '_manuscript.eps'],'fontsize',2.5,'color','rgb')
        saveas(gcf,['J_comp' IC_str '_' num2str(num_meth) '_manuscript.fig'])
        
    end
end