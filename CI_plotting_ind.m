clear all; clc

IC_str = '_front';
stat_meth = '_autor';

num_meth=1;


load(['CI' IC_str stat_meth '_all.mat'])
load(['advection_art_data' IC_str '_all.mat'])

xnsize = [21,41,81,161,321,641,2*640+1];

num_meth_cell = cell(4,1);
num_meth_cell{1} = 'Upwind';
num_meth_cell{2} = 'Lax-Friedrich';
num_meth_cell{3} = 'Lax-Wendroff';
num_meth_cell{4} = 'Beam warming';


num_meth_short_cell = cell(4,1);
num_meth_short_cell{1} = 'UW';
num_meth_short_cell{2} = 'LF';
num_meth_short_cell{3} = 'LW';
num_meth_short_cell{4} = 'BW';
num_meth_short_cell{5} = 'UWFL';

for i = 1:length(xd)
    xndata(i) = length(xd{i});
end
for i = 1:length(eta)
    eta_vec(i) = eta(i);
end

i = 4;
xnstr = zeros(1,numel(data));
eta_str = zeros(1,numel(data));


for m = 1:numel(data)
    

    xdi = ceil(m/length(eta_vec));
    sigmaj = mod(m,length(eta_vec));
    
    if sigmaj == 0
        sigmaj = length(eta_vec);
    end
    
    xnstr(m) = xndata(xdi);
    eta_str(m) = eta_vec(sigmaj);
    
end


% figure
%%make video?
title_m = ['CI_advec_sim' num_meth_short_cell{num_meth} stat_meth '_' num2str(i) '.avi'];
f1 = figure();
vid = VideoWriter(title_m); %%title here
vid.Quality = 100;
vid.FrameRate = 2;
open(vid);



for j = 1

    
    
        hold on
%         plot(log2(q0(1)),log2(q0(2)),'k*')
        plot(q0(1),q0(2),'k*')
        
        ax_plot = gca;


        xlabel('$\alpha$','interpreter','latex','fontsize',15)
        ylabel('$\beta$','interpreter','latex','fontsize',15)

        if strcmp(stat_meth,'_OLS')
            title({strcat('95\% OLS CI,  ',num_meth_short_cell{j},' method') ; strcat(' $N$ = ',num2str(xnstr(i)),', $\eta^2 = $',...
                    num2str(eta_str(i)^2))},'interpreter','latex','fontsize',15)
        elseif strcmp(stat_meth,'_autor')
            title({strcat('95\% auto CI,  ',num_meth_short_cell{j},' method') ; strcat(' $N$ = ',num2str(xnstr(i)),', $\eta^2 = $',...
                num2str(eta_str(i)^2))},'interpreter','latex','fontsize',15)
        end
            
        caxis([1/xnsize(end) 1/xnsize(1)])
        c=colorbar('units','normalized','position',[.8 .4 .04 .25]);
        xns = length(xnsize);
        set(c,'colormap',[1-(1/(xns):(1/(xns)):1)' zeros(xns,1) (1/(xns):1/(xns):1)'])
        set(c,'ticks',[]);
        ax = axes('Position',[0 0 1 1],'visible','off');
        axes(ax)

        a=text(.8125,.67,'$h$','units','normalized','interpreter','latex','fontsize',10);
        text(.84,.635,'$(10\cdot2^0)^{-1}$','units','normalized','interpreter','latex','fontsize',10)
        text(.84,.52,'$(10\cdot2^3)^{-1}$','units','normalized','interpreter','latex','fontsize',10)
        text(.84,.415,'$(10\cdot2^6)^{-1}$','units','normalized','interpreter','latex','fontsize',10)

        axes(ax_plot)
            
        if strcmp(IC_str,'_front')
%             if j == 1
%                 axis([-3 -1.5 -1.1 -.7])
%             elseif j == 3
%                 axis([-2.5 .5 -1.7 -.7])
%             elseif j == 4
%                 axis([-2.7 -1.5 -1.3 -.8])
%             end
            axis([.12 .45 .44 .62])
        elseif strcmp(IC_str,'_gauss')
            if j == 1
                if mod(i,length(eta)) <= 5
                    axis([-3.5 -1.5 -1.4 -.75])
                else
                    axis([-5 0 -2 0])
                end
            elseif j == 2
                axis([-2 -1.5 -1.35 -1.26])
            elseif j ==3
                if mod(i,length(eta))<=5
                    axis([-2.1 -1 -1.5 -1.2])
                elseif mod(i,length(eta))>=6
                    axis([-3.2 -.7 -1.8 -.9])
                elseif mod(i,length(eta))==0
                    axis([-5 -1 -1.5 -.8])
                end
            elseif j == 4
                if mod(i,length(eta))>=6
                    axis([-3.2 -.7 -1.6 -1])
                elseif mod(i,length(eta)) == 0
                    axis([-4 1 -2.5 -.5])
                else
                    axis([-2.4 -1.65 -1.4 -1.15])
                end
            end
        end
        
        
        
%     exportfig(gcf,['CI_h_plot' IC_str '_' num2str(j) stat_meth '.eps'],'fontsize',1.15,'color','rgb')
%     saveas(gcf,['CI_h_plot' IC_str '_' num2str(j) stat_meth '.fig'])

    axes(ax)
    writeVideo(vid, getframe(f1));
    axes(ax_plot)
    for l = 1:7
            C = CI{l,i,j};
            if (~any(isnan(C(:))))&&(~isempty(C))
%                 rectangle('position',[log2(C(1,1)) log2(C(2,1)) log2(C(1,2))-log2(C(1,1)) log2(C(2,2))-log2(C(2,1))],...
                  rectangle('position',[C(1,1) C(2,1) C(1,2)-C(1,1) C(2,2)-C(2,1)],...
                    'edgecolor',[(l-0)/7 0 1-(l-0)/7],'linewidth',3)%repmat(1-(l+1)/8,1,3))
            end
            
%             pause(1)
            axes(ax)
            writeVideo(vid, getframe(f1));
            axes(ax_plot)
    end

    
end

close(vid)
