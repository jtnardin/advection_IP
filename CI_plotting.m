clear all; clc

IC_str = '_front';
stat_meth = '_autor';

if strcmp(IC_str,'_front')
    load(['CI' IC_str stat_meth '_all.mat'])
elseif strcmp(IC_str,'_gauss')
    load(['CI' IC_str stat_meth '_all.mat'])
end

if strcmp(IC_str,'_gauss')
    load(['advection_art_data' IC_str '_all_3_26.mat'])
elseif strcmp(IC_str,'_front')
    load(['advection_art_data' IC_str '_all.mat'])
end

xnsize = [21,41,81,161,321,641,2*640+1];

if strcmp(IC_str,'_front')
    num_meth_cell = cell(4,1);
    num_meth_cell{1} = 'Upwind';
    num_meth_cell{2} = 'Lax-Friedrich';
    num_meth_cell{3} = 'Lax-Wendroff';
    num_meth_cell{4} = 'Beam warming';
    num_meth_cell{5} = 'Upwind FL';


    num_meth_short_cell = cell(4,1);
    num_meth_short_cell{1} = 'UW';
    num_meth_short_cell{2} = 'LF';
    num_meth_short_cell{3} = 'LW';
    num_meth_short_cell{4} = 'BW';
    num_meth_short_cell{5} = 'UWFL';
elseif strcmp(IC_str,'_gauss')
    num_meth_cell = cell(4,1);
    num_meth_cell{1} = 'Upwind';
    num_meth_cell{2} = 'Lax-Wendroff';
    num_meth_cell{3} = 'Beam warming';
    num_meth_cell{4} = 'Upwind FL';


    num_meth_short_cell = cell(4,1);
    num_meth_short_cell{1} = 'UW';
    num_meth_short_cell{2} = 'LW';
    num_meth_short_cell{3} = 'BW';
    num_meth_short_cell{4} = 'UWFL';

end

for i = 1:length(xd)
    xndata(i) = length(xd{i});
end
for i = 1:length(eta)
    eta_vec(i) = eta(i);
end


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





for j = [1]%3 4 5]

    figure('units','normalized','outerposition',[0 0 1 1])
    
    for i = 1:numel(data)

        subplot(size(data,1),size(data,2),i)
        hold on
        plot(q0(1),q0(2),'k*')



        for l = 1:7
            C = CI{l,i,j};
            if (~any(isnan(C(:))))&&(~isempty(C))
%                 rectangle('position',[log2(C(1,1)) log2(C(2,1)) log2(C(1,2))-log2(C(1,1)) log2(C(2,2))-log2(C(2,1))],...
%                     'edgecolor',[(l-0)/7 0 1-(l-0)/7])%repmat(1-(l+1)/8,1,3))
                rectangle('position',[C(1,1) C(2,1) C(1,2)-C(1,1) C(2,2)-C(2,1)],...
                     'edgecolor',[(l-0)/7 0 1-(l-0)/7])%repmat(1-(l+1)/8,1,3))
            end
        end
        
        xlabel('$\alpha$','interpreter','latex')
        ylabel('$\beta$','interpreter','latex')

         title({['95\% CI, ',num_meth_short_cell{j} ' method'] ; ['$N$ = '...
             num2str(xnstr(i)) ', $\eta^2 = ' sprintf('%1.e',eta_str(i)^2) '$']},'interpreter','latex')


        if strcmp(IC_str,'_front')
            if j == 1
                if strcmp(stat_meth,'_OLS')
                    if mod(i,length(eta)) ~= 0 
                        axis([0 .37 .46 .7])
                    else
                        axis([0 .52 .4 .775])
                    end
                elseif strcmp(stat_meth,'_autor')
                    if mod(i,length(eta)) ~= 0 && mod(i,length(eta)) <= 4
                        axis([.15 .4 .4 .62])
                    elseif mod(i,length(eta)) == 5
                        axis([.15 .5 .4 .65])
                    elseif mod(i,length(eta)) == 6
                        axis([0 .75 .4 .68])
                    elseif mod(i,length(eta)) == 0
                        axis([0 1.75 .2 .68])
                    end
                end
            elseif j == 3
                if mod(i,length(eta)) ~= 0 
                    axis([.15 .65 .35 .6])
                elseif mod(i,length(eta)) == 0
                    axis([0 1 0.25 .8])
                end
            elseif j == 4
                if mod(i,length(eta)) ~= 0 && mod(i,length(eta)) <= 5
                    axis([.1539 .3536 .406 .574])
                elseif mod(i,length(eta)) == 6
                    axis([.1 .3536 .406 .6])
                elseif mod(i,length(eta)) == 0
                    axis([0 .4 .35 .65])
                end
            elseif j == 5
                if mod(i,length(eta)) ~= 0 && mod(i,length(eta)) <= 4
                    axis([.2 .45 .44 .57])
                elseif any(mod(i,length(eta)) == [5 6])
                    axis([0 .5 .425 .6])
                elseif mod(i,length(eta)) == 0
                    axis([0 .6 .4 .81])
                end
            end
            
        elseif strcmp(IC_str,'_gauss')
            if j == 1
                if mod(i,length(eta)) <= 5 && mod(i,length(eta)) ~= 0 
                    axis([.1 .35 .35 .6])
                elseif any(mod(i,length(eta)) == [6 7])
                    axis([0 0.5 .35 .8])
                else
                    axis([0 1 .25 1])
                end
%             elseif j == 2
%                 axis(2.^[-2 -1.5 -1.35 -1.26])
            elseif j ==2
                if mod(i,length(eta)) <= 5 && mod(i,length(eta)) ~= 0 
                    axis([.25 .4 .38 .45])
                elseif mod(i,length(eta))>=6
                    axis([.05 .6 .35 .55])
                elseif mod(i,length(eta))==0
                    axis([0 1 .25 .65])
                end
            elseif j == 3
                if mod(i,length(eta))>=6
                    axis([0 .4 .36 .6])
                elseif mod(i,length(eta)) == 0
                    axis([0 .7 .3 .7])
                else
                    axis([0.2 .31 .39 .44])
                end
            elseif j ==4
                if mod(i,length(eta)) == 0
                    axis([0 1 .25 .6])
                elseif mod(i,length(eta)) >= 6
                    axis([0 .4 .35 .65])
                else
                    axis([0.15 0.35 .35 .5])
                end                    
            end
        end
        
        
    end
    
    caxis([1/xnsize(end) 1/xnsize(1)])
    c=colorbar('units','normalized','position',[.93 .45 .01 .15]);
    xns = length(xnsize);
    set(c,'colormap',[1-(1/(xns):(1/(xns)):1)' zeros(xns,1) (1/(xns):1/(xns):1)'])
    set(c,'ticks',[]);
    ax = axes('Position',[0 0 1 1],'visible','off');
    axes(ax)

    a=text(.9325,.62,'$h$','units','normalized','interpreter','latex');
    text(.94,.59,'$(10\cdot2^0)^{-1}$','units','normalized','interpreter','latex')
    text(.94,.525,'$(10\cdot2^3)^{-1}$','units','normalized','interpreter','latex')
    text(.94,.46,'$(10\cdot2^6)^{-1}$','units','normalized','interpreter','latex')

    exportfig(gcf,['CI_h_plot' IC_str '_' num2str(j) stat_meth '.eps'],'fontsize',1.15,'color','rgb')
    saveas(gcf,['CI_h_plot' IC_str '_' num2str(j) stat_meth '.fig'])

    
end

