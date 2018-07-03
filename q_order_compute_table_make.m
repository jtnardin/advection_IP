%q_order_compute_table_make.m written to estimate the order of
%||\hat{q}-q_0||_2, written 6-29-18 by JTN

clear all; clc



xnsize = [21,41,81,161,321,641,2*640+1];
h=1./(xnsize-1);

for k = 1:2

    if k == 1
        num_meth_cell = cell(5,1);
        num_meth_cell{1} = 'upwind';
        num_meth_cell{2} = 'laxwend';
        num_meth_cell{3} = 'beamwarm';
        num_meth_cell{4} = 'upwindfl';


        num_meth_short_cell = cell(5,1);
        num_meth_short_cell{1} = 'Upwind';
        num_meth_short_cell{2} = 'Lax-Wendroff';
        num_meth_short_cell{3} = 'Beam-Warming';
        num_meth_short_cell{4} = 'Upwind w/ flux limiters';

        num_range = 1:4;
        
    elseif k ==2 
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

        num_range = [1 3:5];
    end

    if k == 1
        IC_str = '_gauss';
    elseif k == 2
        IC_str = '_front';
    end

    model_str = 'root';


        clear eta eta_vec
        %load best-fit params, data, and initial condition
        if strcmp(IC_str,'_gauss')
            load(['advection_rates' IC_str '_IC_all_3_26.mat'])
            load(['advection_art_data' IC_str '_all_3_26.mat'])
        elseif strcmp(IC_str,'_front')
            load(['advection_rates_autoreg' IC_str '_IC_all.mat'])
            load(['advection_art_data' IC_str '_all.mat'])

        end


        q_table_order = zeros(4*length(xd),length(eta));
        data_range = 1:numel(data);
        
        row_ind = 3;

        
        
    
            
            
    count = 1;
    for num_meth = num_range
%         count_col = 1;
        for j = data_range
        
            %get data, noise indices
            xdi = ceil(j/length(eta));
            sigmaj = mod(j,length(eta));

            if sigmaj == 0
                sigmaj = length(eta);
            end

        
            q_norm = zeros(7,1);
            for i = 1:7
                q_norm(i) = norm(q_ols{i,j,num_meth}-q0);
            end

            range = 1:7;%J_order_range_define(num_meth_cell{num_meth},IC_str,sigmaj,xdi);
            
            p = polyfit(log(h(range)),log(q_norm(range))',1);
            
            q_table_order(row_ind*(count-1) + xdi,sigmaj) = p(1);
            
%             count_col = count_col + 1;
        end

        
       
        count = count + 1;

    end
    
    write_latex_table(['q_conv_table' IC_str '_.tex'],q_table_order)
    
end

