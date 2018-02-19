%%%%%% order_estimate.m written 2-5-18 by JTN
%%% mean to find some order estimates for J(h), q(h)

clear all; clc


IC_str = '_front';


%load best-fit params, data, and initial condition
load(['advection_rates' IC_str '_IC.mat'])
load(['advection_art_data' IC_str '.mat'])
load('rel_range_sims.mat')


xnsize = [21,41,81,161,321,641,2*640+1];

num_meth_cell = cell(4,1);
num_meth_cell{1} = 'upwind';
num_meth_cell{2} = 'Lax-Friedrich';
num_meth_cell{3} = 'Lax-Wendroff';
num_meth_cell{4} = 'Beam warming';

%%%%%% J computation

%estimate order using best-fit line
order_table = zeros(4,4);
poly_table = cell(4,4);

for i = 1:4
    for j = 1:4
        
%         %%%%only include points from a decreasing subsequence
%         rel_data = false(length(xnsize),1);
%         rel_data(1) = 1;
%         %smallest cost value to date
%         col_min = J(1,i,j);
%         
%         for k = 2:length(xnsize)
%             if J(k,i,j) < col_min
%                rel_data(k)=1;
%                col_min = J(k,i,j);
%             end
%         end
%         
%        
%         %if J not really decreasing
%         if sum(rel_data) < 3
%             poly_table{i,j} = NaN;
%             order_table(i,j) = NaN;
%         else %consider if decreased at least twice
%             poly_table{i,j} = polyfit(log(1./xnsize(rel_data))',...
%                 log(squeeze(J(rel_data,i,j))),1);
%             order_table(i,j) = poly_table{i,j}(1);
%         end 


            
            if strcmp(IC_str,'_gauss')
                rel_range = rel_range_gauss{j};
            elseif strcmp(IC_str,'_front')
                rel_range = rel_range_front{j};
            end


            %include all points
            poly_table{i,j} = polyfit(log(1./xnsize(rel_range))',log(squeeze(J(rel_range,i,j))),1);
            order_table(i,j) = poly_table{i,j}(1);



    end
end


%estimate order using 2-point methods
R=abs(J(1:end-1,:,:)./J(2:end,:,:));

phat = log2(R);


%%%%%%now do the same for q

%get norm of ||q-q0||_2
q_norm = zeros(7,4,4);

for i = 1:7
    for j = 1:4
        for k = 1:4
            
            if J(i,j,k)>0
                q_norm(i,j,k) = norm(q{i,j,k}-q0);
            end
            
        end
    end
end

%%%%% write to latex table.

write_latex_table(['J_1' IC_str '.tex'],order_table')
write_latex_table(['J_2' IC_str '.tex'],squeeze(mean(phat)))
% 

%estimate order using best-fit line
q_order_table = zeros(4,4);
q_poly_table = cell(4,4);

for i = 1:4
    for j = 1:4
        
       if strcmp(IC_str,'_gauss')
            rel_range = rel_range_gauss{j};
        elseif strcmp(IC_str,'_front')
            rel_range = rel_range_front{j};
        end

        
        q_poly_table{i,j} = polyfit(log(1./xnsize(rel_range))',log(squeeze(q_norm(rel_range,i,j))),1);
        q_order_table(i,j) = q_poly_table{i,j}(1);
             
    end
end


%estimate order using 2-point methods
qR=abs(q_norm(1:end-1,:,:)./q_norm(2:end,:,:));

qphat = log2(qR);


%%%%% write to latex table.

write_latex_table(['q_1' IC_str '.tex'],q_order_table);
write_latex_table(['q_2' IC_str '.tex'],squeeze(mean(qphat)));
