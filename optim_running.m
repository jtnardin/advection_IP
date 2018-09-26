%optim_running written 2-3-18 by JTN to run optimization routine
%several times over
clear all; clc
%q_ols = cell(7,8,4);
%J_ols = zeros(7,8,4);

IC_str = '_front';

model_str = 'root';

q_ols = cell(7,21,5);
J_ols = zeros(7,21,5);
num_it = zeros(7,21,5);
num_it_auto = zeros(7,21);

data_str = 'advection_art_data_front_all_rev.mat';

%load(['advection_rates' IC_str '_IC_' model_str '.mat'])

% q_autoreg = cell(7,8,2);
% J_autoreg = zeros(7,8,2);
% 
% phi1 = cell(7,8,2);
% phi2 = cell(7,8,2);

% A = redo_trials;


for i = 1:4%7
    for j = [8:14]%1:21
        for k = [1 3:5]
            [i,j,k]

            if isempty(q_ols{i,j,k})

                % for l = 1:size(A,1)
                % 
                %     t = num2cell(A(l,:));
                %     [i,j,k] = deal(t{:})

                    if k == 1
                        num_meth = 'upwind';
                    elseif k == 2
                        num_meth = 'laxfried';
                    elseif k == 3
                        num_meth = 'laxwend';
                    elseif k == 4
                        num_meth = 'beamwarm';
                    elseif k ==5
                        num_meth = 'upwindfl';
                    end

                    tic
                     if k > 2
                        [q_ols{i,j,k},J_ols(i,j,k),num_it(i,j,k)]=art_advec_fitting_f(i,j,num_meth,IC_str,model_str,data_str);
                 
                     else   %autoregression for upwind, laxfried
                         [q_ols{i,j,k},J_ols(i,j,k),q_autoreg{i,j,k},J_autoreg(i,j,k)...
                         ,phi1{i,j,k},phi2{i,j,k},num_it(i,j,k),num_it_auto(i,j)]...
                         = autoreg_art_advec_fitting_f(i,j,num_meth,IC_str,model_str,data_str);
                     end
                     toc
            
            end
        end
    end
end
% end


% save(['advection_rates' IC_str '_IC_upwindfl_3_17_6.mat'],'q_ols','J_ols');%,...


% save(['advection_rates_autoreg' IC_str '_IC_3_1.mat'],'q_autoreg',...
%     J_autoreg','phi1','phi2')
