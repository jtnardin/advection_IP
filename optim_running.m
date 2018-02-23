%optim_running written 2-3-18 by JTN to run optimization routine
%several times over

q_ols = cell(7,8,4);
J_ols = zeros(7,8,4);

q_autoreg = cell(7,8,2);
J_autoreg = zeros(7,8,2);

phi1 = cell(7,8,2);
phi2 = cell(7,8,2);

IC_str = '_front';

A = redo_trials;


tic
% for i = 1:7
%     for j = [2 3 6 7]
%         for k = 1:4
%             [i,j,k]

for l = 1:size(A,1)

    t = num2cell(A(l,:));
    [i,j,k] = deal(t{:})
    
    if k == 1
        num_meth = 'upwind';
    elseif k == 2
        num_meth = 'laxfried';
    elseif k == 3
        num_meth = 'laxwend';
    elseif k == 4
        num_meth = 'beamwarm';
    end


    if k > 2
        [q_ols{i,j,k},J_ols(i,j,k)]= art_advec_fitting_f(i,j,num_meth,IC_str);

    else   %autoregression for upwind, laxfried
        [q_ols{i,j,k},J_ols(i,j,k),q_autoreg{i,j,k},J_autoreg(i,j,k)...
        ,phi1{i,j,k},phi2{i,j,k}]...
        = autoreg_art_advec_fitting_f(i,j,num_meth,IC_str);
    end
            
%         end
%     end
% end
end
toc

% save(['advection_rates' IC_str '_IC_redo.mat'],'q_ols','J_ols');%,...


save(['advection_rates_autoreg' IC_str '_IC.mat'],'q_ols','J_ols',...
    'q_autoreg','J_autoreg','phi1','phi2')