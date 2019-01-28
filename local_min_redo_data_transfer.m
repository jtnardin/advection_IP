
%load new data
load('advection_rates_front_IC_2_23.mat')

%rewrite
J = J_ols;
q = q_ols;

J_a = J_autoreg;
q_a = q_autoreg;

phi1a = phi1;
phi2a = phi2;

%load in old data
load('advection_rates_autoreg_front_IC_old_2_23.mat')


A = redo_trials;


for l = 1:size(A,1)

    t = num2cell(A(l,:));
    [i,j,k] = deal(t{:});
    
    
    %if the cost function has improved, then rewrite
    if J(i,j,k) < J_ols(i,j,k)
        
        J_ols(i,j,k) = J(i,j,k);
        q_ols{i,j,k} = q{i,j,k};
        
        if k <= 2
            J_autoreg(i,j,k) = J_a(i,j,k);


            q_autoreg{i,j,k} = q_a{i,j,k};

            phi1{i,j,k} = phi1a{i,j,k};
            phi2{i,j,k} = phi2a{i,j,k};
        end
    else
        disp(['case ' num2str(i) ' ' num2str(j) ' ' num2str() ' not rewritten'])
    end
end


save(['advection_rates_autoreg' IC_str '_IC.mat'],'q_ols','J_ols',...
    'q_autoreg','J_autoreg','phi1','phi2')