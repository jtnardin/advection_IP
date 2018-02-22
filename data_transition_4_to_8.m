load('advection_rates_gauss_IC.mat')
load('advection_rates_autoreg_gauss_IC_2_19_first.mat')

% J_autoreg_old = J_autoreg
% J_ols_old = J_ols;
% phi1_old = phi1;
% phi2_old = phi2;
% q_autoreg_old = q_autoreg;
% q_ols_old = q_ols;


J_ols(:,1:2,:) = J(:,1:2,:);
J_ols(:,5:6,:) = J(:,3:4,:);

J_ols(:,1:2,:) = J(:,1:2,:);
J_ols(:,5:6,:) = J(:,3:4,:);


for i = 1:7
    for j = 1:2
        for k = 1:4
            q_ols{i,j,k} = q{i,j,k};
        end
    end
end


for i = 1:7
    for j = 3:4
        for k = 1:4
            q_ols{i,j+2,k} = q{i,j,k};
        end
    end
end

save('advection_rates_autoreg_gauss_IC.mat','J_ols','q_ols','J_autoreg','q_autoreg','phi1','phi2')