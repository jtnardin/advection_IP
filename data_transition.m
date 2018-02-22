load('advection_rates_autoreg_gauss_IC.mat')
J_ols1 = J_ols;
q_ols1 = q_ols;

load('advection_rates_autoreg_gauss_IC_2_20.mat')
J_ols2 = J_ols;
q_ols2 = q_ols;


J_ols = zeros(7,8,4);
q_ols = cell(7,8,4);

J_ols(:,[1 4 5 8],:) = J_ols1(:,[1 3 5 7],:);

take_vec = [1 3 5 7];
put_vec = [1 4 5 8];

for i = 1:7
    for j = 1:4
        for k = 1:4
            q_ols{i,put_vec(j),k} = q_ols1{i,take_vec(j),k};
        end
    end
end


J_ols(:,[2 3 6 7],:) = J_ols2(:,[2 3 6 7],:);

take_vec = [2 3 6 7];
put_vec = [2 3 6 7];

for i = 1:7
    for j = 1:4
        for k = 1:4
            q_ols{i,put_vec(j),k} = q_ols2{i,take_vec(j),k};
        end
    end
end

% save('advection_rates_front_IC_2_20.mat','J_ols','q_ols')




clear all



load('advection_art_data_gauss.mat')

data1 = data;


load('advection_art_data_gauss2.mat')

data2 = data;

data = cell(2,4);

for i = 1:2
    data{i,1} = data1{i,1};
    data{i,4} = data1{i,3};
    data{i,2} = data2{i,2};
    data{i,3} = data2{i,3};
end

eta = [0 .0001 .001 .01];

save('advection_art_data_gauss3.mat','data','eta','q0','td','xd')

