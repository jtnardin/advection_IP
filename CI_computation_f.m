%Used to simulate sensitivity equations, and then compute confidence intervals
%written 2-22-18 by JTN

% clear all;

%%%%specify which initial condition, method,
%%%%and data set we're compating

function CI_computation_f(m)

    IC_str = '_front';

    stat_meth = 'OLS';

    % load(['CI' IC_str '_' stat_meth '_all.mat'])

    if strcmp(IC_str,'_front')

        load(['advection_rates_autoreg' IC_str '_IC_all.mat'])
        load(['CI' IC_str '_' stat_meth '_newton_3_17_upwindfl.mat'])

    elseif strcmp(IC_str,'_gauss')
        load(['advection_rates' IC_str '_IC_all.mat'])
    end

    load(['advection_art_data' IC_str '_all.mat'])
    phi = IC_spec(IC_str(2:end));



     %grid sizes
    xnsize = [21,41,81,161,321,641,2*640+1];
    lambda = 1/20;

    for i = 1:length(xd)
        xndata(i) = length(xd{i});
    end




    for mn = 1:numel(data)
        xdi = ceil(mn/length(eta));
        sigmaj = mod(mn,length(eta));

        if sigmaj == 0
            sigmaj = length(eta);
        end

        xnstr(mn) = xndata(xdi);
        eta_str(mn) = eta(sigmaj);
    end

    num_meth_cell = cell(8,1);
    num_meth_cell{1} = 'upwindsens';
    num_meth_cell{2} = 'laxfriedsens';
    num_meth_cell{3} = 'laxwendsens';
    num_meth_cell{4} = 'beamwarmsens';
    num_meth_cell{5} = 'upwindflsens';



    for xni = 6:7
%         for m = 1:numel(data)
            for num_meth = 5

                if isempty(CI{xni,m,num_meth})

                %load best-fit params, data, and initial condition

                    [xni,m,num_meth]

                    if strcmp(stat_meth,'OLS')
                        q = q_ols{xni,m,num_meth};
                    elseif strcmp(stat_meth,'autor')
                        q = q_autoreg{xni,m,num_meth};
                    else
                        error('incorrect statistical model')
                    end


                    xdi = ceil(m/length(eta));
                    sigmaj = mod(m,length(eta));

                    if sigmaj == 0
                        sigmaj = length(eta);
                    end


                    %select exp. data
                    cell_data = data{xdi,sigmaj};
                    tdata = td;
                    xdata = xd{xdi};


                    %space
                    x = linspace(0,3,3*xnsize(xni)-2);
                    dx = x(2) - x(1);
                    xn = 3*xnsize(xni)-2;

                    %time
                    dt = lambda*dx;
                    tfin = 10;
                    t = 0:dt:tfin;
                    tn = length(t);

                    %points for computaiton
                    [x_int,xbd_0,xbd_1] = int_bd_def(xn);

                    %initial condition
                    IC = [phi(x) zeros(1,2*xn)];


                    %load matrices for computation
                    [A,Abd] = aMatrixupwind(xn,num_meth_cell{num_meth});


                    tic

                        if strcmp(stat_meth,'OLS')

                            %get model sim
                            [model_sim,s1model,s2model,eta_hat,res] = ...
                                sensitivity_model_sim(cell_data,q,dx,xn,x_int,xbd_0,...
                                xbd_1,dt,tn,IC,A,Abd,x,xdata,num_meth_cell{num_meth},...
                                t,tdata);
                        elseif strcmp(stat_meth,'autor')
                            [model_sim,s1model,s2model,eta_hat,res] = ...
                                sensitivity_model_auto_sim(cell_data,q,dx,xn,x_int,xbd_0,...
                                xbd_1,dt,tn,IC,A,Abd,x,xdata,num_meth_cell{num_meth},...
                                t,tdata,phi1{xni,m,num_meth},phi1{xni,m,num_meth});
                        end

                    toc



                    %compute standard error
                    F = [s1model(:) s2model(:)];
                    Sigma = eta_hat*inv(F'*F);

                    %now construct the 95% confidence interval
                    CI{xni,m,num_meth} = zeros(2,2);

                    T_se = tinv(1-.05/2,numel(cell_data)-length(q));

                    for i = 1:2
                        CI{xni,m,num_meth}(i,:) = [q(i)-T_se*sqrt(Sigma(i,i)) q(i)+T_se*sqrt(Sigma(i,i))];
                    end

                    save(['/scratch/summit/jona8898/chem_fitting/CI' IC_str '_' stat_meth '_newton_3_17_upwindfl_' num2str(xni) '_' num2str(m) '.mat'],'CI')

                end    

            end
%         end
    end

    

end


