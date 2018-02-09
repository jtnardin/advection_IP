
m = 3;
num_meth = 1;

phihat = zeros(6,6);

res_mod = cell(6,1);


for xni = 1:6

    res_mod{xni} = zeros(size(res{xni}));
    
    for tj = 2:6

        a = res{xni}(tj,:);

        [g,sigma,sigma_inv] = advection_rate('root',q{xni+1,m,num_meth}(1),...
            q{xni+1,m,num_meth}(2));

        x_past_shock = a(xdata>sigma_inv(tdata(tj),.3));

        ind = find(a==max(x_past_shock));

        phihat(xni,tj) = sum(a(ind:end-1).*a(ind+1:end))/sum(a(ind:end).^2);
        B = sparse([ind ind+1:50 ind+1:50],[ind ind+1:50 ind:49],[sqrt(1-phihat(xni,tj)^2)...
            1*ones(1,50-ind) -phihat(xni,tj)*ones(1,50-ind)],50,50)/sqrt(1-phihat(xni,tj)^2);

        res_mod{xni}(tj,:) = B*a';

    %     figure
    %     subplot(2,1,1)
    %     plot(a,'.')
    %     axis([ind-1 50 -.1 2])
    % 
    %     subplot(2,1,2)
    %     plot(B*a','.')
    %     axis([ind-1 50 -.1 2])

    end

end

