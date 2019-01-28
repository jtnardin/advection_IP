clear all; clc


load('CI_front_OLS_all.mat')
N = [11,31,51];
hn = [21,41,81,161,321,641,2*640+1];

order1 = zeros(7,1);
order2 = zeros(7,1);

for i = 1:7
    width = zeros(3,2);
    count = 1;
    for j = 2:7:21
        C = CI{i,j,5};
        width(count,:)=diff(C,1,2);
        count = count+1;
    end

    p1 = polyfit(log(1./N)',log(width(:,1)),1);
    order1(i) = p1(1);

    p2 = polyfit(log(1./N)',log(width(:,2)),1);
    order2(i) = p2(1);
end

order1
order2

% order1 = zeros(21,1);
% order2 = zeros(21,1);
% 
% for i = 1:21
% 
%     width = zeros(7,2);
% 
%     count = 1;
%     for j = 1:7
%         C = CI{j,i,4};
%         width(count,:)=diff(C,1,2);
%         count = count+1;
%     end
% 
%     p1 = polyfit(log(1./(hn-1))',log(width(:,1)),1);
%     order1(i) = p1(1);
% 
%     p2 = polyfit(log(1./(hn-1))',log(width(:,2)),1);
%     order2(i) = p2(1);
%     
% end
%     
% order1
% order2