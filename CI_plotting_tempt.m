
figure
subplot(1,2,1)
hold on

for i = 1:7
    
    C = CI{i,15,1};
    
    rectangle('position',[C(1,1) C(2,1) C(1,2)-C(1,1) C(2,2)-C(2,1)],...
                    'edgecolor',repmat(1-(i)/7,1,3))
                
                               
end

plot(q0(1),q0(2),'r*','markersize',15)

axis([.05 .4 .35 .6])

subplot(1,2,2)
hold on

for i = 1:7
    
    C = CI2{i,15,1};
    
    rectangle('position',[C(1,1) C(2,1) C(1,2)-C(1,1) C(2,2)-C(2,1)],...
                    'edgecolor',repmat(1-(i)/7,1,3))
                
end

plot(q0(1),q0(2),'r*','markersize',15)

axis([.05 .4 .35 .6])