a =[


161.1683 92.2013 540.8283
106.4284 44.3115 172.7018
 
];

a =(a-a(1,:) -[0 .0048*1000 0])/1000;
theta = 0 ;
Rot_z = [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1]
a = a*Rot_z
plot3(a(:,1),a(:,2),-a(:,3),'ko','markersize',10,'Markerfacecolor', 'green')
hold on
axis equal
grid on
xlabel( 'X','fontweight','bold','fontsize',20)
ylabel ('Y','fontweight','bold','fontsize',20)
zlabel ('Z','fontweight','bold','fontsize',20)
set(gca,'linewidth',1.5,'FontSize',16)
[x,y,z]=tubeplot([initial_shape(:,1) initial_shape(:,2) initial_shape(:,3)]',.0048)