%% Rcylinder/Ractuator = Radius ratio

clf
c = @cmu.colors;

alpha = 50*pi/180;
lambda1 = [1.15:.01:1.29];
Radius_Ratio = ((-1)+4.*((-1)+lambda1).*(1+lambda1).^(-1).*(4.*((-1)+lambda1) ...
  .^2.*(1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1));

plot(lambda1,Radius_Ratio,'--r+','Linewidth',2)

hold on

alpha = 55*pi/180;
lambda1 = [1.15:.01:1.38];
Radius_Ratio = ((-1)+4.*((-1)+lambda1).*(1+lambda1).^(-1).*(4.*((-1)+lambda1) ...
  .^2.*(1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1));

plot(lambda1,Radius_Ratio,'--ro','Linewidth',2)

hold on

alpha = 60*pi/180;
lambda1 = [1.15:.01:1.4];
Radius_Ratio = ((-1)+4.*((-1)+lambda1).*(1+lambda1).^(-1).*(4.*((-1)+lambda1) ...
  .^2.*(1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1));

plot(lambda1,Radius_Ratio,'--kp','Linewidth',2)

hold on
alpha = 65*pi/180;
lambda1 = [1.15:.01:1.4];
Radius_Ratio = ((-1)+4.*((-1)+lambda1).*(1+lambda1).^(-1).*(4.*((-1)+lambda1) ...
  .^2.*(1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1));

plot(lambda1,Radius_Ratio,'--kh','Linewidth',2)

hold on

alpha = 70*pi/180;
lambda1 = [1.15:.01:1.4];
Radius_Ratio = ((-1)+4.*((-1)+lambda1).*(1+lambda1).^(-1).*(4.*((-1)+lambda1) ...
  .^2.*(1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1));

plot(lambda1,Radius_Ratio,'--bs','Linewidth',2)

alpha = 75*pi/180;
lambda1 = [1.15:.01:1.4];
Radius_Ratio = ((-1)+4.*((-1)+lambda1).*(1+lambda1).^(-1).*(4.*((-1)+lambda1) ...
  .^2.*(1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1));

plot(lambda1,Radius_Ratio,'--bd','Linewidth',2)

hold on

set(gca,'linewidth',2,'FontSize',12)

xlabel ('\lambda_1','FontSize',12,'FontWeight','bold')

y  = ylabel('$\frac{r_{wp}}{r_{act}}$','Interpreter','latex','FontSize',22,'FontWeight','bold','Rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.45, 0]);
legend('50 deg','55 deg','60 deg', '65 deg','70 deg','75 deg')
% title('Variation of ratio of radius of cylinder to radius of actuator ({R_{cyl}}/{r_{act}})','FontSize',18)
grid on

set(gcf, 'Position', [100, 100, 600, 600])

h1 = line([1.15 1.3443],[2.31 2.31],'Color','green','LineStyle','-','linewidth',2);
h2 = line([1.15 1.1696],[6.30 6.30],'Color','green','LineStyle','-','linewidth',2);
h3 = line([1.1696 1.1696],[0 6.30],'Color','green','LineStyle','-','linewidth',2);
h4 = line([1.3443 1.3443],[0 2.31],'Color','green','LineStyle','-','linewidth',2);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%% ratio1 is L_workpiece/r_act;;


close all
clf
% 50 degree runs only till 1.29 and 55 till 1.38 (lambda value)
lambda1_50 = [1.15:.01:1.29];
lambda1_55 = [1.15:.01:1.38];

angles = [50:5:75];

alpha = angles(1)*pi/180;
% 50 deg
lambda1 = lambda1_50;
ratio1_50 = (-0.15708E2).*(4.*((-1)+lambda1).^2.*(1+lambda1).^(-2)+(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha)).^2).^(-1).*(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha));
% 55 deg
alpha = angles(2)*pi/180;
lambda1 = lambda1_55;
ratio1_55 = (-0.15708E2).*(4.*((-1)+lambda1).^2.*(1+lambda1).^(-2)+(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha)).^2).^(-1).*(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha));

lambda1 = [1.15:.02:1.39 1.40];

alpha = angles(3:end)'*pi/180;

ratio1_all = (-0.15708E2).*(4.*((-1)+lambda1).^2.*(1+lambda1).^(-2)+(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha)).^2).^(-1).*(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha));


plot(lambda1_50, ratio1_50,'--r+',lambda1_55, ratio1_55,'--ro',...
    lambda1,ratio1_all(1,:),'--kp',lambda1,ratio1_all(2,:),'--kh',...
    lambda1,ratio1_all(3,:),'--bs',...
    lambda1,ratio1_all(4,:),'--bd','Linewidth',2,'markersize',6)

set(gca,'linewidth',2,'FontSize',12)
xlabel('\lambda_1','fontweight','bold','fontsize',12)
y  = ylabel('$\frac{l_{wp}}{r_{act}}$','Interpreter','latex','FontSize',20,'FontWeight','bold','Rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.45, 0]);
grid on

ylim([18 60])
set(gcf, 'Position', [100, 100, 600, 600])
legend('50 deg', '55 deg', '60 deg','65 deg', '70 deg', '75 deg')

%%  ratio 2 is L_act_min/r_act;  NOTE: r_helix = r_cyl+r_act

close all
clf
% 50 degree runs only till 1.29 and 55 till 1.38 (lambda value)
lambda1_50 = [1.15:.01:1.29];
lambda1_55 = [1.15:.01:1.38];

angles = [50:5:75];

alpha = angles(1)*pi/180;
% 50 deg
lambda1 = lambda1_50;
ratio2_50 = 0.314159E2.*((-1)+lambda1).*(1+lambda1).^(-1).*(1+(1/4).*((-1)+ ...
  lambda1).^(-2).*(1+lambda1).^2.*sec(alpha).^2.*(2.^(1/2).*((-1).* ...
  cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2) ...
  .*sec(alpha)+(-2).*sin(alpha)).^2).^(1/2).*(4.*((-1)+lambda1).^2.* ...
  (1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1);
% 55 deg
alpha = angles(2)*pi/180;
lambda1 = lambda1_55;
ratio2_55 = 0.314159E2.*((-1)+lambda1).*(1+lambda1).^(-1).*(1+(1/4).*((-1)+ ...
  lambda1).^(-2).*(1+lambda1).^2.*sec(alpha).^2.*(2.^(1/2).*((-1).* ...
  cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2) ...
  .*sec(alpha)+(-2).*sin(alpha)).^2).^(1/2).*(4.*((-1)+lambda1).^2.* ...
  (1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1);

lambda1 = [1.15:.02:1.39 1.40];

alpha = angles(3:end)'*pi/180;

ratio2_all = 0.314159E2.*((-1)+lambda1).*(1+lambda1).^(-1).*(1+(1/4).*((-1)+ ...
  lambda1).^(-2).*(1+lambda1).^2.*sec(alpha).^2.*(2.^(1/2).*((-1).* ...
  cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2) ...
  .*sec(alpha)+(-2).*sin(alpha)).^2).^(1/2).*(4.*((-1)+lambda1).^2.* ...
  (1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1);


plot(lambda1_50, ratio2_50,'--r+',lambda1_55, ratio2_55,'--ro',...
    lambda1,ratio2_all(1,:),'--kp',lambda1,ratio2_all(2,:),'--kh',...
    lambda1,ratio2_all(3,:),'--bs',...
    lambda1,ratio2_all(4,:),'--bd','Linewidth',2,'markersize',6)

set(gca,'linewidth',2,'FontSize',12)
xlabel('\lambda_1','fontweight','bold','fontsize',12)
y  = ylabel('$\frac{l_{act}}{r_{act}}$','Interpreter','latex','FontSize',20,'FontWeight','bold','Rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.45, 0]);
grid on

ylim([18 100])
set(gcf, 'Position', [100, 100, 600, 600])
legend('50 deg', '55 deg', '60 deg','65 deg', '70 deg', '75 deg')
h1 = line([1.1696 1.1696],[0 75.9],'Color','green','LineStyle','-','linewidth',2);
h2 = line([1.15 1.1696],[75.9 75.9],'Color','green','LineStyle','-','linewidth',2);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


%% Plots for load versus pressure 

Pressure_40 = [12:2:20]*6.89476;
load_40 = [ 275 475 975 1075 1275];

Pressure_30 = [15:1:22]*6.89476;
load_30 = [40 200 400 700 800 900 950 1050];

plot(Pressure_40,load_40,'r--',Pressure_30,load_30,'b--','Linewidth',2)
grid on
set(gca,'linewidth',2,'FontSize',12)
legend('20 mm', '15 mm')
xlabel('Pressure (in kPa)','fontweight','bold','fontsize',12)
ylabel('Load (in gms)','FontSize',12,'FontWeight','bold');
set(gcf, 'Position', [100, 100, 600, 600])

