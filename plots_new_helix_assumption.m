%% Rcylinder/Ractuator = Radius ratio

clf


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

plot(lambda1,Radius_Ratio,'--gp','Linewidth',2)

hold on
alpha = 65*pi/180;
lambda1 = [1.15:.01:1.4];
Radius_Ratio = ((-1)+4.*((-1)+lambda1).*(1+lambda1).^(-1).*(4.*((-1)+lambda1) ...
  .^2.*(1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1));

plot(lambda1,Radius_Ratio,'--gh','Linewidth',2)

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

set(gca,'linewidth',2,'FontSize',16)

xlabel ('\lambda_1','FontSize',18,'FontWeight','bold')

y  = ylabel('$\frac{R_{cyl}}{r_{act}}$','Interpreter','latex','FontSize',24,'FontWeight','bold','Rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);
legend('50 deg','55 deg','60 deg', '65 deg','70 deg','75 deg')
title('Variation of ratio of radius of cylinder to radius of actuator ({R_{cyl}}/{r_{act}})','FontSize',18)
grid on

%% ratio1 is L_workpiece/r_act;;
% ratio 2 is L_act_min/r_act;  NOTE: r_helix = r_cyl+r_act
figure
close all

clf


lambda1 = [1.15:.01:1.29];

angles = [50:10:70];

alpha = angles(1)*pi/180;

ratio1 = (-0.15708E2).*(4.*((-1)+lambda1).^2.*(1+lambda1).^(-2)+(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha)).^2).^(-1).*(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha));

ratio2 = 0.314159E2.*((-1)+lambda1).*(1+lambda1).^(-1).*(1+(1/4).*((-1)+ ...
  lambda1).^(-2).*(1+lambda1).^2.*sec(alpha).^2.*(2.^(1/2).*((-1).* ...
  cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2) ...
  .*sec(alpha)+(-2).*sin(alpha)).^2).^(1/2).*(4.*((-1)+lambda1).^2.* ...
  (1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1);

yyaxis left
plot(lambda1, ratio1,'--b+','Linewidth',2,'markersize',6)

yyaxis right
plot(lambda1, ratio2,'--r+','Linewidth',2,'markersize',6)

xlabel('\lambda_1','fontweight','bold','fontsize',18)

yyaxis left

yleft = ylabel('$\frac{L_{WP}}{r_{act}}$','Fontsize',30,'fontweight','bold','Interpreter','latex','Rotation',0);
set(yleft, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);
% grid on
yyaxis right
yright = ylabel('$\frac{l_{act}}{r_{act}}$','Fontsize',30,'fontweight','bold','Interpreter','latex','Rotation',0);
set(yright, 'Units', 'Normalized', 'Position', [1.05, 0.5, 0]);



hold on
lambda1 = [1.15:.01:1.4];

alpha = angles(2)*pi/180;

ratio1 = (-0.15708E2).*(4.*((-1)+lambda1).^2.*(1+lambda1).^(-2)+(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha)).^2).^(-1).*(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha));

ratio2 = 0.314159E2.*((-1)+lambda1).*(1+lambda1).^(-1).*(1+(1/4).*((-1)+ ...
  lambda1).^(-2).*(1+lambda1).^2.*sec(alpha).^2.*(2.^(1/2).*((-1).* ...
  cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2) ...
  .*sec(alpha)+(-2).*sin(alpha)).^2).^(1/2).*(4.*((-1)+lambda1).^2.* ...
  (1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1);

yyaxis left
plot(lambda1, ratio1,'--bp','Linewidth',2,'markersize',6)

yyaxis right
plot(lambda1, ratio2,'--rp','Linewidth',2,'markersize',6)

%
alpha = angles(3)*pi/180;

ratio1 = (-0.15708E2).*(4.*((-1)+lambda1).^2.*(1+lambda1).^(-2)+(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha)).^2).^(-1).*(2.^( ...
  1/2).*((-1).*cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.* ...
  alpha))).^(1/2).*sec(alpha).^2+(-2).*tan(alpha));

ratio2 = 0.314159E2.*((-1)+lambda1).*(1+lambda1).^(-1).*(1+(1/4).*((-1)+ ...
  lambda1).^(-2).*(1+lambda1).^2.*sec(alpha).^2.*(2.^(1/2).*((-1).* ...
  cos(alpha).^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2) ...
  .*sec(alpha)+(-2).*sin(alpha)).^2).^(1/2).*(4.*((-1)+lambda1).^2.* ...
  (1+lambda1).^(-2)+(2.^(1/2).*((-1).*cos(alpha).^2.*((-2)+ ...
  lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec(alpha).^2+(-2) ...
  .*tan(alpha)).^2).^(-1);

yyaxis left
plot(lambda1, ratio1,'--bs','Linewidth',2,'markersize',6)

yyaxis right
plot(lambda1, ratio2,'--rs','Linewidth',2,'markersize',6)



set(gca,'linewidth',2,'FontSize',16)
yyaxis left
 ylim([19 90])

yyaxis right
ylim([19 90])

legend('50 deg', '60 deg', '70 deg')
title('Minimum Length of workpiece (L_{WP}) and actuator length(l_{act}) needed ','FontSize',18)

