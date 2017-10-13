 clear all
 clc


% pressure curvature torsion
pressure   = 14; % in psi

alpha = 85*pi/180; beta = 58*pi/180;

option = 0; % 0 for shape without gravity ,1 for shape with gravity and 2 for grasping around a object

L          = 20e-2;
gravity_on = 1;
n_t_1      = 61;%51 for < 20 % 61 for 20
n_f        = 41;%31 for < 20 % 61 for 20
% r_cyl      =  1.3e-2+.0048;
r_cyl      =  .014;
Ro = 13/64*.0254;
Ri = 3/16*.0254;


load Material_fit_6088.mat
Eb = 7e5; % 50 88  (.6134) for 60 88 (.4386) for 70,88
WpL = .0326; % .0332-50, .0326-60, .0345-70
lambda1 = interp1(Material_fit(:,1),Material_fit(:,2),pressure);
delta = L.*Ro.^(-1).*cos(alpha).^3.*cos(beta).^3.*csc(alpha+(-1).*beta) ...
  .^2.*csc(alpha+beta).*(((-1).*lambda1.^2+sec(beta).^2).*tan(alpha) ...
  .^2+sec(alpha).^4.*sec(beta).^4.*(cos(alpha).^6.*((-1)+ ...
  lambda1.^2.*cos(alpha).^2).*cos(beta).^6.*((-1)+lambda1.^2.*cos( ...
  beta).^2).*(tan(alpha)+(-1).*tan(beta)).^4).^(1/2)+(2.*lambda1.^2+ ...
  (-1).*sec(alpha).^2+(-1).*sec(beta).^2).*tan(alpha).*tan(beta)+(( ...
  -1).*lambda1.^2+sec(alpha).^2).*tan(beta).^2);

curvature = (lambda1-1)/(2*Ro);
torsion = -delta/L;
init_g_value = 0; % intial guess value for force

tic
                     
[a,b,initial_shape,d,e]  = getShape(curvature,... % e only for cons3 function
                       torsion,Eb,WpL,...
                       L,...
                       gravity_on,...
                       n_t_1,...
                       n_f,...
                       r_cyl,init_g_value,option);
                   
toc

if option == 0 
    plot3(initial_shape(:,1),initial_shape(:,2),initial_shape(:,3),'r')
else
    plot3(a(:,1),a(:,2),a(:,3),'r')
end
grid on
 axis equal
 xlabel 'X'
 ylabel 'Y'
 zlabel 'Z'
 set(gca,'FontWeight','bold')



