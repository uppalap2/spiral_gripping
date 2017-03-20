function initial_shape = getShapeGravity(pressure, alpha , beta,fitmatFile,L)
% fit file is a string with .mat extension
% alpha and beta in radians, L in meters
% pressure curvature torsion
% pressure   = 18; % in psi

% alpha = 85*pi/180; beta = 58*pi/180;

option = 1; % 1 for shape with gravity and 2 for grasping around a object

% L          = 45e-2;
gravity_on = 1;
n_t_1      = 61;%51 for < 20 % 61 for 20
n_f        = 61;%31 for < 20 % 61 for 20
% r_cyl      =  1.3e-2+.0048;
r_cyl      =  1.3e-2;
Ro = 13/64*.0254;
Ri = 3/16*.0254;


% load Material_fit_5088.mat
load(fitmatFile)
Eb = 6*.4386e6; %  (.6134) for 60 88 (.4386) for 70,88 (.6029) 50,88

lambda1 = interp1(Material_fit(:,1),Material_fit(:,2),pressure);
delta = L.*Ro.^(-1).*cos(alpha).^3.*cos(beta).^3.*csc(alpha+(-1).*beta) ...
  .^2.*csc(alpha+beta).*(((-1).*lambda1.^2+sec(beta).^2).*tan(alpha) ...
  .^2+sec(alpha).^4.*sec(beta).^4.*(cos(alpha).^6.*((-1)+ ...
  lambda1.^2.*cos(alpha).^2).*cos(beta).^6.*((-1)+lambda1.^2.*cos( ...
  beta).^2).*(tan(alpha)+(-1).*tan(beta)).^4).^(1/2)+(2.*lambda1.^2+ ...
  (-1).*sec(alpha).^2+(-1).*sec(beta).^2).*tan(alpha).*tan(beta)+(( ...
  -1).*lambda1.^2+sec(alpha).^2).*tan(beta).^2);

curvature = (lambda1-1)/(2*Ro);
torsion = delta/L;


tic
                     
[~,~,initial_shape,~,~]  = getShape(curvature,... % e only for cons3 function
                       torsion,Eb,...
                       L,...
                       gravity_on,...
                       n_t_1,...
                       n_f,...
                       r_cyl,option);
                   
toc

% plot3(initial_shape(:,1),initial_shape(:,2),initial_shape(:,3),'r')
% grid on
%  axis equal
%  xlabel 'X'
%  ylabel 'Y'
%  zlabel 'Z'
%  set(gca,'FontWeight','bold')



