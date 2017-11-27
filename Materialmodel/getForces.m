function [sxint,F] = getForces(pressure, alpha , beta,fitmatFile,L,Eb,WpL,r_cyl, init_g_value)
% fit file is a string with .mat extension
% alpha and beta in radians, L in meters
% pressure curvature torsion
% pressure   = 18; % in psi

% alpha = 85*pi/180; beta = 58*pi/180;

% option = 1; % 0 for shape w/o gravtyty 
            % 1 for shape with gravity 
            % 2 for grasping around a object

% L          = 45e-2;

option     = 2; % getting forces
gravity_on = 1;

n_t_1      = 81;%51 for < 20 % 61 for 20
n_f        = 61;%31 for < 20 % 61 for 20
Ro         = 13/64*.0254;
Ri         = 3/16*.0254;
% init_g_value = 10;


% load Material_fit_5088.mat
load(fitmatFile)
% Eb = 6*.4386e6; %  (.6134) for 60 88 (.4386) for 70,88 (.6029) 50,88

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
                     
[sxint,...
          F,...
          initial_shape,...
          cyl_center,...
          const]  = getShape(curvature,... % e only for cons3 function
                       torsion,Eb,WpL,...
                       L,...
                       gravity_on,...
                       n_t_1,...
                       n_f,...
                       r_cyl,init_g_value,option);
                   
plot3(sxint(:,1),sxint(:,2),sxint(:,3), 'r')                   
trapz([0:L/(n_f-1):L], F)                   
toc

% plot3(initial_shape(:,1),initial_shape(:,2),initial_shape(:,3),'r')
% grid on
%  axis equal
%  xlabel 'X'
%  ylabel 'Y'
%  zlabel 'Z'
%  set(gca,'FontWeight','bold')



