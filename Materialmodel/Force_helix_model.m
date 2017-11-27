%% Basicallly we estimate the forces on the actuator after it has reached its critical pressure.
%At critical pressure it conforms to the cylinder and has its shape fixed,
%when increase pressure it is applying more force. We are calculating this
%force. ASSUMPTION: Circular helix assumption and no change of the
%actuator orientation

lambda1 = [1.12 1.14];
L = 50e-2;
l = L;
r =  0.013/2;
alpha = 55*pi/180;
kappa1 = ((-1)+lambda1).*(1+lambda1).^(-1).*r.^(-1);
tau1 = (-1/2).*l.^(-1).*r.^(-2).*(2.^(1/2).*l.*((-1).*r.^2.*cos(alpha) ...
  .^2.*((-2)+lambda1.^2+lambda1.^2.*cos(2.*alpha))).^(1/2).*sec( ...
  alpha).^2+(-2).*l.*r.*tan(alpha));



u = [kappa1(1) 0 tau1(1)];
u0 = [kappa1(2:end)' zeros(length(lambda1)-1,1) tau1(2:end)'];
% u = [curv_critical 0 tor_critical];
% u0 = [new_curvature 0 new_torsion];
Eb = 7e5;
WpL = .0326;

[~,~,initial_shape,~,~,C] = getShape(u(1),...
                                  u(3),Eb,WpL,...
                                  L,...
                                  0,...
                                  100,...
                                  10,...
                                  0,10,0);
                           
% function [sxint,...
%           F,...
%           initial_shape,...
%           cyl_center,...
%           const,C] = getShape(curvature,...
%                                   torsion,Eb,WpL,...
%                                   L,...
%                                   gravity_on,...
%                                   n_t,...
%                                   n_f,...
%                                   r_cyl,init_g_value,option)                           
v = [0 0 1]';
R_g = zeros(3,3,length(initial_shape));
F = zeros(3,length(initial_shape));
for i = 1:length(initial_shape)
R_g(:,:,i) = reshape(initial_shape(i,4:12),[3 3])';

F(:,i) = (skew(v)*R_g(:,:,i)')\ skew(u')*diag(C)*(u0-u)';
end


% F = inv(skew(v)*R_g')\ skew(u')*diag(C)*(u0-u)';


function ahat = skew(a)
ahat = [0 -a(3) a(2);a(3) 0 -a(1); -a(2) a(1) 0];
end
                              