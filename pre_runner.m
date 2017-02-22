
%% This file used to get the curvature of cylinder which can be grasped at a given pressure and length
clear all
% clc

% for 10 psi 6.7034 13.5069
% 12PSI      8.9947 17.1582
% for 14 psi 10.7905 20.8017   for 16 psi 13.0044 25.4410
% for 18 psi 15.3885 29.7957   for 20 psi 18.1132 34.9487
% for 22 psi 23.0827 43.0181
% the above are different test cases it is in following order
% pressure curvature torsion

L = 60e-2;
for i = 1:length(L)
pressure   = 16; % in psi
curvature  = 13.0044;
torsion    = 25.4410;
% pause
  
gravity_on = 1;
n_t_1      = 101;%51 for < 20 % 61 for 20
n_f        = 101;%31 for < 20 % 61 for 20
% r_cyl      =  1.3e-2+.0048;
r_cyl      =  1.3e-2;

% tic
                     
a  = AFTER_GRAVITY(curvature,... % e only for cons3 function
                       torsion,...
                       L,...
                       gravity_on,...
                       n_t_1,...
                       n_f,...
                       r_cyl);
% toc

TangVec = angle_detect(a);

CylCurv(i) = a(end,13)/(sind(TangVec(end)))^2;

MinActGrasp(i) = 5/2*pi*CylCurv^(-1)/sind(mean(TangVec));

[L(i) CylCurv(i) MinActGrasp(i)]
% pause
end