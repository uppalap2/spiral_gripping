clear all
% clc

% for 10 psi 6.7034 13.5069
% 12PSI      8.9947 17.1582
% for 14 psi 10.7905 20.8017   for 16 psi 13.0044 25.4410
% for 18 psi 15.3885 29.7957   for 20 psi 18.1132 34.9487
% for 22 psi 23.0827 43.0181
% the above are different test cases it is in following order
% pressure curvature torsion
pressure   = 12; % in psi
curvature  = 8.9947;
torsion    = 17.1582;
L          = 57e-2;
gravity_on = 1;
n_t_1      = 61;%51 for < 20 % 61 for 20
n_f        = 61;%31 for < 20 % 61 for 20
% r_cyl      =  1.3e-2+.0048;
r_cyl      =  1.3e-2;

tic
                     
[a,b,c,d,e]  = FINAL_TILL_CURL(curvature,... % e only for cons3 function
                       torsion,...
                       L,...
                       gravity_on,...
                       n_t_1,...
                       n_f,...
                       r_cyl);
toc
%         load 22psi_guess.mat
%         load 22psi_force_guess.mat
%         
% [a,b,c,d,e]  = opti_SE_4_main(curvature,... % e only for cons3 function
%                        torsion,...
%                        L,...
%                        gravity_on,...
%                        n_t_1,...
%                        n_f,...
%                        r_cyl, a, b);
toc

