

lambda1 = [1.000:.01:1.5];

for n = 1:2

switch n
    case 1
%% For 88 60 angles

C1 =  0.551501509158298   ; C2 = 0.000000005245973; % for default
Pvslambda = lambda1.^(-3).*(C2.*(0.113499E9+(-0.925926E0).*lambda1.^12+ ...
  lambda1.^6.*(0.589272E8+(-0.481951E7).*((-1).*((-4)+lambda1.^2).*( ...
  2+lambda1.^2.*((-1)+cos((1/45).*pi)))).^(1/2))+lambda1.^10.*( ...
  0.133686E4+(-0.187604E2).*((-1).*((-4)+lambda1.^2).*(2+ ...
  lambda1.^2.*((-1)+cos((1/45).*pi)))).^(1/2))+lambda1.^8.*(( ...
  -0.55452E6)+0.193475E5.*((-1).*((-4)+lambda1.^2).*(2+lambda1.^2.*( ...
  (-1)+cos((1/45).*pi)))).^(1/2))+lambda1.^2.*((-0.28788E8)+ ...
  0.140044E7.*((-1).*((-4)+lambda1.^2).*(2+lambda1.^2.*((-1)+cos(( ...
  1/45).*pi)))).^(1/2))+lambda1.^4.*((-0.201449E9)+0.174455E9.*((-1) ...
  .*((-4)+lambda1.^2).*(2+lambda1.^2.*((-1)+cos((1/45).*pi)))).^( ...
  1/2))+(-0.147214E9).*((-1).*((-4)+lambda1.^2).*(2+lambda1.^2.*(( ...
  -1)+cos((1/45).*pi)))).^(1/2))+C1.*(0.497884E8+0.833144E2.* ...
  lambda1.^10+lambda1.^6.*(0.286096E8+(-0.13927E7).*((-1).*((-4)+ ...
  lambda1.^2).*(2+lambda1.^2.*((-1)+cos((1/45).*pi)))).^(1/2))+ ...
  lambda1.^8.*((-0.103106E6)+0.168805E4.*((-1).*((-4)+lambda1.^2).*( ...
  2+lambda1.^2.*((-1)+cos((1/45).*pi)))).^(1/2))+lambda1.^2.*(( ...
  -0.187616E8)+0.614331E6.*((-1).*((-4)+lambda1.^2).*(2+lambda1.^2.* ...
  ((-1)+cos((1/45).*pi)))).^(1/2))+lambda1.^4.*((-0.104873E9)+ ...
  0.146009E9.*((-1).*((-4)+lambda1.^2).*(2+lambda1.^2.*((-1)+cos(( ...
  1/45).*pi)))).^(1/2))+(-0.126711E9).*((-1).*((-4)+lambda1.^2).*(2+ ...
  lambda1.^2.*((-1)+cos((1/45).*pi)))).^(1/2))).*(0.109471E4+0.1E1.* ...
  lambda1.^4+lambda1.^2.*((-0.550023E3)+0.202612E2.*((-1).*((-4)+ ...
  lambda1.^2).*(2+lambda1.^2.*((-1)+cos((1/45).*pi)))).^(1/2))+( ...
  -0.278604E4).*((-1).*((-4)+lambda1.^2).*(2+lambda1.^2.*((-1)+cos(( ...
  1/45).*pi)))).^(1/2)).^(-1).*(0.412518E3+(-0.1E1).*lambda1.^2+( ...
  -0.202612E2).*((-1).*((-4)+lambda1.^2).*(2+lambda1.^2.*((-1)+cos(( ...
  1/45).*pi)))).^(1/2)).^(-2);
plot(-Pvslambda*145.08,lambda1,'r','linewidth',2)
hold on

    case 2
        C1 = 0.613441350420766 ; C2 =  0.000000002169605;% for corrected
Pvslambda =  lambda1.^(-3).*(C2.*(0.645074E1+(-0.168523E-4).*lambda1.^12+ ...
  lambda1.^6.*(0.542854E1+(-0.129347E1).*((2+lambda1.^2.*((-1)+cos(( ...
  1/18).*pi))).*(2+lambda1.^2.*((-1)+sin((13/90).*pi)))).^(1/2))+ ...
  lambda1.^10.*(0.398746E-2+(-0.182441E-3).*((2+lambda1.^2.*((-1)+ ...
  cos((1/18).*pi))).*(2+lambda1.^2.*((-1)+sin((13/90).*pi)))).^(1/2) ...
  )+lambda1.^8.*((-0.277467E0)+0.308342E-1.*((2+lambda1.^2.*((-1)+ ...
  cos((1/18).*pi))).*(2+lambda1.^2.*((-1)+sin((13/90).*pi)))).^(1/2) ...
  )+lambda1.^2.*((-0.195588E1)+0.297931E0.*((2+lambda1.^2.*((-1)+ ...
  cos((1/18).*pi))).*(2+lambda1.^2.*((-1)+sin((13/90).*pi)))).^(1/2) ...
  )+lambda1.^4.*((-0.142941E2)+0.926529E1.*((2+lambda1.^2.*((-1)+ ...
  cos((1/18).*pi))).*(2+lambda1.^2.*((-1)+sin((13/90).*pi)))).^(1/2) ...
  )+(-0.55518E1).*((2+lambda1.^2.*((-1)+cos((1/18).*pi))).*(2+ ...
  lambda1.^2.*((-1)+sin((13/90).*pi)))).^(1/2))+C1.*(0.230495E1+ ...
  0.195115E-3.*lambda1.^10+lambda1.^6.*(0.186627E1+(-0.285597E0).*(( ...
  2+lambda1.^2.*((-1)+cos((1/18).*pi))).*(2+lambda1.^2.*((-1)+sin(( ...
  13/90).*pi)))).^(1/2))+lambda1.^8.*((-0.395713E-1)+0.211229E-2.*(( ...
  2+lambda1.^2.*((-1)+cos((1/18).*pi))).*(2+lambda1.^2.*((-1)+sin(( ...
  13/90).*pi)))).^(1/2))+lambda1.^2.*((-0.997159E0)+0.106455E0.*((2+ ...
  lambda1.^2.*((-1)+cos((1/18).*pi))).*(2+lambda1.^2.*((-1)+sin(( ...
  13/90).*pi)))).^(1/2))+lambda1.^4.*((-0.558977E1)+0.522843E1.*((2+ ...
  lambda1.^2.*((-1)+cos((1/18).*pi))).*(2+lambda1.^2.*((-1)+sin(( ...
  13/90).*pi)))).^(1/2))+(-0.359838E1).*((2+lambda1.^2.*((-1)+cos(( ...
  1/18).*pi))).*(2+lambda1.^2.*((-1)+sin((13/90).*pi)))).^(1/2))).*( ...
  0.156267E3+0.1E1.*lambda1.^4+lambda1.^2.*((-0.901381E2)+ ...
  0.108259E2.*((2+lambda1.^2.*((-1)+cos((1/18).*pi))).*(2+ ...
  lambda1.^2.*((-1)+sin((13/90).*pi)))).^(1/2))+(-0.243957E3).*((2+ ...
  lambda1.^2.*((-1)+cos((1/18).*pi))).*(2+lambda1.^2.*((-1)+sin(( ...
  13/90).*pi)))).^(1/2)).^(-1).*(0.288411E0+(-0.42662E-2).* ...
  lambda1.^2+(-0.461855E-1).*((2+lambda1.^2.*((-1)+cos((1/18).*pi))) ...
  .*(2+lambda1.^2.*((-1)+sin((13/90).*pi)))).^(1/2)).^(-2);
plot(-Pvslambda*145.08,lambda1,'g','linewidth',2)
Material_fit = [-Pvslambda'*145.08 lambda1'];
save('Material_fit_6088.mat','Material_fit');
    otherwise
        disp('enter 1 for actual angles, and 2 for corrected angles')
end
end

hold on
P_exp = [0:3:24];
lambda1_exp =  [1
    1.0250
    1.0400
    1.0700
    1.1000
    1.1450
    1.2200
    1.3000
    1.4000]';
plot(P_exp,lambda1_exp,'bo:','linewidth',2)
h = legend('88,60 degrees fit','Corrected angles (85,58)','Experimental');
h.FontSize = 10;
h.FontWeight = 'bold';
xlim([0 30])
grid on
title ('Mooney Rivlin fit','Fontweight','bold','fontsize',12)
set(gca,'linewidth',2)
