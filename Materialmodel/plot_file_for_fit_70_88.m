

lambda1 = [1.000:.01:1.34];






        C1 = 0.4386 ; C2 =  0;% for corrected
Pvslambda =  lambda1.^(-3).*(0.148534E3+(-0.1E1).*lambda1.^2+(-0.20813E2).*(4+( ...
  -0.685783E0).*lambda1.^2+0.230851E-2.*lambda1.^4).^(1/2)).^(-2).*( ...
  0.577572E3+0.1E1.*lambda1.^4+(-0.103047E4).*(4+(-0.685783E0).* ...
  lambda1.^2+0.230851E-2.*lambda1.^4).^(1/2)+lambda1.^2.*(( ...
  -0.198045E3)+0.20813E2.*(4+(-0.685783E0).*lambda1.^2+0.230851E-2.* ...
  lambda1.^4).^(1/2))).^(-1).*(C1.*((-0.241032E7)+(-0.249144E2).* ...
  lambda1.^10+0.430037E7.*(4+(-0.685783E0).*lambda1.^2+0.230851E-2.* ...
  lambda1.^4).^(1/2)+lambda1.^4.*(0.605986E7+(-0.612037E7).*(4+( ...
  -0.685783E0).*lambda1.^2+0.230851E-2.*lambda1.^4).^(1/2))+ ...
  lambda1.^2.*(0.61986E6+(-0.579044E5).*(4+(-0.685783E0).* ...
  lambda1.^2+0.230851E-2.*lambda1.^4).^(1/2))+lambda1.^8.*( ...
  0.111019E5+(-0.518543E3).*(4+(-0.685783E0).*lambda1.^2+ ...
  0.230851E-2.*lambda1.^4).^(1/2))+lambda1.^6.*((-0.114015E7)+ ...
  0.154042E6.*(4+(-0.685783E0).*lambda1.^2+0.230851E-2.*lambda1.^4) ...
  .^(1/2)))+C2.*((-0.653129E7)+0.925926E0.*lambda1.^12+0.628398E7.*( ...
  4+(-0.685783E0).*lambda1.^2+0.230851E-2.*lambda1.^4).^(1/2)+ ...
  lambda1.^4.*(0.138163E8+(-0.978227E7).*(4+(-0.685783E0).* ...
  lambda1.^2+0.230851E-2.*lambda1.^4).^(1/2))+lambda1.^2.*( ...
  0.116373E7+(-0.156904E6).*(4+(-0.685783E0).*lambda1.^2+ ...
  0.230851E-2.*lambda1.^4).^(1/2))+lambda1.^8.*(0.731022E5+( ...
  -0.715607E4).*(4+(-0.685783E0).*lambda1.^2+0.230851E-2.* ...
  lambda1.^4).^(1/2))+lambda1.^10.*((-0.481359E3)+0.192713E2.*(4+( ...
  -0.685783E0).*lambda1.^2+0.230851E-2.*lambda1.^4).^(1/2))+ ...
  lambda1.^6.*((-0.304356E7)+0.654445E6.*(4+(-0.685783E0).* ...
  lambda1.^2+0.230851E-2.*lambda1.^4).^(1/2))));
plot(Pvslambda*145.08,lambda1,'g','linewidth',2)
hold on
% Material_fit = [-Pvslambda'*145.08 lambda1'];
% save('Material_fit_6088.mat','Material_fit');




P_exp = [0:2:18];
lambda1_exp =  [  1.0000    1.0222    1.0667    1.1000    1.1111    1.1333   ...
                    1.1778    1.2222    1.2556    1.3222];
plot(P_exp,lambda1_exp,'bo:','linewidth',2)
h = legend('Corrected angles (65,86)','Experimental');
h.FontSize = 10;
h.FontWeight = 'bold';
xlim([0 max(P_exp)])
grid on
title ('Mooney Rivlin fit','Fontweight','bold','fontsize',12)
set(gca,'linewidth',2)