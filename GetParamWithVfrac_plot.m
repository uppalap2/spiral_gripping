% Plot file for volume frac versus lambda1 and torsion(i.e delta/l)
close all
alpha = 88;
beta = 65;
l = .159;
r = .0048;

vfrac = 0:.01:.3;
params = zeros(length(vfrac),4);
for i = 1:length(vfrac)

params(i,:) = [vfrac(i) GetParamWithVfrac(alpha, beta, l, r, vfrac(i))];


end


params_changed = [params(:,1) (params(:,2)-1)/(2*r) params(:,3) -params(:,4)/l];
% plot wrt curvatures and torsion
figure
plot(vfrac,(params(:,2)-1)/(2*r),'r','linewidth', 2) % vol frac versus curvature
xlabel 'Volume fraction \Delta V/V_0'
ylabel 'Curvature'
grid on


figure

plot(vfrac,-params(:,4)/l,'g','linewidth', 2)
xlabel 'Volume fraction \Delta V/V_0'
ylabel 'Torsion'
grid on

% for 10 psi 6.7034 13.5069
% 12PSI      8.9947 17.1582
% for 14 psi 10.7905 20.8017   for 16 psi 13.0044 25.4410
% for 18 psi 15.3885 29.7957   for 20 psi 18.1132 34.9487
% for 22 psi 23.0827 43.0181

% plot(10:2:22,[6.7034 8.9947 10.7905 13.0044 15.3885 18.1132 23.0827],'g-')

boxplot(params_changed(2:end,4)./params_changed(2:end,2))