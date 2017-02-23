L12 = [40 45 50 55 57 60];
K12 = [32.8377 32.5772 32.8503 33.0244 33.0478 33.0549];

L14 = [25:5:60];
K14 = [60 46.2769 39.8818 40.1999 41.6797 42.0915 41.9215 41.8];

L16 = [25:5:60];
K16 = [64.8072 46.5043 50.3479 55.62222 55.0203 53.8268 53.8436 54.0486];

plot(L16,K16,'g','linewidth',2)
grid on
xlabel 'Length (in cm)'
ylabel 'Curvature'
title '14 psi'