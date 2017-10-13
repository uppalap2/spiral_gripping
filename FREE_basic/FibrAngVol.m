function [deltaV,alphan,betan] = FibrAngVol(alpha,beta,volfrac,h)
% This function gives the updated fiber angles with change in volume
% fraction till the given volfrac is reached
% INPUTS:
% alpha, beta : Initial fiber angles (in degrees)
% volfrac : Volume fraction we want to actuate to
% h       : no of steps (linspace(0,volfrac,h))
% OUTPUTS:
% deltaV  :linspace(0,volfrac,h)
% alphan, betan : new fiber angles at each step

l = 10e-2;
r = 5e-3;
deltaV = linspace(0,volfrac,h);

options = optimoptions('fsolve','Display','iter','MaxIterations',1e3,'MaxFunctionEvaluations',1e4);

% Initialize alphan, betan
alphan = zeros(length(deltaV),1);
betan = zeros(length(deltaV),1);
alphan(1) = alpha;
betan(1) = beta;
for i = 2:length(deltaV)
    Lambdas = fsolve(@(x)getLambda(x,alphan(1),betan(1),l,r,deltaV(i)),[1 1 0],options);
    alphan(i)=atand((Lambdas(2)/Lambdas(1))*(tand(alphan(1))+r*Lambdas(3)/l));
    betan(i)=atand((Lambdas(2)/Lambdas(1))*(tand(betan(1))+r*Lambdas(3)/l));
end

plot(alphan,betan,'r*')
grid on


end

function f = getLambda(x,alpha,beta,l,r,vFrac)
lambda1 = x(1);
lambda2 = x(2);
delta = x(3);

theta=l/r*tand(alpha);phi=l/r*tand(beta);

f(1) = lambda1^2*cosd(alpha)^2+lambda2^2*sind(alpha)^2*((theta+delta)/theta)^2-1;
f(2) = lambda1^2*cosd(beta)^2+lambda2^2*sind(beta)^2*((phi+delta)/phi)^2-1;
f(3) = lambda2^2*lambda1-1-vFrac;
end