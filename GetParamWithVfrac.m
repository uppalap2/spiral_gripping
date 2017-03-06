function [params,alphan, betan] = GetParamWithVfrac(alpha, beta, l, r, vfrac)
% takes in alpha, beta, lenth and radius of FREE along with volume fraction
% to give the new parameters and the final angles should check if has
% reached its lock configuration from angles.
% One more point of concern is after certain volfrac FREEs which extend and
% rotate change the angles to second quadrant, there have to do in two
% steps
inputs =  nargin;
if inputs < 5
    alpha = 88;
    beta = 60;
    l = .159;
    r = .0048;
    vfrac = .3;
end



%initial lambda1,lambda2,delta condition
x0=[1 1 0];
options=optimoptions('fsolve','Display','iter','MaxFunEvals',100000,'MaxIter',100000);

x = fsolve(@(x)finallambda(x,alpha,beta,l,r,vfrac),x0,options);
% [x] = fsolve(@(x)maxl1(x,alpha,beta,l,r),x0,options);
lambda1 = x(1); lambda2 = x(2); delta = x(3);
alphan = atand((lambda2/lambda1)*(tand(alpha)+r*delta/l));
betan = atand((lambda2/lambda1)*(tand(beta)+r*delta/l));
params = [x(1) x(2) x(3)]; 
end

function f= finallambda(x,alpha,beta,l,r,vfrac)
lambda1=x(1);
lambda2=x(2);
delta=x(3);

theta=tand(alpha)*l/r;
phi=tand(beta)*l/r;

f(1)=lambda1^2*cosd(alpha)^2+lambda2^2*sind(alpha)^2*(theta+delta)^2/theta^2-1;
f(2)=lambda1^2*cosd(beta)^2+lambda2^2*sind(beta)^2*(phi+delta)^2/phi^2-1;
f(3)=vfrac+1-lambda1*lambda2^2;
end

% function f = maxl1(x,alpha,beta,l,r)
% lambda1 = x(1);
% lambda2 = x(2);
% delta = x(3);
% % lambda2=(sign(beta)*sqrt(1-cosd(beta)^2*lambda1^2)*cosd(alpha)...
% %     -sign(alpha)*sqrt(1-cosd(alpha)^2*lambda1^2)*cosd(beta))/(sind(alpha-beta));
% % f=-lambda2^2*lambda1+1;
% theta = l/r*tand(alpha);phi = l/r*tand(beta);
% f(1) = lambda1^2*cosd(alpha)^2+lambda2^2*sind(alpha)^2*((theta+delta)/theta)^2-1;
% f(2) = lambda1^2*cosd(beta)^2+lambda2^2*sind(beta)^2*((phi+delta)/phi)^2-1;
% alphan = atand((lambda2/lambda1)*(tand(alpha)+r*delta/l));
% betan = atand((lambda2/lambda1)*(tand(beta)+r*delta/l));
% f(3) = 1+2*cotd(alphan)*cotd(betan);    
% end

