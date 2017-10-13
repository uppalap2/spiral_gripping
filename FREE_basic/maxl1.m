function f=maxl1(x,alpha,beta,l,r)
lambda1=x(1);
lambda2=x(2);
delta=x(3);

theta=l/r*tand(alpha);phi=l/r*tand(beta);
f(1)=lambda1^2*cosd(alpha)^2+lambda2^2*sind(alpha)^2*((theta+delta)/theta)^2-1;
f(2)=lambda1^2*cosd(beta)^2+lambda2^2*sind(beta)^2*((phi+delta)/phi)^2-1;

alphan=atand((lambda2/lambda1)*(tand(alpha)+r*delta/l));
betan=atand((lambda2/lambda1)*(tand(beta)+r*delta/l));

% For second quadrant of design space
% f(3)=1+2*cotd(alphan)*cotd(betan); 

% For first quadrant of design space (till reaches 0 fiber angle)
f(3) = alphan; 
end
