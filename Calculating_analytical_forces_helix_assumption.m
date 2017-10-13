%% This code is written for estimating the rotation matrix from the shape of the FREE
close all
% rad = cylinder radius ; curv = curvature; tor = torsion
curv = 1;
tor = 1;
rad = curv/(curv^2+tor^2);
pitch = tor/(curv^2+tor^2);
t = 0:2.25*pi/100:2.25*pi;
x = rad*cos(t);
y = rad*sin(t);
z = pitch*t;
%% check the spiral shape
plot3(x,y,z,'m*');
axis equal
grid on

r = [x' y' z'];

rdot = diff(r)/(2.25*pi/100); 
rdot_square = cumsum(rdot.^2,2); 
Tangent = rdot./sqrt(rdot_square(:,3));

% check the tangents directions

hold on
quiver3(x(1:100)',y(1:100)',z(1:100)',...
    Tangent(:,1)/10,Tangent(:,2)/10,Tangent(:,3)/10)

Tdot = diff(Tangent)/(2.25*pi/100);
Tdot_square =  cumsum(Tdot.^2,2); 

Normal = Tdot./sqrt(Tdot_square(:,3));

% check the normals direction

hold on
quiver3(x(1:99)',y(1:99)',z(1:99)',...
    Normal(:,1)/10,Normal(:,2)/10,Normal(:,3)/10)

Binormal = cross(Tangent(1:99,:),Normal(1:99,:));

% chcek the binormals direction

hold on
quiver3(x(1:99)',y(1:99)',z(1:99)',...
    Binormal(:,1)/10,Binormal(:,2)/10,Binormal(:,3)/10)

R = [Normal Binormal Tangent(1:99,:)]; % after reshaping take transpose for getting original R matrix

Rdot = diff(R)/(2.25*pi/100); 

Rotdot = [Rdot(:,1:3);Rdot(:,4:6);Rdot(:,7:9)];