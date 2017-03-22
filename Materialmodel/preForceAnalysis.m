% Overview of algorithm

% First get the shape with gravity for a given length(L)
% This is input to angle detect which outputs final angle and mean angle
% Get cylinder curvature and minimum length(L_min) needed check if L
% >=L_min then add a 1 on the next row.
function preForceAnalysis
clear all % needed 
L = 44e-2; %%%%%%%CHANGE THIS

Pressure = 8:18;
result_7088_26 = zeros(length(Pressure),5);
for i = 1:length(Pressure)
    alpha = 86.6470*pi/180;%%%%%%%%CHANGE THIS
    beta = 65.8020*pi/180;%%%%%%%%CHANGE THIS
    
    fitmatFile = 'Material_fit_6088.mat'; %Material_fit_6088.mat Material_fit_7088.mat
    
    initial_shape = getShapeGravity(Pressure(i), alpha , beta,fitmatFile,L);
    
    plot3(initial_shape(:,1),initial_shape(:,2),initial_shape(:,3),'r');
    axis equal
    grid on
    
    [TangVec, bias_angle] = angle_detect(initial_shape);
    CylCurv = initial_shape(end,13)/(sind(TangVec(end)-bias_angle))^2;
    MinLength = 2.5*pi*CylCurv^(-1)/sind(mean(TangVec)-bias_angle);
    check = MinLength>=L;
    
    result_7088_26(i,:) = [L Pressure(i) CylCurv^(-1) MinLength check];%%%%%%%%CHANGE THIS
end

% save('C:\Users\Naveen\Box Sync\git\spiral_gripping\Materialmodel\prototypes\result_7088.mat'...
%     ,'result_7088_26','-append'); % use append after first save only.%%%%%%%%CHANGE THIS

end





function [TangVec, bias_angle] = angle_detect(a)
CentralAxis = [0 0 1];

for i = 1:length(a)
    
    R(:,:,i) = reshape(a(i,4:12),3,3); % this should be transposed default is column wise
    TangVec(i) = atan2d(norm(cross(R(3,:,i)',CentralAxis)),dot(R(3,:,i)',CentralAxis));
end
bias_angle = atan2d(norm(cross([mean(a(:,1));mean(a(:,2));mean(a(:,3))],CentralAxis)),...
    dot([mean(a(:,1));mean(a(:,2));mean(a(:,3))],CentralAxis));
end