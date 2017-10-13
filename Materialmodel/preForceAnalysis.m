% Overview of algorithm

% First get the shape with gravity for a given length(L)
% This is input to angle detect which outputs final angle and mean angle
% Get cylinder curvature and minimum length(L_min) needed check if L
% >=L_min then add a 1 on the next row.
function result = preForceAnalysis
clear all % needed 
L = 45e-2; %%%%%%%CHANGE THIS

Pressure = 8:18;
result_5088_45 = zeros(length(Pressure),5);
Eb = 7e5; %  (.6134) for 60 88 (.4386) for 70,88 (.6029) 50,88
WpL = .0332; % .0332-50, .0326-60, .0345-70
for i = 1:length(Pressure)
    alpha = 86.88*pi/180;%%%%%%%%CHANGE THIS
    beta = 45.49*pi/180;%%%%%%%%CHANGE THIS
    
    fitmatFile = 'Material_fit_5088.mat'; %Material_fit_6088.mat Material_fit_7088.mat
    
    initial_shape = getShapeGravity(Pressure(i), alpha , beta,fitmatFile,L,Eb,WpL);
    
    plot3(initial_shape(:,1),initial_shape(:,2),-initial_shape(:,3),'c');
    axis equal
    grid on
    hold on
    
    [TangVec, bias_angle] = angle_detect(initial_shape);
    CylCurv = initial_shape(end,13)/(sind(TangVec(end)-bias_angle*0))^2;
    MinLength = 2.5*pi*CylCurv^(-1)/sind(mean(TangVec)-bias_angle*0);
    check = MinLength>=L;
    
    result_6088_45(i,:) = [L Pressure(i) CylCurv^(-1) MinLength check];%%%%%%%%CHANGE THIS
end

% save('C:\Users\Naveen\Box Sync\git\spiral_gripping\Materialmodel\prototypes\result_5088_Eb.mat'...
%     ,'result_5088_40','-append'); % use append after first save only.%%%%%%%%CHANGE THIS

result  = result_6088_45;
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