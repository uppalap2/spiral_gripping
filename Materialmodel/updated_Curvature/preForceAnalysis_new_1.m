% Overview of algorithm

% First get the shape with gravity for a given length(L)
% This is input to angle detect which outputs final angle and mean angle
% Get cylinder curvature and minimum length(L_min) needed check if L
% >=L_min then add a 1 on the next row.
function result = preForceAnalysis_new_1
clear all % needed 
% clf
L = 50e-2; %%%%%%%CHANGE THIS

Pressure = 10:1:20;
result_6088_45 = zeros(length(Pressure),7);
Eb = 10e5; %  (.6134) for 60 88 (.4386) for 70,88 (.6029) 50,88
WpL = .0332; % .0332-50, .0326-60, .0345-70
option =  1; % -0 w/o g -1 w g - 2 grasping
for i = 1:length(Pressure)
    alpha = 87*pi/180;%%%%%%%%CHANGE THIS
    beta = 46*pi/180;%%%%%%%%CHANGE THIS
    
    fitmatFile = 'Material_fit_5088.mat'; %Material_fit_6088.mat Material_fit_7088.mat
    
    [initial_shape,lambda1] = getShapeGravity_1(Pressure(i), alpha , beta,fitmatFile,L,Eb,WpL,option);
    
    plot3(initial_shape(:,1),initial_shape(:,2),initial_shape(:,3),'c');
    axis equal
    grid on
    hold on
%     Microscribe_data_plots
    [TangVec] = angle_detect(initial_shape);
    CylCurv = initial_shape(end,13)/(sind(TangVec(end)))^2;
    MinLength = 2.5*pi*CylCurv^(-1)/sind(mean(TangVec));
    check = MinLength>=L;
    Ract = 13/64*.0254;
    result_6088_45(i,:) = [L Pressure(i) CylCurv^(-1)-0*Ract MinLength check TangVec(end) lambda1];%%%%%%%%CHANGE THIS
end

% save('C:\Users\Naveen\Box Sync\git\spiral_gripping\Materialmodel\prototypes\result_5088_Eb.mat'...
%     ,'result_5088_40','-append'); % use append after first save only.%%%%%%%%CHANGE THIS

result  = result_6088_45;
end





function TangVec = angle_detect(a)
% works perfect for without gravity
% CentralAxis = [-.5489 -.0058 .8359];% 50 cm
x = a(end-5:end,6);
y = a(end-5:end,9);
z = a(end-5:end,12);
 plot3(x/10,y/10,z/10,'r.')
 hold on
 
% plot3(a(:,1),a(:,2),a(:,3),'k')

 A=[x-mean(x), y-mean(y), z-mean(z)];

[~,~,V]=svd(A,0);

CentralAxis = V(:,end); %Approximately
plot3([0 CentralAxis(1)/3],[0 CentralAxis(2)/3 ], [0 CentralAxis(3)/3],'k-')
axis equal 
grid on
% CentralAxis = [1 0 1];
for i = 1:length(a)
    
    R(:,:,i) = reshape(a(i,4:12),3,3); % this should be transposed default is column wise
    TangVec(i) = atan2d(norm(cross(R(3,:,i)',CentralAxis)),dot(R(3,:,i)',CentralAxis));
end

end