function TangVec = angle_detect(a)
% works perfect for without gravity
% CentralAxis = [-.5489 -.0058 .8359];% 50 cm
x = a(:,6);
y = a(:,9);
z = a(:,12);
 plot3(x/10,y/10,z/10,'r.')
 hold on
 
plot3(a(:,1),a(:,2),a(:,3),'k')
axis equal 
grid on
 A=[x-mean(x), y-mean(y), z-mean(z)];

[~,~,V]=svd(A,0);

CentralAxis = V(:,end); %Approximately
plot3([0 CentralAxis(1)/3],[0 CentralAxis(2)/3 ], [0 CentralAxis(3)/3],'k-')
% CentralAxis = [1 0 1];
for i = 1:length(a)
    
    R(:,:,i) = reshape(a(i,4:12),3,3); % this should be transposed default is column wise
    TangVec(i) = atan2d(norm(cross(R(3,:,i)',CentralAxis)),dot(R(3,:,i)',CentralAxis));
end

end