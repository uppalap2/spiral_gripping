% clear all

% load 20psiDia1.mat
% load initial_shape_20.mat
% a = initial_shape;
% a = c;
CentralAxis = [0 0 1];

for i = 1:length(a)
    
    R(:,:,i) = reshape(a(i,4:12),3,3); % this should be transposed default is column wise
    TangVec(i) = atan2d(norm(cross(R(3,:,i)',CentralAxis)),dot(R(3,:,i)',CentralAxis));
end

% kcyl = 1/1.3e-2;
kcyl = 1/1.7e-2;

klarge = kcyl.*sind(TangVec).^2;

% length(1.7e-2./sind(TangVec)')
% length(3.4e-2./sind(TangVec)')
% length(1./a(:,13))
% [1./a(:,13) 1.7e-2./sind(TangVec)' 3.4e-2./sind(TangVec)']

% % giving tangent data for a helix and finding the radius of curvature.
% 
% Tangent = [xtan; ytan ;ztan];
% CentralAxis = [0 0 1];
% 
% for i = 1:length(Tangent)
%     angleTan(i) = atan2d(norm(cross(Tangent(:,i),CentralAxis')),dot(Tangent(:,i),CentralAxis'));
% end
% 
% HelixCur = 1.3/(1.3^2+1^2);
% R_HelixCur = 1/HelixCur;