% Image processing to get angles

% prototype 1
function [average1,average2,angle_center] = getfiberangle(image)
prompt = 'Number of points? ';
x = input(prompt);
alpha = zeros(x,1);
beta = zeros(x,1);
true = imread(image);
imshow(true);
gray_image = rgb2gray(true);
imshow(gray_image);
BW = imbinarize(gray_image,.65);
% get angle of the centerline first for refrence)
fprintf (' Select points for horizontal axis \n');


[pc1,pc2] = rbline;
% [pc1,pc2] = getline;
 angle_center = rad2deg(atan((pc2(2)-pc1(2))/(pc2(1)-pc1(1))));

fprintf (' Select points for measuring angle \n');
 
for i = 1:x
    [alpha(i),beta(i)] = getfiberangle_instance;
end

average1 = mean(alpha);
average2 = mean(beta);

end
function [angle1,angle2] = getfiberangle_instance
% image is a string that is like 'pic.jpg'




[p1,p2] = rbline;
[p3,p4] = rbline;

 angle1 = rad2deg(atan((p2(2)-p1(2))/(p2(1)-p1(1))));
 angle2 = rad2deg(atan((p4(2)-p3(2))/(p4(1)-p3(1))));
end