clear all
close all


for i = 2:2:26
    length = 24 + i;
    len = num2str(length);
%     file_name = strcat('result_',len,'cm.mat');
     load(['result_' num2str(len) 'cm.mat'])
    file = strcat('result_',len);
    result = eval(file);
    
    plot(result(:,2), result(:,3),'-.', 'linewidth',2,'color', rand(1,3));
    hold on
    legendInfo{i/2} = [num2str(length) ]; % or whatever is appropriate
    
end
xlabel 'Pressure (in psi)'
ylabel 'Cylinder radius (in m)'
legend(legendInfo)
grid on