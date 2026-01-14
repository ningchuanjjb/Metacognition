close all
clear all
load template.mat

input_image_raw = gpuArray(single(template1));
backgroundColor = median(input_image_raw,'all');

I1_rotate = imrotate(input_image_raw,-3.444,'nearest','crop');
zero_index = (I1_rotate == 0);
I1_rotate = I1_rotate.*(~zero_index) + backgroundColor.*zero_index;


theta = -10:0.01:10;
% rho = zeros(1,length(theta),'single','gpuArray');
rho = zeros(1,length(theta));
for tempi=1:length(theta)
    temp_I_rotate = imrotate(input_image_raw,theta(tempi),'nearest','crop');
    rho(tempi) = corr2(temp_I_rotate,I1_rotate);
    a = 1;
end
[~,max_index] = max(rho);
theta_optimal = theta(max_index); 

I2_rotate = imrotate(input_image_raw,theta_optimal,'nearest','crop');
zero_index = (I2_rotate == 0);
I2_rotate = I2_rotate.*(~zero_index) + backgroundColor.*zero_index;


fig1 = figure('Name','Fig1, I1_rotate','NumberTitle','off');
set(gcf,'Position',[5 115+300 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(I1_rotate.^(1/3));
colormap(gray);
axis off equal

fig2 = figure('Name','Fig2, I2_rotate','NumberTitle','off');
set(gcf,'Position',[5+300 115+300 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(I2_rotate.^(1/3));
colormap(gray);
axis off equal

% I = zeros(100,100);
% I(25:75, 25:75) = 1;
% 
% % theta = 0:179;
% theta = 5;
% [R,xp] = radon(I,theta);

