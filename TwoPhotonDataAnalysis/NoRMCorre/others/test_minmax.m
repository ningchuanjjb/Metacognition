% clear 
close all

image_length = 500;

fig1 = figure(1);
set(gcf,'Position',[5+image_length*0 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
    
X1 = Ytm(:,:,1);

imagesc(X1.^(1/3));

fig2 = figure(2);
set(gcf,'Position',[5+image_length*1 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
    
X2 = Mf(:,:,1);

imagesc(X2.^(1/3));