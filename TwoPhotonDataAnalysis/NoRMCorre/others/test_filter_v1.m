close all

Y_test_single = Y_raw(:,:,1);
Y_test_double= mean(Y_raw(:,:,1:2),3);
Y_test_mean = mean(Y_raw,3);

f_test_single = fftshift(fft2(Y_test_single));

h = fspecial('average',3);
% h = fspecial('disk',3);
% h = fspecial('gaussian',3,2);
% h = fspecial('laplacian',0.2 ;
Y_test_single_filter = imfilter(Y_test_single,h);
Y_test_double_filter = imfilter(Y_test_double,h);
Y_test_mean_filter = imfilter(Y_test_mean,h);

% sigma = 0.25;
% Y_test_single_filter = imgaussfilt(Y_test_single,sigma,'FilterSize',3,'FilterDomain','frequency');
% Y_test_double_filter = imgaussfilt(Y_test_double,sigma,'FilterSize',3,'FilterDomain','frequency');
% Y_test_mean_filter = imgaussfilt(Y_test_mean,sigma,'FilterSize',3,'FilterDomain','frequency');


% pad = [3,3];
% Y_test_single_filter = medfilt2(Y_test_single,pad);
% Y_test_double_filter = medfilt2(Y_test_double,pad);
% Y_test_mean_filter = medfilt2(Y_test_mean,pad);

% pad = [3,3];
% Y_test_single_filter = wiener2(Y_test_single,pad);
% Y_test_double_filter = wiener2(Y_test_double,pad);
% Y_test_mean_filter = wiener2(Y_test_mean,pad);

Y_test_single_filter(Y_test_single_filter<0) = 0;
Y_test_double_filter(Y_test_double_filter<0) = 0;
Y_test_mean_filter(Y_test_mean_filter<0) = 0;

fig1 = figure('Name','Fig1, Y_test_single','NumberTitle','off');
set(gcf,'Position',[5 115+320 320 320]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(Y_test_single.^(1/3));
colormap(gray);
axis off equal

fig2 = figure('Name','Fig2, Y_test_double','NumberTitle','off');
set(gcf,'Position',[5+320 115+320 320 320]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(Y_test_double.^(1/3));
colormap(gray);
axis off equal

fig3 = figure('Name','Fig3, Y_test_mean','NumberTitle','off');
set(gcf,'Position',[5+320*2 115+320 320 320]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(Y_test_mean.^(1/3));
colormap(gray);
axis off equal

fig4 = figure('Name','Fig4, Y_test_single_filter','NumberTitle','off');
set(gcf,'Position',[5 35 320 320]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(Y_test_single_filter.^(1/3));
colormap(gray);
axis off equal

fig5 = figure('Name','Fig5, Y_test_double_filter','NumberTitle','off');
set(gcf,'Position',[5+320 35 320 320]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(Y_test_double_filter.^(1/3));
colormap(gray);
axis off equal

fig6 = figure('Name','Fig6, Y_test_mean_filter','NumberTitle','off');
set(gcf,'Position',[5+320*2 35 320 320]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(Y_test_mean_filter.^(1/3));
colormap(gray);
axis off equal

fig7 = figure('Name','Fig7, f_test_single','NumberTitle','off');
set(gcf,'Position',[5 115+320 320 320]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(log(1+abs(f_test_single)));
colormap(gray);
axis off equal