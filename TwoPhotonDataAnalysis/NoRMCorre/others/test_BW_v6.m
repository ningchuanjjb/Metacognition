close all
% input_image = gather(Y_raw(:,:,1));
% input_image = gather(mean(Y_raw(:,:,1:2),3));
numFrames = 1;
input_image_raw = Y_raw(:,:,1:numFrames);
input_image = reshape(input_image_raw,[512,512*numFrames]);

BW_raw = gather(input_image > 200);
% BW_raw = imbinarize(input_image,200);
BW = bwpack(BW_raw);

nhood = [1 1];
BW1 = imerode(BW,nhood,'ispacked',512);
BW1_2 = imdilate(BW1,nhood,'ispacked');

nhood = [1; 1];
BW2 = imerode(BW,nhood,'ispacked',512);
BW2_2 = imdilate(BW2,nhood,'ispacked');

nhood = [0 1; 1 0];
BW3 = imerode(BW,nhood,'ispacked',512);
BW3_2 = imdilate(BW3,nhood,'ispacked');

nhood = [1 0; 0 1];
BW4 = imerode(BW,nhood,'ispacked',512);
BW4_2 = imdilate(BW4,nhood,'ispacked');

% closeBW = BW1_2 | BW2_2;
% closeBW = BW1_2 | BW2_2 | BW3_2 | BW4_2;
closeBW = bitor(bitor(bitor(BW1_2,BW2_2),BW3_2),BW4_2);
% closeBW = bitor(BW1_2,BW2_2);

% tempBW = bwunpack(BW2,512);


se = strel('rectangle',[3 3]);
closeBW_dilate = bwunpack(imdilate(closeBW,se,'ispacked'),512);



% backgroundColor = 0;
backgroundColor = median(input_image,'all');
close_dilate = input_image;
% close_dilate(~closeBW_dilate) = backgroundColor;
close_dilate = close_dilate.*closeBW_dilate + backgroundColor.*(~closeBW_dilate);
close_dilate = reshape(close_dilate,[512,512,numFrames]);
close_dilate = mean(close_dilate,3);

fig1 = figure('Name','Fig1, input_image','NumberTitle','off');
set(gcf,'Position',[5 115+300 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(mean(input_image_raw,3).^(1/3));
% imagesc(input_image.^(1/3));
% imagesc(input_image);
colormap(gray);
axis off equal

fig2 = figure('Name','Fig2, BW','NumberTitle','off');
set(gcf,'Position',[5+300 115+300 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(mean(reshape(BW_raw,[512,512,numFrames]),3));
% imagesc(BW_raw);
colormap(gray);
axis off equal

fig3 = figure('Name','Fig3, closeBW_dilate','NumberTitle','off');
set(gcf,'Position',[5+300*2 115+300 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(mean(reshape(closeBW_dilate,[512,512,numFrames]),3));
% imagesc(closeBW_dilate);
colormap(gray);
axis off equal

fig4 = figure('Name','Fig4, close_dilate','NumberTitle','off');
set(gcf,'Position',[5+300*3 115+300 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(close_dilate.^(1/3));
% imagesc(close_dilate);
colormap(gray);
axis off equal

% fig5 = figure('Name','Fig5, tempBW','NumberTitle','off');
% set(gcf,'Position',[5+300*4 115+300 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
% set (gca,'position',[0.01,0.01,0.98,0.98] )
% imagesc(tempBW);
% colormap(gray);
% axis off equal
