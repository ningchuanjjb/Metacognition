clear all %#ok<CLALL>
close all

% raw = imread('C:\Users\JJB\Downloads\ZelkuPhoto2_3.jpg');
raw = imread('C:\Users\JJB\Downloads\pic_20240416B.png');

I = rgb2gray(raw);
% A = watershed(D);
% imshow(A);

fig1 = figure(1);
% BW = im2bw(D, thresh);
% BW = imbinarize(D, thresh);
% imshowpair(D,BW,'montage')
% imshow(BW);
% A = watershed(BW);
% 
% fig2 = figure(2);
% imshow(A);

% %线性滤波
% h=fspecial('motion',50,45);
% filtered=imfilter(I,h);

% %中值滤波
% J=imnoise(I,'salt & pepper',0.1);
% filtered=medfilt2(J);

% 均值滤波
h = fspecial('average',5*[1 1]);
filtered=imfilter(I,h);

% % 高斯滤波
% filtered = imgaussfilt(I,2);

% filtered = I;


% 黑白二值化
% thresh = graythresh(filtered);
% BW = im2bw(filtered, thresh);
BW = imbinarize(filtered, 'adaptive');
% BW = ~BW; 


% % imopen 开运算，先腐蚀再膨胀
% se = strel('disk',5);
% filtered2 = imopen(BW,se);

% imclose 闭运算，先膨胀再腐蚀
se = strel('disk',5);
filtered2 = imclose(BW,se);

% % bwskel 骨架化
% filtered3 = bwskel(filtered2,'MinBranchLength',1);

% %求边界
% filtered3 = bwperim(BW,8);

imshowpair(I,BW,'montage')

