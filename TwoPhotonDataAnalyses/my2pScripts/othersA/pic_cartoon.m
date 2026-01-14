clear all %#ok<CLALL>
close all

% raw = imread('C:\Users\JJB\Downloads\ZelkuPhoto2_4.jpg');
raw = imread('C:\Users\JJB\Downloads\pic_20240416B.png');
% raw = imread('C:\Users\JJB\Downloads\pic_20240416C.png');

% fi_img = imbilatfilt(raw,100,20); %双边滤波
fi_img = imbilatfilt(raw,1000,20); %双边滤波
ed_img = 1-edge((rgb2gray(raw)),'canny',[0.07,0.1]); %边缘提取
carton_img = fi_img.*(uint8(repmat(ed_img,[1,1,3]))); %合并


carton_img = im2gray(carton_img);

T1 = 0;
T2 = 255;

carton_img(carton_img<T1) = 0;
carton_img(carton_img>=T2) = 255;


fig1 = figure(1);
% imshowpair(raw,carton_img,'montage')
imshow(carton_img);
% imshow(fi_img);

% exportgraphics(fig1, 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\ZelkuPhoto2.emf', 'BackgroundColor', 'none');