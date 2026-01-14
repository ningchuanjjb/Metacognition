close all;clear all;clc;
% f = imread('C:\Users\JJB\Downloads\ZelkuPhoto2_4.jpg');
% f = imread('C:\Users\JJB\Downloads\pic_20240416B.png');
f = imread('C:\Users\JJB\Downloads\pic_20240416D.png');

[VG,A,PPG] = colorgrad(f);
ppg = im2uint8(PPG);
ppgf = 255 - ppg;
[M,N] = size(ppgf);

T=200;%200

ppgf1 = zeros(M,N);
for ii = 1:M
    for jj = 1:N
        if ppgf(ii,jj)<T
            ppgf1(ii,jj)=0;
        else
            ppgf1(ii,jj)=235/(255-T)*(ppgf(ii,jj)-T);
        end
    end
end
ppgf1 = uint8(ppgf1);

% T2 = 170;%140,170
% ppgf1(ppgf1>=T2) = 255;
% ppgf1(ppgf1<T2) = 0;

figure;
subplot(221);imshow(ppgf);
subplot(222);imshow(ppgf1);
subplot(223);imhist(ppgf);
subplot(224);imhist(ppgf1);
 
figure;imshow(ppgf1);