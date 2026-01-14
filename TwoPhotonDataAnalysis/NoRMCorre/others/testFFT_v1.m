clear 
close all

image_length = 210;
d1 = 225;
d2 = 225;
if_FTpad = 1;

fig1 = figure(1);
set(gcf,'Position',[5+image_length*0 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
    
%P = peaks(20);
%X = repmat(P,[5 10]);
X = rgb2gray(im2double(imread('lena.jpg')));
% X = X(:,randperm(225));

% q = 0.97405;
q = 0.973105;
X = (1-q)*X + q*rand(225);
imagesc(X)


fig2 = figure(2);
set(gcf,'Position',[5+image_length*1 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
Y = fft2(X);
% imagesc(abs(fftshift(Y)))
% imagesc(real(fftshift(Y)))
imagesc(imag(fftshift(Y)))

fig22 = figure(22);
set(gcf,'Position',[5+image_length*1 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
% Y = fft2(X);
imagesc(real(fftshift(sign(Y))))

fig23 = figure(23);
set(gcf,'Position',[5+image_length*1 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
% Y = fft2(X);
imagesc(imag(fftshift(sign(Y))))

fig3 = figure(3);
set(gcf,'Position',[5+image_length*2 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
X2 = ifft2(Y);
imagesc(X2)

fig4 = figure(4);
set(gcf,'Position',[5+image_length*3 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
X3 = real(ifft2(sign(Y)));
imagesc(X3)

fig5 = figure(5);
set(gcf,'Position',[5+image_length*4 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
% X_offset = [X(:,51:225) zeros(225,50)];
X_offset = [median(X,'all')*ones(225,50) X(:,1:175)];
X_offset = (1-q)*X_offset + q*rand(225);

% X_offset = [X(1:175,:); zeros(50,225)];
% X_offset = [zeros(50,225); X(1:175,:)];
imagesc(X_offset)

fig6 = figure(6);
set(gcf,'Position',[5+image_length*5 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
X_offset;
X;
Y_offset = fft2(X_offset);
imagesc(abs(fftshift(Y_offset)))

Y_prod = Y.*conj(Y_offset);% To compute image cross correlation in complex domain

if_FTpad = 1;
if if_FTpad == 1
    FTpad_coeff = 2;
    Y_pad = FTpad(Y_prod,ceil(FTpad_coeff*[d1,d2]));
else
    FTpad_coeff = 1;
    Y_pad = Y_prod;
end
    
Y_pad = sign(Y_pad);
CC = ifft2(Y_pad);
CCabs = abs(CC);

[MMM, III] = max(CCabs, [], 'all','linear'); %#ok<*ASGLU>
[row_shift,col_shift] = ind2sub(size(CCabs),III);
Nr2 = ifftshift(-fix(d1):ceil(d1)-1);
Nc2 = ifftshift(-fix(d2):ceil(d2)-1);
% Nr2 = (1:size(Y_pad,1))-1;
% Nc2 = (1:size(Y_pad,2))-1;

row_shift = ceil(single(Nr2(row_shift)/FTpad_coeff));
col_shift = ceil(single(Nc2(col_shift)/FTpad_coeff));
fprintf('row_shift = %d, col_shift = %d.\n',row_shift,col_shift);    

fig7 = figure(7);
set(gcf,'Position',[5+image_length*6 35 image_length image_length]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] ), colormap(gray)
axis off equal
imagesc(CCabs)

a = 1;



