close all

% MATLAB Code | Butterworth bandstop Filter
	
% Reading input image : input_image
input_image = Y_raw(:,:,1);
% input_image = mean(Y_raw(:,:,1:2),3);
f_input_image = fftshift(fft2(input_image));

% Saving the size of the input_image in pixels-
% M : no of rows (height of the image)
% N : no of columns (width of the image)
[M, N] = size(input_image);

% Getting Fourier Transform of the input_image
% using MATLAB library function fft2 (2D fast fourier transform)
FT_img = fftshift(fft2(double(input_image)));

% Assign the order value
n_low = 30; % one can change this value accordingly
n_high = 30;

% Assign Cut-off Frequency
D0_low = 25; % one can change this value accordingly
D0_high = 20; 

if_low = 1;
if_high = 1;

% Designing filter
u = 0:(M-1);
v = 0:(N-1);
idx = find(u > M/2);
u(idx) = u(idx) - M;
idy = find(v > N/2);
v(idy) = v(idy) - N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy of v
% and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);

% Calculating Euclidean Distance
D = sqrt(U.^2 + V.^2);

% determining the filtering mask
H_low = 1-1./(1 + (D./D0_low).^(2*n_low));
H_low = fftshift(H_low);

H_high = 1-1./(1 + (D0_high./D).^(2*n_high));
H_high = fftshift(H_high);
% Convolution between the Fourier Transformed
% image and the mask
% G = H_low.*FT_img.*H_high;
if if_low == 1 && if_high == 1
    %H = H_low.*H_high;
    H = H_low + H_high;
elseif if_low == 0 && if_high == 1
    H = H_high;
elseif if_low == 1 && if_high == 0
    H = H_low;
end
G = H.*FT_img;

% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function
% ifft2 (2D inverse fast fourier transform)
output_image = real(ifft2(double(ifftshift(G))));
output_image(output_image<0) = 0;

% Displaying Input Image and Output Image
fig1 = figure('Name','Fig1, input_image','NumberTitle','off');
set(gcf,'Position',[5 115+300 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(input_image.^(1/3));
colormap(gray);
axis off equal

fig2 = figure('Name','Fig2, output_image','NumberTitle','off');
set(gcf,'Position',[5+300 115+300 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(output_image.^(1/3));
colormap(gray);
axis off equal

fig3 = figure('Name','Fig3, F','NumberTitle','off');
set(gcf,'Position',[5 35 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(log(1+abs(f_input_image)));
colormap(gray);
axis off equal

fig4 = figure('Name','Fig4, H_low','NumberTitle','off');
set(gcf,'Position',[5+300*1 35 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(abs(H_low));
colormap(gray);
axis off equal

fig5 = figure('Name','Fig5, H_high','NumberTitle','off');
set(gcf,'Position',[5+300*2 35 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(abs(H_high));
colormap(gray);
axis off equal

fig6 = figure('Name','Fig6, H','NumberTitle','off');
set(gcf,'Position',[5+300*3 35 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(abs(H));
colormap(gray);
axis off equal

fig7 = figure('Name','Fig7, G','NumberTitle','off');
set(gcf,'Position',[5+300*4 35 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(log(1+abs(G)));
colormap(gray);
axis off equal
