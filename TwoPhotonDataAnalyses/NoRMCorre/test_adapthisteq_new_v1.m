% figure(1);
% set(gcf,'Position',[35+400 35+0 850 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
% set (gca,'position',[0.01,0.01,0.98,0.98])
% template_n11n = template_in./max(template_in,[],'all');
%
% imshow(template_n11n.^0.4);
%
% colormap(gray);
% axis off equal

gather_template_in; % This is current template.
temp_M1_final = M1_final; %#ok<*VUNUS> % This is current template with correction
% temp_gather_template_n11n = gather_template_n11n; % This is existing template.
temp_gather_template_n11n = temp_reg;

% optimal_NumTiles_histeq = 32;%8
% optimal_nBins_histeq = 4096*2;%256
% optimal_alpha_histeq = 0.11777;%0.4,0.11777
% ClipLimit = 0.1;%0.1
% 
% temp_template_n11n = adapthisteq(initial_template_n11n,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,...
%     'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq,'ClipLimit',ClipLimit);
% 
% gather_template_n11n = rescale(temp_template_n11n,min_template_in_raw,max_template_in_raw); %#ok<*NASGU>
% 
% temp_template_in = adapthisteq(initial_template_in,'NumTiles',[1,1]*optimal_NumTiles_histeq,'NBins',optimal_nBins_histeq,...
%     'Range','original','Distribution','rayleigh','Alpha',optimal_alpha_histeq,'ClipLimit',ClipLimit);
% gather_template_in = rescale(temp_template_in,min_template_in_raw,max_template_in_raw);
% 
% max_gather_template_n11n = max(gather_template_n11n,[],'all');
% max_gather_template_in = max(gather_template_in,[],'all');


% rho_0 = corr2(gather_template_in,gather_template_n11n);

close all

if_plotHist = 0;

temp_max = max_gather_template_in;

figure(1);
set(gcf,'Position',[35+400 35+0 850 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98])
I = gather_template_in;
I = I ./ temp_max;
imshow(I);
colormap(gray);
axis off equal

if if_plotHist == 1
    figure(12)
    I = gather_template_in;
    I = I ./ temp_max;
    imhist(I)
end

figure(2);
set(gcf,'Position',[35+400 35+0 850 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98])
I = temp_M1_final;
I = I ./ temp_max;
imshow(I);
colormap(gray);
axis off equal

if if_plotHist == 1
    figure(22)
    I = temp_M1_final;
    I = I ./ temp_max;
    imhist(I)
end
figure(3);
set(gcf,'Position',[35+400 35+0 850 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98])
I = temp_gather_template_n11n;
I = I ./ temp_max;
imshow(I);
colormap(gray);
axis off equal

if if_plotHist == 1
    figure(32)
    I = temp_gather_template_n11n;
    I = I ./ temp_max;
    imhist(I)
end

