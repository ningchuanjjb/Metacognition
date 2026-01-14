close all

s = load(fullFileName_Fall,'stat');
roi_stats_raw = s.stat;
roi_stats = roi_stats_raw(cellIndex);


tempIndex = 886 + 1;%69
temp_roi_stat = roi_stats{cellIndex==tempIndex};

% tempIndex = 908 + 1;
% temp_roi_stat = roi_stats{tempIndex};
 
% temp_roi_med = zeros(1,2);
% temp_roi_med(1) = median(double(temp_roi_stat.xpix));
% temp_roi_med(2) = median(double(temp_roi_stat.ypix));

temp_roi_npix = double(temp_roi_stat.npix);

temp_roi_xpix = double(temp_roi_stat.xpix);
temp_roi_ypix = double(temp_roi_stat.ypix);

[temp_x_min,temp_x_max] = bounds(temp_roi_xpix);
[temp_y_min,temp_y_max] = bounds(temp_roi_ypix);

a = 1;
temp_roi_xpix;
temp_roi_ypix;



temp_diff = temp_x_max-temp_x_min;
if temp_diff < 5
    temp_x_min = temp_x_min - 2;
    temp_x_max = temp_x_max + 2;
end
temp_diff = temp_y_max-temp_y_min;
if temp_diff < 5
    temp_y_min = temp_y_min - 2;
    temp_y_max = temp_y_max + 2;
end

a = 1;

temp_if_max0_min1 = 0;
%template_path = autoGetFileName_general('template*.tif', output_path,temp_if_max0_min1);
template_path = autoGetFileName_general('maxProjection*.tif', output_path,temp_if_max0_min1);
template = double(loadtiff(template_path));


% fig1 = figure('Name',sprintf('fig1'),'NumberTitle','off');
% set(gcf,'Position',[35+0 35+0 400 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
% 
% scatter(temp_roi_xpix,temp_roi_ypix,60,'s','filled','MarkerFaceColor',[0.4660 0.6740 0.1880]);
% 
% xlim([temp_x_min,temp_x_max]);
% ylim([temp_y_min,temp_y_max]);
% 
% axis equal;
% set(gca, 'YDir','reverse');
% title(sprintf('npix=%d',temp_roi_npix));

fig2 = figure('Name',sprintf('fig2'),'NumberTitle','off');
set(gcf,'Position',[35+0 35+0 400 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
% image(template(temp_x_min:temp_x_max,temp_y_min:temp_y_max).^0.8);

temp_delta = 25;
temp_y_min = max(temp_y_min - temp_delta,1);
temp_y_max = min(temp_y_max + temp_delta,512);
temp_x_min = max(temp_x_min - temp_delta,1);
temp_x_max = min(temp_x_max + temp_delta,512);


temp_template = template(temp_y_min:temp_y_max,temp_x_min:temp_x_max);
temp_template = temp_template / max(temp_template,[],'all');

temp_template2 = (temp_template.^0.7)*255;

image(temp_template2);
hold on
colormap(gray);

temp_I = false(512,512);
for tempi=1:length(temp_roi_xpix)
    temp_I(temp_roi_ypix(tempi),temp_roi_xpix(tempi)) = true;
end
temp_I2 = temp_I(temp_y_min:temp_y_max,temp_x_min:temp_x_max);
B = bwboundaries(temp_I2,'noholes');
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2)+1,boundary(:,1)+1,'w', 'LineWidth', 0.1);
   hold on
end

axis equal;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca, 'YDir','reverse');
a = 1;

fig3 = figure('Name',sprintf('fig3'),'NumberTitle','off');
set(gcf,'Position',[35+0 35+0 500 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

temp_I = false(512,512);
for tempi=1:length(temp_roi_xpix)
    temp_I(temp_roi_ypix(tempi),temp_roi_xpix(tempi)) = true;
end

a = 1;

se = strel('disk',2);
% se = strel('square',2);
% temp_I2 = imopen(temp_I,se);
temp_I2 = imdilate(temp_I,se);

temp_I3 = temp_I2 * 255;

temp_I4 = temp_I3(temp_y_min:temp_y_max,temp_x_min:temp_x_max);


image(temp_I4);
colormap(gray);
axis equal;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca, 'YDir','reverse');

