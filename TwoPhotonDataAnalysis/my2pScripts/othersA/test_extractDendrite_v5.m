clear
close all
load stat_sessionB

padSize = 15;
temp_id = 67 + 1;
temp_stat = stat_sessionB{temp_id};

I = false(512,512);
temp_index = double([temp_stat.xpix;temp_stat.ypix]);
for tempi=1:length(temp_index)
    I(temp_index(2,tempi),temp_index(1,tempi)) = true; % suite2p order is [y,x]
end

[x_min,x_max] = bounds(temp_stat.xpix);
[y_min,y_max] = bounds(temp_stat.ypix);

x_min = max(1,x_min-padSize);
x_max = min(512,x_max+padSize);
y_min = max(1,y_min-padSize);
y_max = min(512,y_max+padSize);

I0 = I(y_min:y_max,x_min:x_max); % suite2p order is [y,x]

% BW0 = bwmorph(I0,'close');
BW0 = I0;
BW1 = imfill(BW0,'holes');

target_npix = 50;
scale = max(1,target_npix/double(temp_stat.npix));
% scale = target_npix/double(temp_stat.npix);
BW2 = imresize(BW1,scale,'Method','nearest');


% regionprops(BW2, 'BoundingBox')
% props = regionprops('table',BW2,'Centroid',...
%     'MajorAxisLength','MinorAxisLength');
% props = regionprops('table',BW2,'MajorAxisLength','MinorAxisLength');

props = regionprops('table',BW2,'MajorAxisLength','MinorAxisLength','MaxFeretProperties','MinFeretProperties','Solidity');


axisRatio = props.MajorAxisLength./props.MinorAxisLength;
% axisRatio = max(axisRatio);
[~,temp_axisRatioIndex] = max(props.MajorAxisLength);
axisRatio = axisRatio(temp_axisRatioIndex);

feretRatio = props.MaxFeretDiameter./props.MinFeretDiameter;
feretRatio = feretRatio(temp_axisRatioIndex);

Solidity = props.Solidity;
somaRatio = double(temp_stat.npix_soma)/double(temp_stat.npix);

fprintf('axisRatio = %.2f.\n',axisRatio);
fprintf('feretRatio = %.2f.\n',feretRatio);
fprintf('Solidity = %.2f.\n',Solidity);
fprintf('somaRatio = %.2f.\n',somaRatio);

BW3 = bwmorph(BW2,'remove');
BW4 = bwmorph(BW2,'skel',Inf);

BW5 = (imdilate(BW4,ones(3,3))) & BW2;
% BW5 = (imdilate(BW4,ones(3,3)));

npix_full = sum(BW2,'all');
npix_contour = sum(BW3,'all');
npix_skel = sum(BW4,'all');
npix_skel_dilate = sum(BW5,'all');

npix_PWM = (npix_skel+npix_contour)/npix_full;
fprintf('npix_PWM = %.2f.\n',npix_PWM);
npix_PWM2 = (npix_skel_dilate)/npix_full;
fprintf('npix_PWM2 = %.2f.\n',npix_PWM2);

npix_raw = double(temp_stat.npix);
axisRatio;
feretRatio;
npix_PWM;

% move from soma to dendrite: npix_PWM>1
% move from dendrite to soma: npix_PWM<0.8

fig1 = figure('Name','I0','NumberTitle','off');
set(gcf,'Position',[35+500 35+125 600 600]);set(gcf,'color',[0.3 0.3 0.3]);set (gca,'position',[0.01,0.01,0.98,0.98]);
imshow(I0);

% fig2 = figure('Name','close','NumberTitle','off');
% set(gcf,'Position',[35+500 35+125 600 600]);set(gcf,'color',[0.3 0.3 0.3]);set (gca,'position',[0.01,0.01,0.98,0.98]);
% imshow(BW0);

% fig21 = figure('Name','holes','NumberTitle','off');
% set(gcf,'Position',[35+500 35+125 600 600]);set(gcf,'color',[0.3 0.3 0.3]);set (gca,'position',[0.01,0.01,0.98,0.98]);
% imshow(BW1);

fig32 = figure('Name','scale','NumberTitle','off');
set(gcf,'Position',[35+500 35+125 600 600]);set(gcf,'color',[0.3 0.3 0.3]);set (gca,'position',[0.01,0.01,0.98,0.98]);
imshow(BW2);

fig3 = figure('Name','contour','NumberTitle','off');
set(gcf,'Position',[35+500 35+125 600 600]);set(gcf,'color',[0.3 0.3 0.3]);set (gca,'position',[0.01,0.01,0.98,0.98]);
imshow(BW3);

fig4 = figure('Name','skel','NumberTitle','off');
set(gcf,'Position',[35+500 35+125 600 600]);set(gcf,'color',[0.3 0.3 0.3]);set (gca,'position',[0.01,0.01,0.98,0.98]);
imshow(BW4);

fig5 = figure('Name','skel','NumberTitle','off');
set(gcf,'Position',[35+500 35+125 600 600]);set(gcf,'color',[0.3 0.3 0.3]);set (gca,'position',[0.01,0.01,0.98,0.98]);
imshow(BW5);

