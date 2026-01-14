close all

color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;


%% Figure 1
fig1 = figure (1);
set(gcf,'Position',[0 50 500 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

x = [0:0.001:1]; %#ok<NBRAK>
y = normpdf(x,0.5,0.1);
% plot(x,y,'Color','k','LineWidth',6);
plot(x,y,'Color',[0 0 0],'LineWidth',6);
axis off
ylim([0 max(y)+0.1*max(y)]);

fileName_fig1 = 'normDistri_v1.pdf';
% to judge whether save figure or not
if false
    exportgraphics(fig1, fileName_fig1, 'BackgroundColor', 'none', 'ContentType', 'vector');
end

%% Figure 2
fig2 = figure (2);
set(gcf,'Position',[0 50 500 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

x1 = [0:0.001:1]; %#ok<NBRAK>
y1 = normpdf(x,0.5,0.1);

% tempSampleNum1 = round(length(x1)/2);

temp_coeff = 0.6;%0.5,0.45,0.55

tempSampleNum1 = round(length(x1)*temp_coeff);
area(x1(1:tempSampleNum1),y1(1:tempSampleNum1),...
    'FaceColor',color_choiceOffload);%[0.75 0.25 0.25]
hold on
area(x1(tempSampleNum1+1:end),y1(tempSampleNum1+1:end),...
    'FaceColor',color_choiceMemory);%[0.25 0.75 0.25]
hold on

plot(x1,y1,'Color','k','LineWidth',6);
hold on

% x2 = [0.5 0.5];
x2 = [1 1]*temp_coeff;
% y2 = [-0.2*max(y1) max(y1)+0.2*max(y1)];
y2 = [0 max(y1)+0.1*max(y1)];
plot(x2,y2,'Color','k','LineWidth',6);
hold on 

axis off
ylim([0 max(y1)+0.1*max(y1)]);

fileName_fig2 = 'normDistri_withThreshold_v1.pdf';
% to judge whether save figure or not
if false
    exportgraphics(fig2, fileName_fig2, 'BackgroundColor', 'none', 'ContentType', 'vector');
end