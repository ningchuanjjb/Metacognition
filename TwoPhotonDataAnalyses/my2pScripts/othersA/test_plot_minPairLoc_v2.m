close all

%% Fig
fig = figure('Name',' ','NumberTitle','off');
set(gcf,'Position',[10 50 1500 450]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% t.Title.String = sprintf('%s, %s, r2 = %.4f',ModelName_v,dataSetName,FittingResults.Rsquared.All);
% t.Title.FontSize = 20;
% t.Title.Interpreter = 'none';

temp_color = [0.0 0.0 0.0;
    0.4 0.4 0.4;
    0.7 0.7 0.7;
    0.9 0.9 0.9];


length_rangeHead = 1;
length_rangeTail = numFrames;

for tempi=2:4
    nexttile
    
    %if tempi == 4
    if tempi == 5
        x2 = minPairLoc_peakLocSeq_eachLength_length123(:,:,2);
        x2 = x2 ./ sum(x2,2);
        x3 = minPairLoc_peakLocSeq_eachLength_length123(:,:,3);
        x3 = x3 ./ sum(x3,2);
        
        a = 1;
        x = (x2+x3)./2;
    else
        x = minPairLoc_peakLocSeq_eachLength_length123(:,:,tempi);
        x = x ./ sum(x,2);
    end
    
    temp_max = max(x);
    %temp_totalMax = max(temp_max);
    temp_totalMax = 1;
    
    temp_plot = x;
    % Plot heat map
    imagesc(temp_plot, [0 temp_totalMax]);
    my_gray = gray;
    my_gray = my_gray(end:-1:1,:);
    colormap(my_gray);
    hold on
    
    
    axis equal
    xlim([0.5 numFrames+0.5]);
    ylim([0.5 numFrames+0.5]);
    
    set(gca, 'FontSize', 20)
    set(gca,'YDir','reverse');% Reverse the direction of y axis
    set(gca,'box','off');% 取消右、上边框
    set(gca,'XTick',[1:1:numFrames]);%设置要显示坐标刻度的范围
    set(gca,'YTick',[1:1:numFrames]);%设置要显示坐标刻度的范围
    
    ylabel('Max tuning location','FontSize',18,'FontWeight','bold');
    xlabel('Worst tuning pair','FontSize',18,'FontWeight','bold');
    if tempi == 4
    else
        title(sprintf('length%d',tempi),'fontsize',20);
    end
    
end
c = colorbar;
c.Label.String = 'proportion';



%% Fig
fig = figure('Name',' ','NumberTitle','off');
set(gcf,'Position',[10 50 500 350]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

temp_p = p234_peakLocSeq_dffSEM;
temp_1 = temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,2);
temp_2 = temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,3);
temp_3 = temp_dffSEM_peakLocSeq_eachLength_length123(tempAnovaBoolIndex,4);

temp1_SEM = std(temp_1)/sqrt(length(temp_1));
temp2_SEM = std(temp_2)/sqrt(length(temp_2));
temp3_SEM = std(temp_3)/sqrt(length(temp_3));

temp_y_min = 0;
temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);

bar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)], ...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
hold on
errorbar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
tempTxt = sprintf('');
if temp_p < 0.001
    tempTxt = sprintf('***');
elseif temp_p < 0.01
    tempTxt = sprintf('**');
elseif temp_p < 0.05
    tempTxt = sprintf('*');
end
text(1,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
    'HorizontalAlignment','center');

ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
set(gca, 'FontSize', 14)
set(gca,'XTickLabel', ["length2"; "length3"; "length4"],'FontSize', 12);%给坐标加标签
set(gca,'box','off');% 取消右、上边框
ylabel('dff SEM', 'FontSize', 15, 'FontWeight', 'bold');
title(sprintf('length 123 roi dff SEM (seq of max tuning location) \n Num=%d, linear regression',sum(tempAnovaBoolIndex)),'fontsize',16);
