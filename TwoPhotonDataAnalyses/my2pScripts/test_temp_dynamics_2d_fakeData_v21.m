% Chuan's 16th script (20251214)
% This script: Time courses of memory strength & meta-memory.
%% Intialization
% clear
close all

if_fakeData = 1;

if_plotA = 0;
if_plotB = 0;
if_plotC = 0;
if_plotC2 = 0;
if_plotC3 = 0;
if_plotC4 = 0;
if_plotC5 = 1;
if_plotD = 0;
if_plotE = 0;

color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;
color_choiceMemoryError = [0.5 0.5 0.5];%[0.3 0.3 0.3],[0.5 0.5 0.5]
% color_choiceMemoryError = [91,155,213]/255;

color_memoryQuality = [252,141,89]/255;
color_meta = [145,191,219]/255;
color_prior = [217,217,162]/255;


temp_markerList = ['o','^','s','x'];
crossTime_interval = [1,35,89,137];
timeStamp_max = 115;%89,121,94,115

timeStamp_memoryUpdate = [35,35-1+36];
timeStamp_metaUpdate = [timeStamp_memoryUpdate(2)+1,timeStamp_memoryUpdate(2)+18];
% timeStamp_metaUpdate = [timeStamp_memoryUpdate(2)+1,timeStamp_memoryUpdate(2)+18];

temp_threshold_memoryPrecision = 0.5;
temp_threshold_meta = 0.5;


%% meta_trialLevel_multiPeriod_3types_collapsed_fake
% 3 types: choice-memory correct, choice-offload, choice-memory error

if if_fakeData == 1
    meta_multiPeriod_choiceMemoryCorrect_collapsed_fake = nan(crossTime_interval(4),1);
    meta_multiPeriod_choiceOffload_collapsed_fake = nan(crossTime_interval(4),1);
    meta_multiPeriod_choiceMemoryError_collapsed_fake = nan(crossTime_interval(4),1); %#ok<*NASGU>
    
    meta_multiPeriod_choiceMemoryCorrect_collapsed_fake(crossTime_interval(1):timeStamp_metaUpdate(1)-1) = 0.6;
    meta_multiPeriod_choiceMemoryCorrect_collapsed_fake(timeStamp_metaUpdate(1):timeStamp_metaUpdate(2)) = linspace(0.6,0.935,timeStamp_metaUpdate(2)-timeStamp_metaUpdate(1)+1);%[0.6,0.9]
    meta_multiPeriod_choiceMemoryCorrect_collapsed_fake(timeStamp_metaUpdate(2)+1:end) = meta_multiPeriod_choiceMemoryCorrect_collapsed_fake(timeStamp_metaUpdate(2));
    
    meta_multiPeriod_choiceOffload_collapsed_fake(crossTime_interval(1):timeStamp_metaUpdate(1)-1) = 0.4;
    meta_multiPeriod_choiceOffload_collapsed_fake(timeStamp_metaUpdate(1):timeStamp_metaUpdate(2)) = linspace(0.4,0.1,timeStamp_metaUpdate(2)-timeStamp_metaUpdate(1)+1);
    meta_multiPeriod_choiceOffload_collapsed_fake(timeStamp_metaUpdate(2)+1:end) = meta_multiPeriod_choiceOffload_collapsed_fake(timeStamp_metaUpdate(2));
    
    % meta_multiPeriod_choiceMemoryError_collapsed_fake = meta_multiPeriod_choiceMemoryCorrect_collapsed_fake;
    meta_multiPeriod_choiceMemoryError_collapsed_fake(crossTime_interval(1):timeStamp_metaUpdate(1)-1) = 0.58;
    meta_multiPeriod_choiceMemoryError_collapsed_fake(timeStamp_metaUpdate(1):timeStamp_metaUpdate(2)) = linspace(0.58,0.825,timeStamp_metaUpdate(2)-timeStamp_metaUpdate(1)+1);
    meta_multiPeriod_choiceMemoryError_collapsed_fake(timeStamp_metaUpdate(2)+1:end) = meta_multiPeriod_choiceMemoryError_collapsed_fake(timeStamp_metaUpdate(2));
    
    meta_trialLevel_multiPeriod_3types_collapsed_fake = ...
        [meta_multiPeriod_choiceMemoryCorrect_collapsed_fake,...
        meta_multiPeriod_choiceOffload_collapsed_fake,...
        meta_multiPeriod_choiceMemoryError_collapsed_fake];
    
elseif if_fakeData == 0
    meta_trialLevel_multiPeriod_3types_collapsed_fake = meta_trialLevel_multiPeriod_3types_length123_collapsed;
end


%% precision_trialLevel_multiPeriod_3types_collapsed_fake

if if_fakeData == 1
    precision_multiPeriod_choiceMemoryCorrect_collapsed_fake = nan(crossTime_interval(4),1);
    precision_multiPeriod_choiceOffload_collapsed_fake = nan(crossTime_interval(4),1);
    precision_multiPeriod_choiceMemoryError_collapsed_fake = nan(crossTime_interval(4),1); %#ok<*NASGU>
    
    precision_multiPeriod_choiceMemoryCorrect_collapsed_fake(crossTime_interval(1):timeStamp_memoryUpdate(1)-1) = 0.1;%0.12
    precision_multiPeriod_choiceMemoryCorrect_collapsed_fake(timeStamp_memoryUpdate(1):timeStamp_memoryUpdate(2)) = linspace(0.1,0.8,timeStamp_memoryUpdate(2)-timeStamp_memoryUpdate(1)+1);
    precision_multiPeriod_choiceMemoryCorrect_collapsed_fake(timeStamp_memoryUpdate(2)+1:end) = precision_multiPeriod_choiceMemoryCorrect_collapsed_fake(timeStamp_memoryUpdate(2));
    
    precision_multiPeriod_choiceOffload_collapsed_fake(crossTime_interval(1):timeStamp_memoryUpdate(1)-1) = 0.1;%0.08
    precision_multiPeriod_choiceOffload_collapsed_fake(timeStamp_memoryUpdate(1):timeStamp_memoryUpdate(2)) = linspace(0.1,0.225,timeStamp_memoryUpdate(2)-timeStamp_memoryUpdate(1)+1);
    precision_multiPeriod_choiceOffload_collapsed_fake(timeStamp_memoryUpdate(2)+1:end) = precision_multiPeriod_choiceOffload_collapsed_fake(timeStamp_memoryUpdate(2));
    
    % precision_multiPeriod_choiceMemoryError_collapsed_fake = precision_multiPeriod_choiceOffload_collapsed_fake;
    precision_multiPeriod_choiceMemoryError_collapsed_fake(crossTime_interval(1):timeStamp_memoryUpdate(1)-1) = 0.1;
    precision_multiPeriod_choiceMemoryError_collapsed_fake(timeStamp_memoryUpdate(1):timeStamp_memoryUpdate(2)) = linspace(0.1,0.30,timeStamp_memoryUpdate(2)-timeStamp_memoryUpdate(1)+1);
    precision_multiPeriod_choiceMemoryError_collapsed_fake(timeStamp_memoryUpdate(2)+1:end) = precision_multiPeriod_choiceMemoryError_collapsed_fake(timeStamp_memoryUpdate(2));
    
    precision_trialLevel_multiPeriod_3types_collapsed_fake = ...
        [precision_multiPeriod_choiceMemoryCorrect_collapsed_fake,...
        precision_multiPeriod_choiceOffload_collapsed_fake,...
        precision_multiPeriod_choiceMemoryError_collapsed_fake];
    
elseif if_fakeData == 0
    precision_trialLevel_multiPeriod_3types_collapsed_fake = precision_trialLevel_multiPeriod_3types_length123_collapsed;
end



%% Plot for 3types (3d), x is time
if if_plotA == 1
    fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
    %set(gcf,'Position',[450 350 343*2 142*2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    %set(gcf,'Position',[450 350 343*2 142*2*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    set(gcf,'Position',[50 560 343*2 142*2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    t = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');
    
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_crossTime_interval = crossTime_interval;
    temp_crossTime_max = timeStamp_max;
    
    
    y = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    z = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    temptemp_markerList = temp_markerList;
    
    x = 1:temp_crossTime_marker(end);
    
    windowSize = 5;%8
    y = smoothdata(y,1,'gaussian',windowSize);
    z = smoothdata(z,1,'gaussian',windowSize);
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:)); %#ok<*ASGLU>
    [temp_zmin,temp_zmax] = bounds(z(:));
    
    %     temp_ymin = 0;%0
    %     temp_ymax = 1;%1
    %     temp_zmin = 0;%0
    %     temp_zmax = 1;%1
    
    temp_xmin_raw = temp_xmin;
    temp_ymin_raw = temp_ymin;
    temp_zmin_raw = temp_zmin;
    temp_xmax_raw = temp_xmax;
    temp_ymax_raw = temp_ymax;
    temp_zmax_raw = temp_zmax;
    
    %temp_xmin = temp_xmin_raw-(temp_xmax_raw-temp_xmin_raw)*0.02;
    %temp_xmax = temp_xmax_raw+(temp_xmax_raw-temp_xmin_raw)*0.02;
    temp_ymin = temp_ymin_raw-(temp_ymax_raw-temp_ymin_raw)*0.12;
    temp_ymax = temp_ymax_raw+(temp_ymax_raw-temp_ymin_raw)*0.12;
    temp_zmin = temp_zmin_raw-(temp_zmax_raw-temp_zmin_raw)*0.12;
    temp_zmax = temp_zmax_raw+(temp_zmax_raw-temp_zmin_raw)*0.12;
    
    
    temp_LineWidth = 0.5;%0.75,0.5
    temp_LineWidth2 = 1.5;%1.5
    temp_Linestyle = ':';
    temp_Linestyle2 = '-';
    
    % 3d plot
    plot3(x, y(:,1),z(:,1), temp_Linestyle2, 'LineWidth', temp_LineWidth2, 'Color', color_choiceMemory);
    hold on
    
    plot3(x, y(:,2),z(:,2), temp_Linestyle2, 'LineWidth', temp_LineWidth2, 'Color', color_choiceOffload);
    hold on
    
    plot3(x, y(:,3),z(:,3), temp_Linestyle2, 'LineWidth', temp_LineWidth2, 'Color', color_choiceMemoryError);
    hold on
    
    %plot3([1 1]*temp_crossTime_marker(end), [y(temp_crossTime_marker(end),2) y(temp_crossTime_marker(end),3)]-0.01,...
    %    [z(temp_crossTime_marker(end),2) z(temp_crossTime_marker(end),3)], temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', [1 1 1]*0.5);
    %hold on
    
    
    % 2d plot in meta-memory panel
    plot3(x, temp_ymin.*ones(size(x)),z(:,1), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemory);
    hold on
    
    plot3(x, temp_ymin.*ones(size(x)),z(:,2), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceOffload);
    hold on
    
    plot3(x, temp_ymin.*ones(size(x)),z(:,3), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryError);
    hold on
    
    % 2d plot in memory panel
    plot3(x, y(:,1),temp_zmin.*ones(size(x)), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemory);
    hold on
    
    plot3(x, y(:,2),temp_zmin.*ones(size(x)), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceOffload);
    hold on
    
    plot3(x, y(:,3),temp_zmin.*ones(size(x)), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryError);
    hold on
    
    
    XYBinLimits = [0 1];
    
    temp_offset = 0.003;%0.02
    temp_offset2 = 0.015;%0.015
    temp_lineWidth = 1;%4,1.5
    
    temp_threshold_memoryPrecision = temp_threshold_memoryPrecision; %#ok<*ASGSL>
    temp_threshold_meta = temp_threshold_meta;
    
    temp_facealpha = 0.65;%0.3,0.9
    temp_facealpha2 = 0.5;
    
    % Final time stamp plane
    temp2 = [temp_ymin temp_ymin temp_ymax temp_ymax];
    temp3 = [temp_zmin temp_zmax temp_zmax temp_zmin];
    temp1 = [1 1 1 1].*(timeStamp_metaUpdate(2)+1);
    h = fill3(temp1,temp2,temp3,[1 1 1]*0.9);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % Final time stamp plane
    temp2 = [temp_ymin temp_ymin temp_ymax temp_ymax];
    temp3 = [temp_zmin temp_zmax temp_zmax temp_zmin];
    temp1 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,[1 1 1]*0.9);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    
    % Prior time zone in meta panel
    temp2 = [temp_ymin temp_ymin temp_ymin temp_ymin];% memory axis
    temp3 = [temp_zmin temp_zmax temp_zmax temp_zmin];% meta axis
    temp1 = [temp_crossTime_marker(1) temp_crossTime_marker(1) timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(1)];% time axis
    h = fill3(temp1,temp2,temp3,color_prior);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    % Prior time zone in memory panel
    temp2 = [temp_ymin temp_ymax temp_ymax temp_ymin];% memory axis
    temp3 = [temp_zmin temp_zmin temp_zmin temp_zmin];% meta axis
    temp1 = [temp_crossTime_marker(1) temp_crossTime_marker(1) timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(1)];% time axis
    h = fill3(temp1,temp2,temp3,color_prior);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    % Memory update time zone in meta panel
    temp2 = [temp_ymin temp_ymin temp_ymin temp_ymin];% memory axis
    temp3 = [temp_zmin temp_zmax temp_zmax temp_zmin];% meta axis
    temp1 = [timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(2)+1 timeStamp_memoryUpdate(2)+1];% time axis
    h = fill3(temp1,temp2,temp3,color_memoryQuality);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    % Memory update time zone in memory panel
    temp2 = [temp_ymin temp_ymax temp_ymax temp_ymin];% memory axis
    temp3 = [temp_zmin temp_zmin temp_zmin temp_zmin];% meta axis
    temp1 = [timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(2)+1 timeStamp_memoryUpdate(2)+1];% time axis
    h = fill3(temp1,temp2,temp3,color_memoryQuality);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    % Meta update time zone in meta panel
    temp2 = [temp_ymin temp_ymin temp_ymin temp_ymin];% memory axis
    temp3 = [temp_zmin temp_zmax temp_zmax temp_zmin];% meta axis
    temp1 = [timeStamp_metaUpdate(1) timeStamp_metaUpdate(1) timeStamp_metaUpdate(2)+1 timeStamp_metaUpdate(2)+1];% time axis
    h = fill3(temp1,temp2,temp3,color_meta);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    % Meta update time zone in memory panel
    temp2 = [temp_ymin temp_ymax temp_ymax temp_ymin];% memory axis
    temp3 = [temp_zmin temp_zmin temp_zmin temp_zmin];% meta axis
    temp1 = [timeStamp_metaUpdate(1) timeStamp_metaUpdate(1) timeStamp_metaUpdate(2)+1 timeStamp_metaUpdate(2)+1];% time axis
    h = fill3(temp1,temp2,temp3,color_meta);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    
    % Plot marker
    for temptempi=1:length(temp_crossTime_marker)
        
        %if temptempi==3
        %    continue
        %end
        
        temp_marker = temptemp_markerList(temptempi);
        
        %temptemp_markerTime = temp_crossTime_marker(temptempi);
        
        if temptempi == 1
            temptemp_markerTime = temp_crossTime_marker(1);
            %temp_marker = temptemp_markerList(1);
            %temp_marker = '.';
        elseif temptempi == 2
            temptemp_markerTime = timeStamp_memoryUpdate(1);
            %temp_marker = temptemp_markerList(2);
        elseif temptempi == 3
            temptemp_markerTime = timeStamp_metaUpdate(1);
            %temp_marker = temptemp_markerList(2);
        elseif temptempi == 4
            temptemp_markerTime = temp_crossTime_marker(3);
            %temp_marker = temptemp_markerList(1);
        end
        
        temp_size = 13*1.5*1.3;%20
        if strcmp(temp_marker,'x')
            temp_size = 50*1.5;%20
        end
        temp_alpha = 1;
        temp_LineWidth = 1;%1.5,1.75
        
        %         % 3d plot
        %         scatter3(temptemp_markerTime,y(temptemp_markerTime,1),z(temptemp_markerTime,1),...
        %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'LineWidth',temp_LineWidth,...
        %             'MarkerEdgeAlpha',temp_alpha);
        %         hold on
        %         scatter3(temptemp_markerTime,y(temptemp_markerTime,2),z(temptemp_markerTime,2),...
        %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'LineWidth',temp_LineWidth,...
        %             'MarkerEdgeAlpha',temp_alpha);
        %         hold on
        %         scatter3(temptemp_markerTime,y(temptemp_markerTime,3),z(temptemp_markerTime,3),...
        %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth,...
        %             'MarkerEdgeAlpha',temp_alpha);
        %         hold on
        
        
        temp_LineWidth2 = 1;%1.75
        
        % 2d plot in meta-memory panel
        scatter3(temptemp_markerTime,temp_ymin,z(temptemp_markerTime,1),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temptemp_markerTime,temp_ymin,z(temptemp_markerTime,2),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temptemp_markerTime,temp_ymin,z(temptemp_markerTime,3),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        
        
        % 2d plot in memory panel
        scatter3(temptemp_markerTime,y(temptemp_markerTime,1),temp_zmin,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temptemp_markerTime,y(temptemp_markerTime,2),temp_zmin,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temptemp_markerTime,y(temptemp_markerTime,3),temp_zmin,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        
        
        %         % 2d plot in memory panel
        %         scatter3(temptemp_markerTime,0.975,0.025,...
        %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth,...
        %             'MarkerEdgeAlpha',temp_alpha);
        %         hold on
        
    end
    
    
    
    temp1 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
    yticks([temp1(2) temp1(4)]);
    
    temp2 = temp_zmin:(temp_zmax-temp_zmin)/4:temp_zmax;
    zticks([temp2(2) temp2(4)]);
    
    yticklabels({'Low','High'});
    zticklabels({'Low','High'});
    
    xtickangle(0);
    ytickangle(0);
    ztickangle(0);
    
    %xticks([]);
    xticks(temp_crossTime_marker(1:3));
    xticklabels({'Baseline-on','Sample-on','Delay-on'});
    
    %grid on
    
    set(gca,'ydir','reverse');
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.02 temp_xmax+(temp_xmax-temp_xmin)*0.02]);%0.08
    %     ylim([temp_ymin-(temp_ymax-temp_ymin)*0.02 temp_ymax+(temp_ymax-temp_ymin)*0.02]);
    %     zlim([temp_zmin-(temp_zmax-temp_zmin)*0.02 temp_zmax+(temp_zmax-temp_zmin)*0.02]);
    %     ylim([0 1]);
    %     zlim([0 1]);
    ylim([temp_ymin temp_ymax]);
    zlim([temp_zmin temp_zmax]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)%9
    set(gca,'box','off');% 取消右、上边框
    ylabel('Memory quality', 'FontSize', 9, 'FontWeight', 'normal');
    zlabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
    xlabel('Time', 'FontSize', 9, 'FontWeight', 'normal');
    
    %view([-37.5 30])
    view([-12 20])
    
end


%% Plot for 3types (3d), z is time
if if_plotB == 1
    fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
    set(gcf,'Position',[150 50 343*1.02 142*2*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    t = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');
    
    
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_crossTime_interval = crossTime_interval;
    temp_crossTime_max = timeStamp_max;
    
    
    y = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    x = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    temptemp_markerList = temp_markerList;
    
    z = 1:temp_crossTime_marker(end);
    
    windowSize = 5;%8
    y = smoothdata(y,1,'gaussian',windowSize);
    x = smoothdata(x,1,'gaussian',windowSize);
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:)); %#ok<*ASGLU>
    [temp_zmin,temp_zmax] = bounds(z(:));
    
    %temp_xmin-(temp_xmax-temp_xmin)*0.02
    
    temp_xmin_raw = temp_xmin;
    temp_ymin_raw = temp_ymin;
    temp_zmin_raw = temp_zmin;
    temp_xmax_raw = temp_xmax;
    temp_ymax_raw = temp_ymax;
    temp_zmax_raw = temp_zmax;
    
    temp_xmin = temp_xmin_raw-(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_xmax = temp_xmax_raw+(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_ymin = temp_ymin_raw-(temp_ymax_raw-temp_ymin_raw)*0.12;
    temp_ymax = temp_ymax_raw+(temp_ymax_raw-temp_ymin_raw)*0.12;
    %temp_zmin = temp_zmin_raw-(temp_zmax_raw-temp_zmin_raw)*0.02;
    %temp_zmax = temp_zmax_raw+(temp_zmax_raw-temp_zmin_raw)*0.02;
    
    
    %     temp_ymin = 0;%0
    %     temp_ymax = 1;%1
    %     temp_xmin = 0;%0
    %     temp_xmax = 1;%1
    
    
    temp_LineWidth = 0.9;%0.75,0.5
    temp_LineWidth2 = 1.5;%1.5
    temp_Linestyle = ':';%:
    temp_Linestyle2 = '-';
    
    % 3d plot
    plot3(x(:,1), y(:,1),z, temp_Linestyle2, 'LineWidth', temp_LineWidth2, 'Color', color_choiceMemory);
    hold on
    
    plot3(x(:,2), y(:,2),z, temp_Linestyle2, 'LineWidth', temp_LineWidth2, 'Color', color_choiceOffload);
    hold on
    
    plot3(x(:,3), y(:,3),z, temp_Linestyle2, 'LineWidth', temp_LineWidth2, 'Color', color_choiceMemoryError);
    hold on
    
    %plot3([1 1]*temp_crossTime_marker(end), [y(temp_crossTime_marker(end),2) y(temp_crossTime_marker(end),3)]-0.01,...
    %    [z(temp_crossTime_marker(end),2) z(temp_crossTime_marker(end),3)], temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', [1 1 1]*0.5);
    %hold on
    
    
    % 2d plot in meta-memory panel
    plot3(x(:,1), temp_ymin.*ones(size(z)),z, temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemory);
    hold on
    
    plot3(x(:,2), temp_ymin.*ones(size(z)),z, temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceOffload);
    hold on
    
    plot3(x(:,3), temp_ymin.*ones(size(z)),z, temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryError);
    hold on
    
    % 2d plot in memory panel
    plot3(temp_xmin.*ones(size(z)), y(:,1), z, temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemory);
    hold on
    
    plot3(temp_xmin.*ones(size(z)), y(:,2), z, temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceOffload);
    hold on
    
    plot3(temp_xmin.*ones(size(z)), y(:,3), z, temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryError);
    hold on
    
    
    temp_offset = 0.003;%0.02
    temp_offset2 = 0.015;%0.015
    temp_lineWidth = 1;%4,1.5
    
    temp_threshold_memoryPrecision = temp_threshold_memoryPrecision; %#ok<*ASGSL>
    temp_threshold_meta = temp_threshold_meta;
    
    temp_facealpha = 0.65;%0.3,0.9
    temp_facealpha2 = 0.4;%0.5
    temp_facealpha3 = 0.25;
    
    %     % Memory quality updata time stamp plane
    %     temp2 = [temp_ymin temp_ymin temp_ymax temp_ymax];
    %     temp1 = [temp_xmin temp_xmax temp_xmax temp_xmin];
    %     temp3 = [1 1 1 1].*timeStamp_memoryUpdate(1);
    %     h = fill3(temp1,temp2,temp3,color_memoryQuality);
    %     set(h,'edgealpha',0,'facealpha',temp_facealpha3)
    
    %     % Meta-memory updata time stamp plane
    %     temp2 = [temp_ymin temp_ymin temp_ymax temp_ymax];
    %     temp1 = [temp_xmin temp_xmax temp_xmax temp_xmin];
    %     temp3 = [1 1 1 1].*(timeStamp_metaUpdate(2)+1);
    %     h = fill3(temp1,temp2,temp3,color_meta);
    %     set(h,'edgealpha',0,'facealpha',temp_facealpha3)
    
    % Final time stamp plane
    temp2 = [temp_ymin temp_ymin temp_ymax temp_ymax];
    temp1 = [temp_xmin temp_xmax temp_xmax temp_xmin];
    %temp3 = [1 1 1 1].*temp_crossTime_marker(end);
    temp3 = [1 1 1 1].*(timeStamp_metaUpdate(2)+1);
    h = fill3(temp1,temp2,temp3,[1 1 1]*0.9);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % Final time stamp plane
    temp2 = [temp_ymin temp_ymin temp_ymax temp_ymax];
    temp1 = [temp_xmin temp_xmax temp_xmax temp_xmin];
    temp3 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,[1 1 1]*0.9);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % Prior time zone in meta panel
    temp2 = [temp_ymin temp_ymin temp_ymin temp_ymin];% memory axis
    temp1 = [temp_xmin temp_xmax temp_xmax temp_xmin];% meta axis
    temp3 = [temp_crossTime_marker(1) temp_crossTime_marker(1) timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(1)];% time axis
    h = fill3(temp1,temp2,temp3,color_prior);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    % Prior time zone in memory panel
    temp2 = [temp_ymin temp_ymax temp_ymax temp_ymin];% memory axis
    temp1 = [temp_xmin temp_xmin temp_xmin temp_xmin];% meta axis
    temp3 = [temp_crossTime_marker(1) temp_crossTime_marker(1) timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(1)];% time axis
    h = fill3(temp1,temp2,temp3,color_prior);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    % Memory update time zone in meta panel
    temp2 = [temp_ymin temp_ymin temp_ymin temp_ymin];% memory axis
    temp1 = [temp_xmin temp_xmax temp_xmax temp_xmin];% meta axis
    temp3 = [timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(2)+1 timeStamp_memoryUpdate(2)+1];% time axis
    h = fill3(temp1,temp2,temp3,color_memoryQuality);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    % Memory update time zone in memory panel
    temp2 = [temp_ymin temp_ymax temp_ymax temp_ymin];% memory axis
    temp1 = [temp_xmin temp_xmin temp_xmin temp_xmin];% meta axis
    temp3 = [timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(1) timeStamp_memoryUpdate(2)+1 timeStamp_memoryUpdate(2)+1];% time axis
    h = fill3(temp1,temp2,temp3,color_memoryQuality);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    % Meta update time zone in meta panel
    temp2 = [temp_ymin temp_ymin temp_ymin temp_ymin];% memory axis
    temp1 = [temp_xmin temp_xmax temp_xmax temp_xmin];% meta axis
    temp3 = [timeStamp_metaUpdate(1) timeStamp_metaUpdate(1) timeStamp_metaUpdate(2)+1 timeStamp_metaUpdate(2)+1];% time axis
    h = fill3(temp1,temp2,temp3,color_meta);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    % Meta update time zone in memory panel
    temp2 = [temp_ymin temp_ymax temp_ymax temp_ymin];% memory axis
    temp1 = [temp_xmin temp_xmin temp_xmin temp_xmin];% meta axis
    temp3 = [timeStamp_metaUpdate(1) timeStamp_metaUpdate(1) timeStamp_metaUpdate(2)+1 timeStamp_metaUpdate(2)+1];% time axis
    h = fill3(temp1,temp2,temp3,color_meta);
    set(h,'edgealpha',0,'facealpha',temp_facealpha2)
    
    
    % Plot marker
    for temptempi=1:length(temp_crossTime_marker)
        
        %if temptempi==3
        %    continue
        %end
        
        %temp_marker = temptemp_markerList(temptempi);
        %temp_marker = temptemp_markerList(1);
        temp_marker = '.';
        
        %temptemp_markerTime = temp_crossTime_marker(temptempi);
        
        if temptempi == 1
            temptemp_markerTime = temp_crossTime_marker(1);
            %temp_marker = temptemp_markerList(1);
        elseif temptempi == 2
            temptemp_markerTime = timeStamp_memoryUpdate(1);
            %temp_marker = temptemp_markerList(2);
        elseif temptempi == 3
            temptemp_markerTime = timeStamp_metaUpdate(1);
            %temp_marker = temptemp_markerList(2);
        elseif temptempi == 4
            temptemp_markerTime = temp_crossTime_marker(3);
            %temp_marker = temptemp_markerList(1);
        end
        
        temp_size = 13*1.5*1.3;%20
        if strcmp(temp_marker,'x')
            temp_size = 50*1.5;%20
        end
        temp_alpha = 1;
        temp_LineWidth = 1;%1.5,1.75
        
        %         % 3d plot
        %         scatter3(x(temptemp_markerTime,1),y(temptemp_markerTime,1),temptemp_markerTime,...
        %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'LineWidth',temp_LineWidth,...
        %             'MarkerEdgeAlpha',temp_alpha);
        %         hold on
        %         scatter3(x(temptemp_markerTime,2),y(temptemp_markerTime,2),temptemp_markerTime,...
        %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'LineWidth',temp_LineWidth,...
        %             'MarkerEdgeAlpha',temp_alpha);
        %         hold on
        %         scatter3(x(temptemp_markerTime,3),y(temptemp_markerTime,3),temptemp_markerTime,...
        %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth,...
        %             'MarkerEdgeAlpha',temp_alpha);
        %         hold on
        
        
        temp_LineWidth2 = 1;%1.75
        
        % 2d plot in meta-memory panel
        scatter3(x(temptemp_markerTime,1),temp_ymin,temptemp_markerTime,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(x(temptemp_markerTime,2),temp_ymin,temptemp_markerTime,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(x(temptemp_markerTime,3),temp_ymin,temptemp_markerTime,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        
        
        % 2d plot in memory panel
        scatter3(temp_xmin,y(temptemp_markerTime,1),temptemp_markerTime,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temp_xmin,y(temptemp_markerTime,2),temptemp_markerTime,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temp_xmin,y(temptemp_markerTime,3),temptemp_markerTime,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        
        
        %         % 2d plot in memory panel
        %         scatter3(temptemp_markerTime,0.975,0.025,...
        %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth,...
        %             'MarkerEdgeAlpha',temp_alpha);
        %         hold on
        
    end
    
    
    
    temp1 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
    yticks([temp1(2) temp1(4)]);
    
    temp2 = temp_xmin:(temp_xmax-temp_xmin)/4:temp_xmax;
    xticks([temp2(2) temp2(4)]);
    
    yticklabels({'Low','High'});
    xticklabels({'Low','High'});
    
    xtickangle(0);
    ytickangle(0);
    ztickangle(0);
    
    %zticks([]);
    zticks(temp_crossTime_marker(1:3));
    zticklabels({'Baseline-on','Sample-on','Delay-on'});
    
    %grid on
    
    set(gca,'xdir','reverse');
    set(gca,'ydir','reverse');
    set(gca,'zdir','reverse');
    
    %zlim([temp_zmin temp_zmax]);
    zlim([temp_zmin-(temp_zmax-temp_zmin)*0.02 temp_zmax+(temp_zmax-temp_zmin)*0.02]);%0.08
    %     ylim([temp_ymin-(temp_ymax-temp_ymin)*0.02 temp_ymax+(temp_ymax-temp_ymin)*0.02]);
    %     xlim([temp_xmin-(temp_xmax-temp_xmin)*0.02 temp_xmax+(temp_xmax-temp_xmin)*0.02]);
    %     ylim([0 1]);
    %     xlim([0 1]);
    ylim([temp_ymin temp_ymax]);
    xlim([temp_xmin temp_xmax]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)%9
    set(gca,'box','off');% 取消右、上边框
    ylabel('Memory quality', 'FontSize', 9, 'FontWeight', 'normal');
    xlabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
    zlabel('Time', 'FontSize', 9, 'FontWeight', 'normal');
    
    %view([-37.5 30])
    view([-17 11])
    
end



%% Plot for 3types (2d)
if if_plotC == 1
    fig = figure('Name','asd','NumberTitle','off'); %#ok<*UNRCH>
    set(gcf,'Position',[810 290 300 300]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_crossTime_interval = crossTime_interval;
    temp_crossTime_max = timeStamp_max;
    
    
    x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(2):temp_crossTime_max,:);
    y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(2):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    %temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(2) + 1;
    temp_crossTime_marker(1) = [];
    
    temptemp_markerList = temp_markerList;
    
    windowSize = 50;%8,5
    x(5:end,:) = smoothdata(x(5:end,:),1,'gaussian',windowSize);
    y(5:end,:) = smoothdata(y(5:end,:),1,'gaussian',windowSize);
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:)); %#ok<*ASGLU>
    
    temp_xmin_raw = temp_xmin;
    temp_ymin_raw = temp_ymin;
    temp_xmax_raw = temp_xmax;
    temp_ymax_raw = temp_ymax;
    
    temp_xmin = temp_xmin_raw-(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_xmax = temp_xmax_raw+(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_ymin = temp_ymin_raw-(temp_ymax_raw-temp_ymin_raw)*0.12;
    temp_ymax = temp_ymax_raw+(temp_ymax_raw-temp_ymin_raw)*0.12;
    
    
    
    temp_LineWidth = 2.5;%0.75,1
    plot(x(:,1), y(:,1), '-', 'LineWidth', temp_LineWidth, 'Color', color_choiceMemory);
    hold on
    
    plot(x(:,3), y(:,3), '-', 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryError);
    hold on
    
    plot(x(:,2), y(:,2), '-', 'LineWidth', temp_LineWidth, 'Color', color_choiceOffload);
    hold on
    
    
    
    % Plot marker
    for temptempi=1:length(temp_crossTime_marker)
        
        temp_marker = '.';
        
        if temptempi == 1
            temptemp_markerTime = temp_crossTime_marker(1);
            %temp_marker = temptemp_markerList(1);
        elseif temptempi == 2
            temptemp_markerTime = timeStamp_memoryUpdate(1);
            %temp_marker = temptemp_markerList(2);
        elseif temptempi == 3
            temptemp_markerTime = timeStamp_metaUpdate(1);
            %temp_marker = temptemp_markerList(2);
        elseif temptempi == 4
            temptemp_markerTime = temp_crossTime_marker(3);
            %temp_marker = temptemp_markerList(1);
        end
        
        temp_size = 300;%100
        temp_alpha = 1;
        temp_LineWidth2 = 1;%1.75
        
        % 2d plot in meta-memory panel
        scatter(x(temptemp_markerTime,1),y(temptemp_markerTime,1),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter(x(temptemp_markerTime,2),y(temptemp_markerTime,2),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter(x(temptemp_markerTime,3),y(temptemp_markerTime,3),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
    end
    
    
    le = legend('Memory-correct','Memory-error','Offload',...
        'Location','northwest','fontsize',6.5);%10
    le.ItemTokenSize = ones(1,3)*10;
    legend('boxoff');
    
    
    temp1 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
    yticks([temp1(2) temp1(4)]);
    
    temp2 = temp_xmin:(temp_xmax-temp_xmin)/4:temp_xmax;
    xticks([temp2(2) temp2(4)]);
    
    yticklabels({'Low','High'});
    xticklabels({'Low','High'});
    
    
    ylim([temp_ymin temp_ymax]);
    xlim([temp_xmin temp_xmax]);
    
    
    set(gca,'linewidth',1.5);
    set(gca, 'FontSize', 8);%10
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM quality'), 'FontSize', 9, 'FontWeight', 'normal');%10
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');%10
    
    
end



%% Plot for 3types (2d)
if if_plotC2 == 1
    fig = figure('Name','asd','NumberTitle','off'); %#ok<*UNRCH>
    %set(gcf,'Position',[810 290 300 300]);
    set(gcf,'Position',[510 290 300*0.84 300*0.84]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_crossTime_interval = crossTime_interval;
    temp_crossTime_max = timeStamp_max;
    
    
    x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(2):temp_crossTime_max,:);
    y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(2):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    %temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(2) + 1;
    temp_crossTime_marker(1) = [];
    
    temptemp_markerList = temp_markerList;
    
    windowSize = 50;%8,5
    x(5:end,:) = smoothdata(x(5:end,:),1,'gaussian',windowSize);
    y(5:end,:) = smoothdata(y(5:end,:),1,'gaussian',windowSize);
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:)); %#ok<*ASGLU>
    
    temp_xmin_raw = temp_xmin;
    temp_ymin_raw = temp_ymin;
    temp_xmax_raw = temp_xmax;
    temp_ymax_raw = temp_ymax;
    
    temp_xmin = temp_xmin_raw-(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_xmax = temp_xmax_raw+(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_ymin = temp_ymin_raw-(temp_ymax_raw-temp_ymin_raw)*0.12;
    temp_ymax = temp_ymax_raw+(temp_ymax_raw-temp_ymin_raw)*0.12;
    
    
    plot([temp_xmin temp_xmax],[1 1].*(temp_ymin+temp_ymax)./2, '--', 'LineWidth', 1, 'Color', [0 0 0]);
    hold on
    plot([1 1].*(temp_xmin+temp_xmax)./2,[temp_ymin temp_ymax], '--', 'LineWidth', 1, 'Color', [0 0 0]);
    hold on
    
    temp1 = [temp_xmin (temp_xmin+temp_xmax)./2 (temp_xmin+temp_xmax)./2 temp_xmin];
    temp2 = [temp_ymin temp_ymin (temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2];
    h = fill(temp1,temp2,[1 1 1]*0.75);
    set(h,'edgealpha',0,'facealpha',0.4)
    
    temp1 = [(temp_xmin+temp_xmax)./2 temp_xmax temp_xmax (temp_xmin+temp_xmax)./2];
    temp2 = [(temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2 temp_ymax temp_ymax];
    h = fill(temp1,temp2,[1 1 1]*0.75);
    set(h,'edgealpha',0,'facealpha',0.4)
    
    
    
    h_line = [];
    
    temp_LineWidth = 2.5;%0.75,1
    % choiceMemory
    %temptemp_x = x(temp_crossTime_marker(1),1)+(rand(1,20)-0.5)./10;
    %temptemp_y = y(temp_crossTime_marker(1),1)+(rand(1,20)-0.5)./10;
    temptemp_x = linspace(0.1,x(temp_crossTime_marker(1),1),20);
    temptemp_y = linspace(0.5,y(temp_crossTime_marker(1),1),20);
    h_line = [h_line plot(temptemp_x,temptemp_y, '-', 'LineWidth', temp_LineWidth, 'Color', color_prior)];
    hold on
    h_line = [h_line plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality)];
    hold on
    h_line = [h_line plot(x(timeStamp_memoryUpdate(1):timeStamp_metaUpdate(1),1),...
        y(timeStamp_memoryUpdate(1):timeStamp_metaUpdate(1),1), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta)];
    hold on
    
    % choiceMemoryError
    %temptemp_x = x(temp_crossTime_marker(1),3)+(rand(1,20)-0.5)./10;
    %temptemp_y = y(temp_crossTime_marker(1),3)+(rand(1,20)-0.5)./10;
    temptemp_x = linspace(0.1,x(temp_crossTime_marker(1),3),20);
    temptemp_y = linspace(0.5,y(temp_crossTime_marker(1),3),20);
    plot(temptemp_x,temptemp_y, '-', 'LineWidth', temp_LineWidth, 'Color', color_prior);
    hold on
    plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality);
    hold on
    plot(x(timeStamp_memoryUpdate(1):timeStamp_metaUpdate(1),3),...
        y(timeStamp_memoryUpdate(1):timeStamp_metaUpdate(1),3), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta);
    hold on
    
    % choiceOffload
    %temptemp_x = x(temp_crossTime_marker(1),2)+(rand(1,20)-0.5)./10;
    %temptemp_y = y(temp_crossTime_marker(1),2)+(rand(1,20)-0.5)./10;
    temptemp_x = linspace(0.1,x(temp_crossTime_marker(1),2),20);
    temptemp_y = linspace(0.5,y(temp_crossTime_marker(1),2),20);
    plot(temptemp_x,temptemp_y, '-', 'LineWidth', temp_LineWidth, 'Color', color_prior);
    hold on
    plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality);
    hold on
    plot(x(timeStamp_memoryUpdate(1):timeStamp_metaUpdate(1),2),...
        y(timeStamp_memoryUpdate(1):timeStamp_metaUpdate(1),2), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta);
    hold on
    
    
    h_line = [h_line scatter(0.1,0.5,...
        25,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',1,...
        'MarkerEdgeAlpha',1)];
    hold on
    
    
    % Plot marker
    if false
        %for temptempi=1:length(temp_crossTime_marker)
        if true
            temptempi = 3;
            
            temp_marker = 'o';%'.'
            
            if temptempi == 1
                temptemp_markerTime = temp_crossTime_marker(1);
                %temp_marker = temptemp_markerList(1);
            elseif temptempi == 2
                temptemp_markerTime = timeStamp_memoryUpdate(1);
                %temp_marker = temptemp_markerList(2);
            elseif temptempi == 3
                temptemp_markerTime = timeStamp_metaUpdate(1);
                %temp_marker = temptemp_markerList(2);
            elseif temptempi == 4
                temptemp_markerTime = temp_crossTime_marker(3);
                %temp_marker = temptemp_markerList(1);
            end
            
            temp_size = 25;%100,300
            temp_alpha = 1;
            temp_LineWidth2 = 1;%1.75
            
            % 2d plot in meta-memory panel
            h_line = [h_line scatter(x(temptemp_markerTime,1),y(temptemp_markerTime,1),...
                temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'MarkerFaceColor',color_choiceMemory,'LineWidth',temp_LineWidth2,...
                'MarkerEdgeAlpha',temp_alpha)];
            %         scatter(x(temptemp_markerTime,1),y(temptemp_markerTime,1),...
            %            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            %            'MarkerEdgeAlpha',temp_alpha);
            hold on
            
            h_line = [h_line scatter(x(temptemp_markerTime,3),y(temptemp_markerTime,3),...
                temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'MarkerFaceColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
                'MarkerEdgeAlpha',temp_alpha)];
            %         scatter(x(temptemp_markerTime,3),y(temptemp_markerTime,3),...
            %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            %             'MarkerEdgeAlpha',temp_alpha);
            hold on
            
            h_line = [h_line scatter(x(temptemp_markerTime,2),y(temptemp_markerTime,2),...
                temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'MarkerFaceColor',color_choiceOffload,'LineWidth',temp_LineWidth2,...
                'MarkerEdgeAlpha',temp_alpha)];
            %         scatter(x(temptemp_markerTime,2),y(temptemp_markerTime,2),...
            %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            %             'MarkerEdgeAlpha',temp_alpha);
            hold on
            
        end
    end
    
    
    %le = legend('Memory-correct','Memory-error','Offload','Location','northwest','fontsize',6.5);%10
    le = legend(h_line,'Prior update','WM quality update','Meta-WM update','Trial start','Location','southeast','fontsize',6.5);%10
    %le = legend(h_line,'Prior','WM quality update','Meta-WM update',...
    %    'Trial start','Memory-correct','Memory-error','Offload',...
    %    'Location','southeast','fontsize',6.5);%10
    le.ItemTokenSize = ones(1,3)*10;
    legend('boxoff');
    
    %icons = findobj(le,'type','line');
    %set(icons(end-3:end),'MarkerSize',11);
    %set(icons(end),'MarkerSize',11);
    
    temp1 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
    yticks([temp1(2) temp1(4)]);
    
    temp2 = temp_xmin:(temp_xmax-temp_xmin)/4:temp_xmax;
    xticks([temp2(2) temp2(4)]);
    
    yticklabels({'Low','High'});
    xticklabels({'Low','High'});
    
    
    ylim([temp_ymin temp_ymax]);
    xlim([temp_xmin temp_xmax]);
    
    
    set(gca,'linewidth',1.5);
    set(gca, 'FontSize', 8);%10
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM quality'), 'FontSize', 9, 'FontWeight', 'normal');%10
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');%10
    
    
end



%% Plot for 3types (2d)
if if_plotC3 == 1
    fig = figure('Name','asd','NumberTitle','off'); %#ok<*UNRCH>
    %set(gcf,'Position',[810 290 300 300]);
    set(gcf,'Position',[510 290 300*0.84 300*0.84]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_crossTime_interval = crossTime_interval;
    temp_crossTime_max = timeStamp_max;
    
    
    x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    
    temptemp_markerList = temp_markerList;
    
    
    
    temptemp_y1 = linspace(0.5,y(temp_crossTime_marker(1),1),temp_crossTime_interval(2));
    temptemp_y2 = linspace(0.5,y(temp_crossTime_marker(1),2),temp_crossTime_interval(2));
    temptemp_y3 = linspace(0.5,y(temp_crossTime_marker(1),3),temp_crossTime_interval(2));
    
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),1) = temptemp_y1;
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),2) = temptemp_y2;
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),3) = temptemp_y3;
    
    
    windowSize = 50;%8,5
    x(40:end,:) = smoothdata(x(40:end,:),1,'gaussian',windowSize);
    y(40:end,:) = smoothdata(y(40:end,:),1,'gaussian',windowSize);
    
    windowSize2 = 10;
    x(30:40,:) = smoothdata(x(30:40,:),1,'gaussian',windowSize2);
    
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:)); %#ok<*ASGLU>
    
    temp_xmin_raw = temp_xmin;
    temp_ymin_raw = temp_ymin;
    temp_xmax_raw = temp_xmax;
    temp_ymax_raw = temp_ymax;
    
    temp_xmin = temp_xmin_raw-(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_xmax = temp_xmax_raw+(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_ymin = temp_ymin_raw-(temp_ymax_raw-temp_ymin_raw)*0.12;
    temp_ymax = temp_ymax_raw+(temp_ymax_raw-temp_ymin_raw)*0.12;
    
    
    
    %     temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_collapsed_fake;
    %     temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_collapsed_fake;
    %     temp_crossTime_interval = crossTime_interval;
    %     temp_crossTime_max = timeStamp_max;
    %
    %
    %     x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(2):temp_crossTime_max,:);
    %     y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(2):temp_crossTime_max,:);
    %     temp_crossTime_marker = temp_crossTime_interval;
    %     temp_crossTime_marker(end) = temp_crossTime_max;
    %
    %     %temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    %     temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(2) + 1;
    %     temp_crossTime_marker(1) = [];
    %
    %     temptemp_markerList = temp_markerList;
    %
    %     windowSize = 50;%8,5
    %     x(5:end,:) = smoothdata(x(5:end,:),1,'gaussian',windowSize);
    %     y(5:end,:) = smoothdata(y(5:end,:),1,'gaussian',windowSize);
    %
    %     [temp_xmin,temp_xmax] = bounds(x(:));
    %     [temp_ymin,temp_ymax] = bounds(y(:)); %#ok<*ASGLU>
    %
    %     temp_xmin_raw = temp_xmin;
    %     temp_ymin_raw = temp_ymin;
    %     temp_xmax_raw = temp_xmax;
    %     temp_ymax_raw = temp_ymax;
    %
    %     temp_xmin = temp_xmin_raw-(temp_xmax_raw-temp_xmin_raw)*0.12;
    %     temp_xmax = temp_xmax_raw+(temp_xmax_raw-temp_xmin_raw)*0.12;
    %     temp_ymin = temp_ymin_raw-(temp_ymax_raw-temp_ymin_raw)*0.12;
    %     temp_ymax = temp_ymax_raw+(temp_ymax_raw-temp_ymin_raw)*0.12;
    
    
    plot([temp_xmin temp_xmax],[1 1].*(temp_ymin+temp_ymax)./2, '--', 'LineWidth', 1, 'Color', [0 0 0]);
    hold on
    plot([1 1].*(temp_xmin+temp_xmax)./2,[temp_ymin temp_ymax], '--', 'LineWidth', 1, 'Color', [0 0 0]);
    hold on
    
    temp1 = [temp_xmin (temp_xmin+temp_xmax)./2 (temp_xmin+temp_xmax)./2 temp_xmin];
    temp2 = [temp_ymin temp_ymin (temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2];
    h = fill(temp1,temp2,[1 1 1]*0.75);
    set(h,'edgealpha',0,'facealpha',0.4)
    
    temp1 = [(temp_xmin+temp_xmax)./2 temp_xmax temp_xmax (temp_xmin+temp_xmax)./2];
    temp2 = [(temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2 temp_ymax temp_ymax];
    h = fill(temp1,temp2,[1 1 1]*0.75);
    set(h,'edgealpha',0,'facealpha',0.4)
    
    
    
    h_line = [];
    
    temp_LineWidth = 2.5;%0.75,1
    % choiceMemory
    h_line = [h_line plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1), '-', 'LineWidth', temp_LineWidth, 'Color', color_prior)];
    hold on
    h_line = [h_line plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality)];
    hold on
    h_line = [h_line plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta)];
    hold on
    
    
    % choiceMemoryError
    plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3), '-', 'LineWidth', temp_LineWidth, 'Color', color_prior);
    hold on
    plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),3),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),3), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality);
    hold on
    plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,3),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,3), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta);
    hold on
    
    
    % choiceOffload
    plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2), '-', 'LineWidth', temp_LineWidth, 'Color', color_prior);
    hold on
    plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality);
    hold on
    plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta);
    hold on
    
    
    h_line = [h_line scatter(0.1,0.5,...
        25,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',1,...
        'MarkerEdgeAlpha',1)];
    hold on
    
    
    % Plot marker
    if false
        %for temptempi=1:length(temp_crossTime_marker)
        if true
            temptempi = 3;
            
            temp_marker = 'o';%'.'
            
            if temptempi == 1
                temptemp_markerTime = temp_crossTime_marker(1);
                %temp_marker = temptemp_markerList(1);
            elseif temptempi == 2
                temptemp_markerTime = timeStamp_memoryUpdate(1);
                %temp_marker = temptemp_markerList(2);
            elseif temptempi == 3
                temptemp_markerTime = timeStamp_metaUpdate(1);
                %temp_marker = temptemp_markerList(2);
            elseif temptempi == 4
                temptemp_markerTime = temp_crossTime_marker(3);
                %temp_marker = temptemp_markerList(1);
            end
            
            temp_size = 25;%100,300
            temp_alpha = 1;
            temp_LineWidth2 = 1;%1.75
            
            % 2d plot in meta-memory panel
            h_line = [h_line scatter(x(temptemp_markerTime,1),y(temptemp_markerTime,1),...
                temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'MarkerFaceColor',color_choiceMemory,'LineWidth',temp_LineWidth2,...
                'MarkerEdgeAlpha',temp_alpha)];
            %         scatter(x(temptemp_markerTime,1),y(temptemp_markerTime,1),...
            %            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            %            'MarkerEdgeAlpha',temp_alpha);
            hold on
            
            h_line = [h_line scatter(x(temptemp_markerTime,3),y(temptemp_markerTime,3),...
                temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'MarkerFaceColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
                'MarkerEdgeAlpha',temp_alpha)];
            %         scatter(x(temptemp_markerTime,3),y(temptemp_markerTime,3),...
            %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            %             'MarkerEdgeAlpha',temp_alpha);
            hold on
            
            h_line = [h_line scatter(x(temptemp_markerTime,2),y(temptemp_markerTime,2),...
                temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'MarkerFaceColor',color_choiceOffload,'LineWidth',temp_LineWidth2,...
                'MarkerEdgeAlpha',temp_alpha)];
            %         scatter(x(temptemp_markerTime,2),y(temptemp_markerTime,2),...
            %             temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth2,...
            %             'MarkerEdgeAlpha',temp_alpha);
            hold on
            
        end
    end
    
    
    %le = legend('Memory-correct','Memory-error','Offload','Location','northwest','fontsize',6.5);%10
    le = legend(h_line,'Prior update','WM quality update','Meta-WM update','Trial start','Location','southeast','fontsize',6.5);%10
    %le = legend(h_line,'Prior','WM quality update','Meta-WM update',...
    %    'Trial start','Memory-correct','Memory-error','Offload',...
    %    'Location','southeast','fontsize',6.5);%10
    le.ItemTokenSize = ones(1,3)*10;
    legend('boxoff');
    
    %icons = findobj(le,'type','line');
    %set(icons(end-3:end),'MarkerSize',11);
    %set(icons(end),'MarkerSize',11);
    
    temp1 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
    yticks([temp1(2) temp1(4)]);
    
    temp2 = temp_xmin:(temp_xmax-temp_xmin)/4:temp_xmax;
    xticks([temp2(2) temp2(4)]);
    
    yticklabels({'Low','High'});
    xticklabels({'Low','High'});
    
    %set(gca,'TickLength',[0 0]);
    
    ylim([temp_ymin temp_ymax]);
    xlim([temp_xmin temp_xmax]);
    
    
    set(gca,'linewidth',1.5);
    set(gca, 'FontSize', 8);%10
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM quality'), 'FontSize', 9, 'FontWeight', 'normal');%10
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');%10
    
    
end



%% Plot for 3types (2d)
if if_plotC4 == 1
    fig = figure('Name','asd','NumberTitle','off'); %#ok<*UNRCH>
    set(gcf,'Position',[510 290 300*0.84*0.97 300*0.84*0.9*0.97]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_crossTime_interval = crossTime_interval;
    temp_crossTime_max = timeStamp_max;
    
    
    x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    
    temptemp_markerList = temp_markerList;
    
    
    
    temptemp_y1 = linspace(0.5,y(temp_crossTime_marker(1),1),temp_crossTime_interval(2));
    temptemp_y2 = linspace(0.5,y(temp_crossTime_marker(1),2),temp_crossTime_interval(2));
    temptemp_y3 = linspace(0.5,y(temp_crossTime_marker(1),3),temp_crossTime_interval(2));
    
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),1) = temptemp_y1;
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),2) = temptemp_y2;
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),3) = temptemp_y3;
    
    
    windowSize = 50;%8,5
    x(40:end,:) = smoothdata(x(40:end,:),1,'gaussian',windowSize);
    y(40:end,:) = smoothdata(y(40:end,:),1,'gaussian',windowSize);
    
    windowSize2 = 10;
    x(30:40,:) = smoothdata(x(30:40,:),1,'gaussian',windowSize2);
    
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:)); %#ok<*ASGLU>
    
    temp_xmin_raw = temp_xmin;
    temp_ymin_raw = temp_ymin;
    temp_xmax_raw = temp_xmax;
    temp_ymax_raw = temp_ymax;
    
    temp_xmin = temp_xmin_raw-(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_xmax = temp_xmax_raw+(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_ymin = temp_ymin_raw-(temp_ymax_raw-temp_ymin_raw)*0.12;
    temp_ymax = temp_ymax_raw+(temp_ymax_raw-temp_ymin_raw)*0.12;
    
    temp_xmin = 0;
    temp_xmax = 1;
    temp_ymin = 0;
    temp_ymax = 1;
    
    
    plot([temp_xmin temp_xmax],[1 1].*(temp_ymin+temp_ymax)./2, '--', 'LineWidth', 1, 'Color', [0 0 0]);
    hold on
    plot([1 1].*(temp_xmin+temp_xmax)./2,[temp_ymin temp_ymax], '--', 'LineWidth', 1, 'Color', [0 0 0]);
    hold on
    
    temp1 = [temp_xmin (temp_xmin+temp_xmax)./2 (temp_xmin+temp_xmax)./2 temp_xmin];
    temp2 = [temp_ymin temp_ymin (temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2];
    h = fill(temp1,temp2,[1 1 1]*0.75);
    set(h,'edgealpha',0,'facealpha',0.4)
    
    temp1 = [(temp_xmin+temp_xmax)./2 temp_xmax temp_xmax (temp_xmin+temp_xmax)./2];
    temp2 = [(temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2 temp_ymax temp_ymax];
    h = fill(temp1,temp2,[1 1 1]*0.75);
    set(h,'edgealpha',0,'facealpha',0.4)
    
    
    
    h_line = [];
    
    temp_LineWidth = 2.5;%0.75,1
    % choiceMemory
    h_line = [h_line plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1), '-', 'LineWidth', temp_LineWidth, 'Color', color_prior)];
    hold on
    h_line = [h_line plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality)];
    hold on
    h_line = [h_line plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta)];
    hold on
    
    
    % choiceMemoryError
    plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3), '-', 'LineWidth', 1.5, 'Color', color_prior);
    hold on
    plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),3),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),3), '-', 'LineWidth', 1.5, 'Color', color_memoryQuality);
    hold on
    plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,3),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,3), '-', 'LineWidth', 1.5, 'Color', color_meta);
    hold on
    
    
    % choiceOffload
    plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2), '-', 'LineWidth', temp_LineWidth, 'Color', color_prior);
    hold on
    plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality);
    hold on
    plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta);
    hold on
    
    
    h_line = [h_line scatter(0.1,0.5,...
        25,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',1,...
        'MarkerEdgeAlpha',1)];
    hold on
    
    
    
    le = legend(h_line,'Prior update','WM quality update','Meta-WM update','Trial start','Location','southeast','fontsize',6);%10
    le.ItemTokenSize = ones(1,3)*10;
    legend('boxoff');
    
    
    hAxis=gca;
    hAxis.YAxis.FirstCrossoverValue  = hAxis.XLim(2);
    hAxis.YAxis.SecondCrossoverValue = hAxis.XLim(2);
    
    
    temp1 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
    yticks([temp1(2) temp1(4)]);
    
    temp2 = temp_xmin:(temp_xmax-temp_xmin)/4:temp_xmax;
    xticks([temp2(2) temp2(4)]);
    
    yticklabels({'Low','High'});
    xticklabels({'Low','High'});
    
    set(gca,'TickLength',[0 0]);
    
    ylim([temp_ymin temp_ymax]);
    xlim([temp_xmin temp_xmax]);
    
    
    set(gca,'linewidth',1.5);
    set(gca, 'FontSize', 8);%10
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM quality'), 'FontSize', 9, 'FontWeight', 'normal');%10
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');%10
    
    
    %%  Proj
    if true
        
        temp_meta_mean = temp_meta_trialLevel_collapsed(end,:);
        temp_meta_std = [1,1,1]*0.20;
        
        temp_N = 1000000;
        proj1 = normrnd(temp_meta_mean(1),temp_meta_std(1),[1,temp_N])';
        proj2 = normrnd(temp_meta_mean(2),temp_meta_std(2),[1,temp_N])';
        proj3 = normrnd(temp_meta_mean(3),temp_meta_std(3),[1,temp_N])';
        
        proj123 = [proj1;proj2;proj3];
        [proj123_min,proj123_max] = bounds(proj123);
        
        tempColor_1 = color_choiceMemory;
        tempColor_2 = color_choiceOffload;
        tempColor_3 = color_choiceMemoryError;
        
        if false
            
            fig = figure('Name','Proj of Mixed (3) two-dimension','NumberTitle','off');
            %set(gcf,'Position',[210 460 100*1.38*1.07*0.84*1.03*0.89 120.2*0.46*1.32]);
            set(gcf,'Position',[210 460 300*0.84*0.97*0.8*1.064*0.97 120.2*0.46*1.32]);
            
            t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
            
            nexttile
            
            
            x1 = proj1;
            n=100;
            n=2^ceil(log2(n)); % round up n to the next power of 2;
            [pdf1,xmesh1,bandwidth1] = ksdensity(x1','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
            temp_ratio1 = 1;
            y1 = pdf1*temp_ratio1;
            h1 = plot(xmesh1,y1,'LineWidth',1.5,'color',tempColor_1);
            hold on
            
            
            x2 = proj2;
            n=100;
            n=2^ceil(log2(n)); % round up n to the next power of 2;
            [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
            temp_ratio2 = 1;
            y2 = pdf2*temp_ratio2;
            h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',tempColor_2);
            hold on
            
            
            x3 = proj3;
            n=100;
            n=2^ceil(log2(n)); % round up n to the next power of 2;
            [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
            temp_ratio3 = 1;
            y3 = pdf3*temp_ratio3;
            h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',tempColor_3); %#ok<*NASGU>
            hold on
            
            
            h = fill(xmesh1,y1,tempColor_1);
            set(h,'edgealpha',0,'facealpha',0.25)
            
            h = fill(xmesh2,y2,tempColor_2);
            set(h,'edgealpha',0,'facealpha',0.25)
            
            h = fill(xmesh3,y3,tempColor_3);
            set(h,'edgealpha',0,'facealpha',0.25)
            
            
            [x1_min,x1_max] = bounds(xmesh1);
            [x2_min,x2_max] = bounds(xmesh2);
            [x3_min,x3_max] = bounds(xmesh3);
            x_min = min([x1_min,x2_min,x3_min]);
            x_max = max([x1_max,x2_max,x3_max]);
            
            [y1_min,y1_max] = bounds(y1);
            [y2_min,y2_max] = bounds(y2);
            [y3_min,y3_max] = bounds(y3);
            y_min = min([y1_min,y2_min,y3_min]);
            y_max = max([y1_max,y2_max,y3_max]);
            
            
            %yticks([0 1 2]);
            
            axis off
            
            set(gca,'xdir','reverse')
            
            set(gca,'linewidth',1.5)
            %xlim([x_min+(x_max-x_min)*0.01 x_max-(x_max-x_min)*0.01]);
            xlim([0 1])
            ylim([y_min y_max+(y_max-y_min)*0.4]);
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            %xlabel('Proj', 'FontSize', 9);
            ylabel('Pdf', 'FontSize', 9);
            
        end
        
    end
    
    
    
    %%  Proj
    if true
        
        temp_precision_mean = temp_memoryPrecision_trialLevel_collapsed(end,:);
        temp_precision_std = [1,1,1]*0.20;
        
        temp_N = 1000000;
        proj1 = normrnd(temp_precision_mean(1),temp_precision_std(1),[1,temp_N])';
        proj2 = normrnd(temp_precision_mean(2),temp_precision_std(2),[1,temp_N])';
        proj3 = normrnd(temp_precision_mean(3),temp_precision_std(3),[1,temp_N])';
        
        proj123 = [proj1;proj2;proj3];
        [proj123_min,proj123_max] = bounds(proj123);
        
        tempColor_1 = color_choiceMemory;
        tempColor_2 = color_choiceOffload;
        tempColor_3 = color_choiceMemoryError;
        
        tempColor_12 = [0.3010 0.7450 0.9330]+[0 -0 -0];%[0.3010 0.7450 0.9330]
        
        if false
            
            fig = figure('Name','Proj of Mixed (3) two-dimension','NumberTitle','off');
            set(gcf,'Position',[210 260 300*0.84*0.97*0.8*1.064 120.2*0.46*1.32]);
            
            t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
            
            nexttile
            
            
            x1 = proj1;
            n=100;
            n=2^ceil(log2(n)); % round up n to the next power of 2;
            [pdf1,xmesh1,bandwidth1] = ksdensity(x1','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
            temp_ratio1 = 1;
            y1 = pdf1*temp_ratio1;
            h1 = plot(xmesh1,y1,'LineWidth',1.5,'color',tempColor_1);
            hold on
            
            
            x2 = proj2;
            n=100;
            n=2^ceil(log2(n)); % round up n to the next power of 2;
            [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
            temp_ratio2 = 1;
            y2 = pdf2*temp_ratio2;
            h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',tempColor_2);
            hold on
            
            
            x3 = proj3;
            n=100;
            n=2^ceil(log2(n)); % round up n to the next power of 2;
            [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
            temp_ratio3 = 1;
            y3 = pdf3*temp_ratio3;
            h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',tempColor_3); %#ok<*NASGU>
            hold on
            
            h = fill(xmesh1,y1,tempColor_1);
            set(h,'edgealpha',0,'facealpha',0.25)
            
            h = fill(xmesh2,y2,tempColor_2);
            set(h,'edgealpha',0,'facealpha',0.25)
            
            h = fill(xmesh3,y3,tempColor_3);
            set(h,'edgealpha',0,'facealpha',0.25)
            
            
            [x1_min,x1_max] = bounds(xmesh1);
            [x2_min,x2_max] = bounds(xmesh2);
            [x3_min,x3_max] = bounds(xmesh3);
            x_min = min([x1_min,x2_min,x3_min]);
            x_max = max([x1_max,x2_max,x3_max]);
            
            [y1_min,y1_max] = bounds(y1);
            [y2_min,y2_max] = bounds(y2);
            [y3_min,y3_max] = bounds(y3);
            y_min = min([y1_min,y2_min,y3_min]);
            y_max = max([y1_max,y2_max,y3_max]);
            
            
            %yticks([0 1 2]);
            
            axis off
            
            set(gca,'linewidth',1.5)
            %xlim([x_min+(x_max-x_min)*0.01 x_max-(x_max-x_min)*0.01]);
            xlim([0 1])
            ylim([y_min y_max+(y_max-y_min)*0.4]);
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            %xlabel('Proj', 'FontSize', 9);
            ylabel('Pdf', 'FontSize', 9);
            
        end
        
    end
    
end



%% Plot for 3types (2d)
if if_plotC5 == 1
    fig = figure('Name','asd','NumberTitle','off'); %#ok<*UNRCH>
    set(gcf,'Position',[510 290 300*0.84*0.97 300*0.84*0.9*0.97]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_crossTime_interval = crossTime_interval;
    temp_crossTime_max = timeStamp_max;
    
    
    x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    
    temptemp_markerList = temp_markerList;
    
    
    
    temptemp_y1 = linspace(0.5,y(temp_crossTime_marker(1),1),temp_crossTime_interval(2));
    temptemp_y2 = linspace(0.5,y(temp_crossTime_marker(1),2),temp_crossTime_interval(2));
    temptemp_y3 = linspace(0.5,y(temp_crossTime_marker(1),3),temp_crossTime_interval(2));
    
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),1) = temptemp_y1;
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),2) = temptemp_y2;
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),3) = temptemp_y3;
    
    
    windowSize = 50;%8,5
    x(40:end,:) = smoothdata(x(40:end,:),1,'gaussian',windowSize);
    %y(40:end,:) = smoothdata(y(40:end,:),1,'gaussian',windowSize);
    
    y(40:end,1) = smoothdata(y(40:end,1),1,'gaussian',windowSize*1.5);
    y(40:end,2) = smoothdata(y(40:end,2),1,'gaussian',windowSize);
    y(40:end,3) = smoothdata(y(40:end,3),1,'gaussian',windowSize);
    
    
    windowSize2 = 10;
    %x(30:40,:) = smoothdata(x(30:40,:),1,'gaussian',windowSize2);
    
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:)); %#ok<*ASGLU>
    
    temp_xmin_raw = temp_xmin;
    temp_ymin_raw = temp_ymin;
    temp_xmax_raw = temp_xmax;
    temp_ymax_raw = temp_ymax;
    
    temp_xmin = temp_xmin_raw-(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_xmax = temp_xmax_raw+(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_ymin = temp_ymin_raw-(temp_ymax_raw-temp_ymin_raw)*0.12;
    temp_ymax = temp_ymax_raw+(temp_ymax_raw-temp_ymin_raw)*0.12;
    
    temp_xmin = 0;
    temp_xmax = 1;
    temp_ymin = 0;
    temp_ymax = 1;
    
    
    plot([temp_xmin temp_xmax],[1 1].*(temp_ymin+temp_ymax)./2, ':', 'LineWidth', 0.5, 'Color', [1 1 1]*0);
    hold on
    plot([1 1].*(temp_xmin+temp_xmax)./2,[temp_ymin temp_ymax], ':', 'LineWidth', 0.5, 'Color', [1 1 1]*0);
    hold on
    
    temp1 = [temp_xmin (temp_xmin+temp_xmax)./2 (temp_xmin+temp_xmax)./2 temp_xmin];
    temp2 = [temp_ymin temp_ymin (temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2];
    h = fill(temp1,temp2,[1 1 1]*0.75);
    set(h,'edgealpha',0,'facealpha',0.4)
    
    temp1 = [(temp_xmin+temp_xmax)./2 temp_xmax temp_xmax (temp_xmin+temp_xmax)./2];
    temp2 = [(temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2 temp_ymax temp_ymax];
    h = fill(temp1,temp2,[1 1 1]*0.75);
    set(h,'edgealpha',0,'facealpha',0.4)
    
    
    
    h_line = [];
    
    temp_LineWidth = 2;%0.75,1,2.5,1
    % choiceMemory
    %     h_line = [h_line plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1),...
    %         y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1), '-', 'LineWidth', temp_LineWidth, 'Color', color_choiceMemory)];
    %     hold on
    h_line = [h_line plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1), ':', 'LineWidth', temp_LineWidth, 'Color', [color_choiceMemory,1])];
    hold on
    h_line = [h_line plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1), ':', 'LineWidth', temp_LineWidth, 'Color', [color_choiceMemory,1])];
    hold on
    
    
%     % choiceMemoryError
%     %     plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3),...
%     %         y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3), '-', 'LineWidth', 1.5, 'Color', color_choiceMemoryError);
%     %     hold on
%     plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),3),...
%         y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),3), '-', 'LineWidth', temp_LineWidth*0.75, 'Color', color_choiceMemoryError);
%     hold on
%     plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,3),...
%         y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,3), '-', 'LineWidth', temp_LineWidth*0.75, 'Color', color_choiceMemoryError);
%     hold on
    
    
    % choiceOffload
    %     plot(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2),...
    %         y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2), '-', 'LineWidth', temp_LineWidth, 'Color', color_choiceOffload);
    %     hold on
    plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2), ':', 'LineWidth', temp_LineWidth, 'Color', [color_choiceOffload,1]);
    hold on
    plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2), ':', 'LineWidth', temp_LineWidth, 'Color', [color_choiceOffload,1]);
    hold on
    
    
    %     h_line = [h_line scatter(0.1,0.5,...
    %         25,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',1,...
    %         'MarkerEdgeAlpha',1)];
    %     hold on
    
    
    
    %le = legend(h_line,'Prior update','WM quality update','Meta-WM update','Trial start','Location','southeast','fontsize',6);%10
    %le.ItemTokenSize = ones(1,3)*10;
    %legend('boxoff');
    
    
    
%     % choiceMemory
%     h_line = [h_line plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1)+0.015,...
%         y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1)-0.18, ':', 'LineWidth', temp_LineWidth, 'Color', color_choiceMemory)];
%     hold on
%     h_line = [h_line plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1)+0.015,...
%         y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1)-0.18, ':', 'LineWidth', temp_LineWidth, 'Color', color_choiceMemory)];
%     hold on    
%     
%     % choiceOffload
%     plot(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2)+0.03,...
%         y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2)+0.17, ':', 'LineWidth', temp_LineWidth, 'Color', color_choiceOffload);
%     hold on
%     plot(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2)+0.03,...
%         y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2)+0.17, ':', 'LineWidth', temp_LineWidth, 'Color', color_choiceOffload);
%     hold on    
    
    hAxis=gca;
    hAxis.YAxis.FirstCrossoverValue  = hAxis.XLim(2);
    hAxis.YAxis.SecondCrossoverValue = hAxis.XLim(2);
    
    
    temp1 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
    yticks([temp1(2) temp1(4)]);
    
    temp2 = temp_xmin:(temp_xmax-temp_xmin)/4:temp_xmax;
    xticks([temp2(2) temp2(4)]);
    
    yticklabels({'Low','High'});
    xticklabels({'Low','High'});
    
    set(gca,'TickLength',[0 0]);
    
    ylim([temp_ymin temp_ymax]);
    xlim([temp_xmin temp_xmax]);
    
    %set(gca, 'visible', 'off')
    
    set(gca,'linewidth',1.5);
    set(gca, 'FontSize', 6.5);%10
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM quality'), 'FontSize', 7, 'FontWeight', 'normal');%10
    ylabel(sprintf('Meta-WM'), 'FontSize', 7, 'FontWeight', 'normal');%10
    
    
end



%% Plot for 3types (3d), prior a distinct axis (z)
if if_plotD == 1
    fig = figure('Name','asd','NumberTitle','off'); %#ok<*UNRCH>
    set(gcf,'Position',[810 290 300*0.84 300*0.84]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_collapsed_fake;
    temp_crossTime_interval = crossTime_interval;
    temp_crossTime_max = timeStamp_max;
    
    
    x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    
    temptemp_markerList = temp_markerList;
    
    
    
    temptemp_y1 = linspace(0.5,y(temp_crossTime_marker(1),1),temp_crossTime_interval(2));
    temptemp_y2 = linspace(0.5,y(temp_crossTime_marker(1),2),temp_crossTime_interval(2));
    temptemp_y3 = linspace(0.5,y(temp_crossTime_marker(1),3),temp_crossTime_interval(2));
    
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),1) = temptemp_y1;
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),2) = temptemp_y2;
    y(temp_crossTime_marker(1):temp_crossTime_marker(2),3) = temptemp_y3;
    
    z = y;% z is prior
    temptemp_size = size(z,1)-temp_crossTime_marker(2)+1;
    z(temp_crossTime_marker(2):end,:) = ones(temptemp_size,1) * z(temp_crossTime_marker(2),:);
    
    
    windowSize = 50;%8,5
    x(40:end,:) = smoothdata(x(40:end,:),1,'gaussian',windowSize);
    y(40:end,:) = smoothdata(y(40:end,:),1,'gaussian',windowSize);
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:)); %#ok<*ASGLU>
    [temp_zmin,temp_zmax] = bounds(z(:));
    
    temp_xmin_raw = temp_xmin;
    temp_ymin_raw = temp_ymin;
    temp_zmin_raw = temp_zmin;
    temp_xmax_raw = temp_xmax;
    temp_ymax_raw = temp_ymax;
    temp_zmax_raw = temp_zmax;
    
    temp_xmin = temp_xmin_raw-(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_xmax = temp_xmax_raw+(temp_xmax_raw-temp_xmin_raw)*0.12;
    temp_ymin = temp_ymin_raw-(temp_ymax_raw-temp_ymin_raw)*0.12;
    temp_ymax = temp_ymax_raw+(temp_ymax_raw-temp_ymin_raw)*0.12;
    temp_zmin = temp_zmin_raw-(temp_zmax_raw-temp_zmin_raw)*0.12;
    temp_zmax = temp_zmax_raw+(temp_zmax_raw-temp_zmin_raw)*0.12;
    
    
    
    plot3([temp_xmin temp_xmax],[1 1].*(temp_ymin+temp_ymax)./2,[1 1].*0.5, '--', 'LineWidth', 1, 'Color', [0 0 0]);
    hold on
    plot3([1 1].*(temp_xmin+temp_xmax)./2,[temp_ymin temp_ymax],[1 1].*0.5, '--', 'LineWidth', 1, 'Color', [0 0 0]);
    hold on
    
    plot3([temp_xmin temp_xmax],[1 1].*(temp_ymin+temp_ymax)./2,[1 1].*0.01, '--', 'LineWidth', 1, 'Color', [0 0 0]);
    hold on
    plot3([1 1].*(temp_xmin+temp_xmax)./2,[temp_ymin temp_ymax],[1 1].*0.01, '--', 'LineWidth', 1, 'Color', [0 0 0]);
    hold on
    
    %
    %temp3 = [1 1 1 1].*0.5;
    temp3 = [1 1 1 1].*0.01;
    
    temp1 = [temp_xmin (temp_xmin+temp_xmax)./2 (temp_xmin+temp_xmax)./2 temp_xmin];
    temp2 = [temp_ymin temp_ymin (temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2];
    h = fill3(temp1,temp2,temp3,[1 1 1]*0.75);
    set(h,'edgealpha',0,'facealpha',0.4)
    
    temp1 = [(temp_xmin+temp_xmax)./2 temp_xmax temp_xmax (temp_xmin+temp_xmax)./2];
    temp2 = [(temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2 temp_ymax temp_ymax];
    h = fill3(temp1,temp2,temp3,[1 1 1]*0.75);
    set(h,'edgealpha',0,'facealpha',0.4)
    
    
    %     %
    %     temp3 = [1 1 1 1].*0.5;
    %
    %     temp1 = [temp_xmin (temp_xmin+temp_xmax)./2 (temp_xmin+temp_xmax)./2 temp_xmin];
    %     temp2 = [temp_ymin temp_ymin (temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2];
    %     h = fill3(temp1,temp2,temp3,[1 1 1]*0.75);
    %     set(h,'edgealpha',0,'facealpha',0.4)
    %
    %     temp1 = [(temp_xmin+temp_xmax)./2 temp_xmax temp_xmax (temp_xmin+temp_xmax)./2];
    %     temp2 = [(temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2 temp_ymax temp_ymax];
    %     h = fill3(temp1,temp2,temp3,[1 1 1]*0.75);
    %     set(h,'edgealpha',0,'facealpha',0.4)
    
    
    
    h_line = [];
    
    temp_LineWidth = 2.5;%0.75,1
    % choiceMemory
    h_line = [h_line plot3(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1),...
        z(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),1), '-', 'LineWidth', temp_LineWidth, 'Color', color_prior)];
    hold on
    h_line = [h_line plot3(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1),...
        z(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),1), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality)];
    hold on
    h_line = [h_line plot3(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1),...
        z(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,1), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta)];
    hold on
    
    % choiceMemoryError
    plot3(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3),...
        z(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),3), '-', 'LineWidth', temp_LineWidth, 'Color', color_prior);
    hold on
    plot3(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),3),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),3),...
        z(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),3), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality);
    hold on
    plot3(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,3),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,3),...
        z(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,3), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta);
    hold on
    
    % choiceOffload
    plot3(x(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2),...
        y(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2),...
        z(temp_crossTime_marker(1):timeStamp_memoryUpdate(1),2), '-', 'LineWidth', temp_LineWidth, 'Color', color_prior);
    hold on
    plot3(x(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2),...
        y(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2),...
        z(temp_crossTime_marker(2):timeStamp_memoryUpdate(2),2), '-', 'LineWidth', temp_LineWidth, 'Color', color_memoryQuality);
    hold on
    plot3(x(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2),...
        y(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2),...
        z(timeStamp_memoryUpdate(2):timeStamp_metaUpdate(2)+10,2), '-', 'LineWidth', temp_LineWidth, 'Color', color_meta);
    hold on
    
    
    h_line = [h_line scatter3(0.1,0.5,0.5,...
        25,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',1,...
        'MarkerEdgeAlpha',1)];
    hold on
    
    
    %le = legend(h_line,'Prior','WM quality update','Meta-WM update','Trial start','Location','southeast','fontsize',6.5);%10
    %le.ItemTokenSize = ones(1,3)*10;
    %legend('boxoff');
    
    
    temp2 = temp_xmin:(temp_xmax-temp_xmin)/4:temp_xmax;
    xticks([temp2(2) temp2(4)]);
    
    temp1 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
    yticks([temp1(2) temp1(4)]);
    zticks([temp1(2) temp1(4)]);
    
    
    xticklabels({'Low','High'});
    yticklabels({'Low','High'});
    zticklabels({'Low','High'});
    
    
    xlim([temp_xmin temp_xmax]);
    ylim([temp_ymin temp_ymax]);
    zlim([temp_ymin temp_ymax]);
    
    %     hAxis=gca;
    %     hAxis.ZAxis.FirstCrossoverValue  = hAxis.XLim(1);
    %     hAxis.ZAxis.SecondCrossoverValue = hAxis.XLim(1);
    %     hAxis.YAxis.FirstCrossoverValue  = hAxis.XLim(1);
    %     hAxis.YAxis.SecondCrossoverValue = hAxis.XLim(1);
    %     hAxis.XAxis.FirstCrossoverValue  = hAxis.YLim(1);
    %     hAxis.XAxis.SecondCrossoverValue = hAxis.YLim(1);
    
    set(gca,'linewidth',1.5);
    set(gca, 'FontSize', 8);%10
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('WM quality'), 'FontSize', 9, 'FontWeight', 'normal');%10
    ylabel(sprintf('Meta-WM'), 'FontSize', 9, 'FontWeight', 'normal');%10
    zlabel(sprintf('Prior'), 'FontSize', 9, 'FontWeight', 'normal');%10
    
    view([-33 21])
    
    
end




if if_plotE == 1
    
    %temp_precision_trialMean_3types_fake = precision_trialLevel_multiPeriod_3types_collapsed_fake(end,:);
    %temp_meta_trialMean_3types_fake = meta_trialLevel_multiPeriod_3types_collapsed_fake(end,:);
    
    temp_precision_trialMean_3types_fake = [0.7,0.25,0.3]; %[0.8,0.225,0.3]
    temp_meta_trialMean_3types_fake = [0.7,0.25,0.65];% [0.9,0.1,0.825]
    
    temp_trialNum_3types = [30,40,10]*1;
    temp_std_3types = [1 1 1].*0.15;%0.1
    
    % memoryPrecision
    temp_mu = temp_precision_trialMean_3types_fake(1);
    temp_std = temp_std_3types(1);
    alpha = ( ((1-temp_mu)/(temp_std^2)) -1/temp_mu )*(temp_mu^2);
    beta = alpha*((1/temp_mu)-1);
    memoryPrecision_trialLevel_CMC_fake = betarnd(alpha,beta,[temp_trialNum_3types(1),1]);
    
    temp_mu = temp_precision_trialMean_3types_fake(2);
    temp_std = temp_std_3types(2);
    alpha = ( ((1-temp_mu)/(temp_std^2)) -1/temp_mu )*(temp_mu^2);
    beta = alpha*((1/temp_mu)-1);
    memoryPrecision_trialLevel_CF_fake = betarnd(alpha,beta,[temp_trialNum_3types(1),1]);
    
    temp_mu = temp_precision_trialMean_3types_fake(3);
    temp_std = temp_std_3types(3);
    alpha = ( ((1-temp_mu)/(temp_std^2)) -1/temp_mu )*(temp_mu^2);
    beta = alpha*((1/temp_mu)-1);
    memoryPrecision_trialLevel_CME_fake = betarnd(alpha,beta,[temp_trialNum_3types(1),1]);
    
    % meta
    temp_mu = temp_meta_trialMean_3types_fake(1);
    temp_std = temp_std_3types(1);
    alpha = ( ((1-temp_mu)/(temp_std^2)) -1/temp_mu )*(temp_mu^2);
    beta = alpha*((1/temp_mu)-1);
    meta_trialLevel_CMC_fake = betarnd(alpha,beta,[temp_trialNum_3types(1),1]);
    
    temp_mu = temp_meta_trialMean_3types_fake(2);
    temp_std = temp_std_3types(2);
    alpha = ( ((1-temp_mu)/(temp_std^2)) -1/temp_mu )*(temp_mu^2);
    beta = alpha*((1/temp_mu)-1);
    meta_trialLevel_CF_fake = betarnd(alpha,beta,[temp_trialNum_3types(1),1]);
    
    temp_mu = temp_meta_trialMean_3types_fake(3);
    temp_std = temp_std_3types(3);
    alpha = ( ((1-temp_mu)/(temp_std^2)) -1/temp_mu )*(temp_mu^2);
    beta = alpha*((1/temp_mu)-1);
    meta_trialLevel_CME_fake = betarnd(alpha,beta,[temp_trialNum_3types(1),1]);
    
    
    proj1 = memoryPrecision_trialLevel_CMC_fake;
    proj2 = memoryPrecision_trialLevel_CF_fake;
    proj3 = memoryPrecision_trialLevel_CME_fake;
    proj123 = [proj1;proj2;proj3];
    [proj123_min,proj123_max] = bounds(proj123);
    
    proj4 = meta_trialLevel_CMC_fake;
    proj5 = meta_trialLevel_CF_fake;
    proj6 = meta_trialLevel_CME_fake;
    proj456 = [proj4;proj5;proj6];
    [proj456_min,proj456_max] = bounds(proj456);
    
    %% Plot 2d WM quality & Meta-WM pdf (beta distri)
    if true
        fig = figure('Name','2d WM quality & Meta-WM pdf (beta distri)','NumberTitle','off');
        %set(gcf,'Position',[810 390 196.4*0.905*0.9*0.89*1.16 227.5*0.905*0.9*0.8*1.02*0.965*1.16*0.96]);
        set(gcf,'Position',[810 390 196.4*0.905*0.9*0.89*1.16 227.5*0.905*0.9*0.8*1.02*0.965*1.16*0.96*0.95]);
        
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        x1 = proj1;
        y1 = proj4;
        
        x2 = proj2;
        y2 = proj5;
        
        x3 = proj3;
        y3 = proj6;
        
        [temp_xmin,temp_xmax] = bounds([x1;x2;x3]);
        [temp_ymin,temp_ymax] = bounds([y1;y2;y3]);
        
        
        plot([temp_xmin temp_xmax],[1 1].*(temp_ymin+temp_ymax)./2, '--', 'LineWidth', 1, 'Color', [0 0 0]);
        hold on
        plot([1 1].*(temp_xmin+temp_xmax)./2,[temp_ymin temp_ymax], '--', 'LineWidth', 1, 'Color', [0 0 0]);
        hold on
        
        temp1 = [temp_xmin (temp_xmin+temp_xmax)./2 (temp_xmin+temp_xmax)./2 temp_xmin];
        temp2 = [temp_ymin temp_ymin (temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2];
        h = fill(temp1,temp2,[1 1 1]*0.75);
        set(h,'edgealpha',0,'facealpha',0.4)
        
        temp1 = [(temp_xmin+temp_xmax)./2 temp_xmax temp_xmax (temp_xmin+temp_xmax)./2];
        temp2 = [(temp_ymin+temp_ymax)./2 (temp_ymin+temp_ymax)./2 temp_ymax temp_ymax];
        h = fill(temp1,temp2,[1 1 1]*0.75);
        set(h,'edgealpha',0,'facealpha',0.4)
        
        
        
        temp_flag_mean0_median1 = 1;
        temp_flag_std0_sem1 = 0;
        
        if temp_flag_mean0_median1 == 0
            x1_mu = mean(x1);%
            y1_mu = mean(y1);
            x2_mu = mean(x2);
            y2_mu = mean(y2);
            x3_mu = mean(x3);
            y3_mu = mean(y3);
        elseif temp_flag_mean0_median1 == 1
            x1_mu = median(x1);
            y1_mu = median(y1);
            x2_mu = median(x2);
            y2_mu = median(y2);
            x3_mu = median(x3);
            y3_mu = median(y3);
        end
        
        if temp_flag_std0_sem1 == 0
            x1_std = std(x1);
            y1_std = std(y1);
            
            x2_std = std(x2);
            y2_std = std(y2);
            
            x3_std = std(x3);
            y3_std = std(y3);
        elseif temp_flag_std0_sem1 == 1
            x1_std = std(x1)/sqrt(length(x1));
            y1_std = std(y1)/sqrt(length(y1));
            
            x2_std = std(x2)/sqrt(length(x2));
            y2_std = std(y2)/sqrt(length(y2));
            
            x3_std = std(x3)/sqrt(length(x3));
            y3_std = std(y3)/sqrt(length(y3));
        end
        
        
        
        temp_size = 5;%3,5,7,35,20,15
        temp_MarkerFaceAlpha = 0.45; %0.5,0.3,0.25,0.125,0.175,0.625
        temp_MarkerEdgeAlpha = 0; %0.7,0
        
        
        temp_h12 = scatter([x1;x2], [y1;y2], ...
            temp_size, 'filled', 'MarkerFaceColor', tempColor_12-[0 0 0], ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);%0.2
        hold on
        
        temp_h3 = scatter(x3, y3, ...
            temp_size, 'filled', 'MarkerFaceColor', [0 0 0], ... %color_choiceMemoryError
            'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha);
        hold on
        
        
        
        mdl = fitglm([x1_mu,x2_mu],[y1_mu,y2_mu]);
        x_fit = temp_xmin:0.001:temp_xmax;
        y_fit = predict(mdl,x_fit')';
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', tempColor_12); %[0.35 0.35 0.35 0.7]
        hold on
        
        temp_y3_diff = y3_mu-predict(mdl,x3_mu);
        y_fit = predict(mdl,x_fit')' + temp_y3_diff;
        %plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0 0 0]); %[0.35 0.35 0.35 0.7]
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', tempColor_3); %[0.35 0.35 0.35 0.7]
        hold on
        
        temp_slope_A = mdl.Coefficients.Estimate(2);
        
        temp_angle_A = atan(temp_slope_A)*180/pi;
        
        temp_angle_B = temp_angle_A + 90;
        
        temp_slope_B = tan(temp_angle_B*pi/180);
        
        
        
        temp_size = 70;%70
        temp_MarkerFaceAlpha = 1;
        temp_MarkerEdgeAlpha = 1;
        
        scatter(x1_mu, y1_mu, ...
            temp_size, 'filled', 'MarkerFaceColor', tempColor_1,'MarkerEdgeColor',[1 1 1], ...
            'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
        hold on
        scatter(x2_mu, y2_mu, ...
            temp_size, 'filled', 'MarkerFaceColor', tempColor_2,'MarkerEdgeColor',[1 1 1], ...
            'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
        hold on
        scatter(x3_mu, y3_mu, ...
            temp_size, 'filled', 'MarkerFaceColor', tempColor_3,'MarkerEdgeColor',[1 1 1], ...
            'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', temp_MarkerEdgeAlpha,'LineWidth',1);
        hold on
        
        
        %xticks([0 1]);
        %yticks([0 1]);
        
        temp1 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
        yticks([temp1(2) temp1(4)]);
        
        temp2 = temp_xmin:(temp_xmax-temp_xmin)/4:temp_xmax;
        xticks([temp2(2) temp2(4)]);
        
        yticklabels({'Low','High'});
        xticklabels({'Low','High'});
        
        set(gca,'TickLength',[0 0]);
        
        %xlim([-0.05 1.05]);
        %ylim([-0.05 1.05]);
        ylim([temp_ymin temp_ymax]);
        xlim([temp_xmin temp_xmax]);
        
        hAxis=gca;
        hAxis.YAxis.FirstCrossoverValue  = hAxis.XLim(2);
        hAxis.YAxis.SecondCrossoverValue = hAxis.XLim(2);
        
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 6.5)
        set(gca,'box','off');% 取消右、上边框
        xlabel(sprintf('WM quality'), 'FontSize', 7, 'FontWeight', 'normal');
        ylabel(sprintf('Meta-WM'), 'FontSize', 7, 'FontWeight', 'normal');
        
        %title(sprintf('asd'), 'FontSize', 9, 'FontWeight', 'normal');
        %subtitle(sprintf('asd'), 'FontSize', 8, 'FontWeight', 'normal');
    end
    
    
    %% Plot WM quality & Meta-WM pdf (beta distri)
    if true
        fig = figure('Name','WM quality & Meta-WM pdf (beta distri)','NumberTitle','off');
        set(gcf,'Position',[210 60 100*1.38*1.07*0.84*1.03*0.89*2*1.6 120.2*0.46*1.32*2*0.8]);
        
        t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        
        x1 = proj1;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf1,xmesh1,bandwidth1] = ksdensity(x1','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
        temp_ratio1 = 1;
        y1 = pdf1*temp_ratio1;
        h1 = plot(xmesh1,y1,'LineWidth',1.5,'color',tempColor_1);
        hold on
        
        
        x2 = proj2;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
        temp_ratio2 = 1;
        y2 = pdf2*temp_ratio2;
        h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',tempColor_2);
        hold on
        
        
        x3 = proj3;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[proj123_min-0.01,proj123_max+0.01],'BoundaryCorrection', 'Reflection');
        temp_ratio3 = 1;
        y3 = pdf3*temp_ratio3;
        h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',tempColor_3); %#ok<*NASGU>
        hold on
        
        
        [x1_min,x1_max] = bounds(xmesh1);
        [x2_min,x2_max] = bounds(xmesh2);
        [x3_min,x3_max] = bounds(xmesh3);
        x_min = min([x1_min,x2_min,x3_min]);
        x_max = max([x1_max,x2_max,x3_max]);
        
        [y1_min,y1_max] = bounds(y1);
        [y2_min,y2_max] = bounds(y2);
        [y3_min,y3_max] = bounds(y3);
        y_min = min([y1_min,y2_min,y3_min]);
        y_max = max([y1_max,y2_max,y3_max]);
        
        
        %yticks([0 1 2]);
        
        %axis off
        
        set(gca,'linewidth',1.5)
        xlim([x_min+(x_max-x_min)*0.01 x_max-(x_max-x_min)*0.01]);
        %xlim([0 1]);
        ylim([y_min y_max+(y_max-y_min)*0.1]);
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('WM quality', 'FontSize', 9);
        ylabel('Pdf', 'FontSize', 9);
        
        
        
        
        nexttile
        
        
        x1 = proj4;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf1,xmesh1,bandwidth1] = ksdensity(x1','NumPoints',n,'Function','pdf','Support',[proj456_min-0.01,proj456_max+0.01],'BoundaryCorrection', 'Reflection');
        temp_ratio1 = 1;
        y1 = pdf1*temp_ratio1;
        h1 = plot(xmesh1,y1,'LineWidth',1.5,'color',tempColor_1);
        hold on
        
        
        x2 = proj5;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[proj456_min-0.01,proj456_max+0.01],'BoundaryCorrection', 'Reflection');
        temp_ratio2 = 1;
        y2 = pdf2*temp_ratio2;
        h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',tempColor_2);
        hold on
        
        
        x3 = proj6;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[proj456_min-0.01,proj456_max+0.01],'BoundaryCorrection', 'Reflection');
        temp_ratio3 = 1;
        y3 = pdf3*temp_ratio3;
        h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',tempColor_3); %#ok<*NASGU>
        hold on
        
        
        [x1_min,x1_max] = bounds(xmesh1);
        [x2_min,x2_max] = bounds(xmesh2);
        [x3_min,x3_max] = bounds(xmesh3);
        x_min = min([x1_min,x2_min,x3_min]);
        x_max = max([x1_max,x2_max,x3_max]);
        
        [y1_min,y1_max] = bounds(y1);
        [y2_min,y2_max] = bounds(y2);
        [y3_min,y3_max] = bounds(y3);
        y_min = min([y1_min,y2_min,y3_min]);
        y_max = max([y1_max,y2_max,y3_max]);
        
        
        %yticks([0 1 2]);
        
        %axis off
        
        set(gca,'linewidth',1.5)
        xlim([x_min+(x_max-x_min)*0.01 x_max-(x_max-x_min)*0.01]);
        %xlim([0 1]);
        ylim([y_min y_max+(y_max-y_min)*0.1]);
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta', 'FontSize', 9);
        ylabel('Pdf', 'FontSize', 9);
        
    end
    
    
    
    
end



%% End
