% Chuan's 22th script (20251214)
% This script: To plot single-neuron's raw trace.
%% Need to run stepB1 first
%% Initialization
close all


F_dff;
cellIndex_suite2p;

% num_frames = 2000;%2000
% num_frames = 10000;%2000

% num_frames = 4000;%2000
num_frames = 3000;%2000
% num_frames = 20000;%2000


% T_start = 42000;%36000,80000,115000,   ,11000,31000
% T_start = 0+20000*12+3000+0000;%1,11,12,15,16
T_start = 0+20000*6+11800;%1,11,12,15,16

%(0-->290), 2,26,27,32,36(same 27),41,57,89,125(low),147(low)

% cellIndex_suite2p_temptempRaw_multi = [2,11,19,26,27,32,36,41,57,89,125,147];
% cellIndex_suite2p_temptempRaw_multi = [2,11,19,26,27,32,36,41,57,89,125,147,0,1,3,4,5,6,7,8,9,10,12,13,14];
% cellIndex_suite2p_temptempRaw_multi = [2,11,32,14];

% cellIndex_suite2p_temptempRaw_multi = [2,11,14];

% cellIndex_suite2p_temptempRaw_multi = [19,25,92];
cellIndex_suite2p_temptempRaw_multi = [19,25,126];
cellIndex_suite2p_temptempRaw_multi = sort(cellIndex_suite2p_temptempRaw_multi,'descend');

% cellIndex_suite2p_temptempRaw = 2;
% cellIndex_suite2p_temptempRaw = 11;% loc3
% cellIndex_suite2p_temptempRaw = 19;% loc5
% cellIndex_suite2p_temptempRaw = 26;
% cellIndex_suite2p_temptempRaw = 27;
% cellIndex_suite2p_temptempRaw = 32;
% cellIndex_suite2p_temptempRaw = 36;
% cellIndex_suite2p_temptempRaw = 41;
% cellIndex_suite2p_temptempRaw = 57;
% cellIndex_suite2p_temptempRaw = 89;
% cellIndex_suite2p_temptempRaw = 125;
% cellIndex_suite2p_temptempRaw = 147;


% cellIndex_suite2p_temptempRaw = 11;% loc3
% T_start = 36000;%37000

% cellIndex_suite2p_temptempRaw = 19;% loc5
% T_start = 251000;%32000,251000

% cellIndex_suite2p_temptempRaw = 2;%
% T_start = 42000;%44000

% cellIndex_suite2p_temptempRaw = 27;
% T_start = 354300;%269000, 101000

% cellIndex_suite2p_temptempRaw = 89;
% T_start = 18000;%37000

% cellIndex_suite2p_temptempRaw = 125;
% cellIndex_suite2p_temptempRaw = 147;

T_end = T_start+num_frames;% -1 mean all frames



backgounrdColor = [1 1 1];
temp_LineWidth = 1;%1


for tempIndex = 1:length(cellIndex_suite2p_temptempRaw_multi)
    
    cellIndex_suite2p_temptempRaw = cellIndex_suite2p_temptempRaw_multi(tempIndex);
    
    temp_id2 = find(cellIndex_suite2p == cellIndex_suite2p_temptempRaw);
    F_dff_exampleRoi = F_dff(temp_id2,:);
    
    % if T_end == -1
    %     F_dff_exampleRoi_valid = F_dff_exampleRoi(T_start:end);
    % else
    %     F_dff_exampleRoi_valid = F_dff_exampleRoi(T_start:T_end);
    % end
    
    
    markerParse_trialLevel;
    
    trialStart_frameIndex = nan(length(markerParse_trialLevel),1);
    for tempi=1:length(trialStart_frameIndex)
        trialStart_frameIndex(tempi) = markerParse_trialLevel{tempi}.currentTrialTotalMarkers_frameIndex(1);
    end
    
    
    x_min = T_start;
    if T_end == -1
        x_max = length(F_dff_exampleRoi);
    else
        x_max = T_end;
    end
    
    y_min = -0.3;%-0.5
    y_max = 3.7;%10
    
    %% Plot
    fig1 = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[35+0 42+0 381 127*0.5*0.8*0.8*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[35+0 42+37*(tempIndex-1) 381 127*0.5*0.8*0.8*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[35+0 42+37*(tempIndex-1) 381*4 127*0.5*0.8*0.8*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','tight','Padding','loose'); %#ok<*NASGU>
    
    plot(F_dff_exampleRoi,'-','LineWidth',temp_LineWidth,'Color',[0 0.4470 0.7410]);%[0.3010 0.7450 0.9330]
    hold on
    
    % plot(trialStart_frameIndex,'-','LineWidth',temp_LineWidth,'Color',[0.75 0.25 0.25]);%[0.3010 0.7450 0.9330]
    % hold on
    
    temptempBoolIndex = trialStart_frameIndex>=x_min & trialStart_frameIndex<=x_max;
    temptempIndex = find(temptempBoolIndex==true);
    
    if sum(temptempBoolIndex) > 0
        for tempi=1:sum(temptempBoolIndex)
            plot(trialStart_frameIndex(temptempIndex(tempi)).*[1 1],[y_min, y_max*0.9],'-','LineWidth',1,'Color',[0.75 0.25 0.25]);%[0.3010 0.7450 0.9330]
            hold on
            a = 1;
        end
    end
    
    % [y_min,y_max] = bounds(F_dff_exampleRoi(x_min:x_max));
    
    xlim([x_min,x_max]);
    ylim([y_min-(y_max-y_min)*0.04,y_max+(y_max-y_min)*0.04]);
    
    % h = gca;
    % h.XAxis.Visible = 'off';
    axis off ;
    
    set(gca,'linewidth',1.5);
    set(gca,'color',backgounrdColor);
    set(gca,'xticklabel',[])
    set(gca,'xtick',[])
    
    set(gca, 'FontSize', 10)
    set(gca,'box','off');% 取消右、上边框
    % ylabel(sprintf('dF/F'), 'FontSize', 10);
    
end