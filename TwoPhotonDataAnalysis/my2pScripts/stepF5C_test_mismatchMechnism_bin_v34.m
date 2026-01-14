% Chuan's 8th script (20251214)
%% Initialization

if_resample_meta = if_resample_meta; %#ok<*ASGSL>

if_plot = 1;

if_plot_A = 0;
if_plot_B = 1;
if_plot_C = 1;
if_plot_D = 1;

% if_testMeta_0baselineBin_1delay1Bin = 1;
if_testMeta_0baselineBin_1delay1Bin = if_trainMeta_0baseline_1delay1;

% memoryPrecision_trialLevel;
% meta_trialLevel;


color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;
color_choiceMemoryError = [0.3 0.3 0.3];

%% meta_trialLevel_test
if if_testMeta_0baselineBin_1delay1Bin == 0
    meta_trialLevel_test = meta_trialLevel_baseline;
elseif if_testMeta_0baselineBin_1delay1Bin == 1
    meta_trialLevel_test = meta_trialLevel_delay1;
end



%%
% %% Test metalegend
% temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
% F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,temptemp_range),3);
% F_dff_baselineBin = double(F_dff_baselineBin);
% F_dff_baselineBin = F_dff_baselineBin + eps;
% 
% F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
% F_dff_decisionBin1 = double(F_dff_decisionBin1);
% F_dff_decisionBin1 = F_dff_decisionBin1 + eps;
% 
% fun_testMeta_currentTime_Name_v = autoGetFunName_myScripts('fun_testMeta_currentTime', [targetPATH '\functions']);
% fun_testMeta_currentTime = str2func(fun_testMeta_currentTime_Name_v);
% 
% options_testMeta = struct;
% options_testMeta.KFold_num = KFold_num_meta;
% options_testMeta.svm_Meta = svm_Meta;
% % options_testMeta.choiceBoolIndex = choiceBoolIndex;
% options_testMeta.choiceBoolIndex = choiceBoolIndex_validLength;
% options_testMeta.noChoiceBoolIndex = noChoiceBoolIndex_validLength;
% 
% % options_testMeta.seqIndex = seqIndex;
% 
% 
% %%  Compute meta-memory at x period (training at x period)
% if if_testMeta_0baselineBin_1delay1Bin == 0
%     F_dff_currentTime = F_dff_baselineBin;
%     testMeta_output = fun_testMeta_currentTime(F_dff_currentTime,options_testMeta);
%     meta_trialLevel_test = testMeta_output.meta_trialLevel_currentTime;
%     
% elseif if_testMeta_0baselineBin_1delay1Bin == 1
%     F_dff_currentTime = F_dff_decisionBin1;
%     testMeta_output = fun_testMeta_currentTime(F_dff_currentTime,options_testMeta);
%     meta_trialLevel_test = testMeta_output.meta_trialLevel_currentTime;
%     
% end
% % end
% meta_trialLevel_test;



%% Split meta_trialLevel_test
% all trials
% meta_trialLevel_x_choiceMemory = meta_trialLevel_test(choiceMemoryBoolIndex);
% meta_trialLevel_x_choiceOffload = meta_trialLevel_test(choiceOffloadBoolIndex);

% meta_trialLevel_x_choiceMemory_raw = meta_trialLevel_test(choiceMemoryBoolIndex);
% meta_trialLevel_x_choiceOffload_raw = meta_trialLevel_test(choiceOffloadBoolIndex);
% meta_trialLevel_x_choiceMemory = meta_trialLevel_test(choiceMemoryBoolIndex&(~isnan(memoryPrecision_trialLevel)'));
% meta_trialLevel_x_choiceOffload = meta_trialLevel_test(choiceOffloadBoolIndex(~isnan(memoryPrecision_trialLevel)'));
meta_trialLevel_x_choiceMemory = meta_trialLevel_test(choiceMemoryBoolIndex_validLength);
meta_trialLevel_x_choiceOffload = meta_trialLevel_test(choiceOffloadBoolIndex_validLength);

memoryPrecision_trialLevel_choiceMemoryCorrect = memoryPrecision_trialLevel(choiceMemoryBoolIndex_validLength & choiceMemoryCorrectBoolIndex);
memoryPrecision_trialLevel_choiceMemoryError = memoryPrecision_trialLevel(choiceMemoryBoolIndex_validLength & choiceMemoryErrorBoolIndex);


% mismatch trials
trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh;
trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow;
trialBoolIndex_memoryPrecisionLowError_choiceMemory;
trialBoolIndex_memoryPrecisionHigh_choiceOffload;

a1 = sum(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh);
a2 = sum(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow);

if if_plot_fineTuningMismatch == 1
    meta_trialLevel_x_overMismatch = meta_trialLevel_test(trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh);
    meta_trialLevel_x_underMismatch = meta_trialLevel_test(trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow);
elseif if_plot_fineTuningMismatch == 0
    meta_trialLevel_x_overMismatch = meta_trialLevel_test(trialBoolIndex_memoryPrecisionLowError_choiceMemory);
    meta_trialLevel_x_underMismatch = meta_trialLevel_test(trialBoolIndex_memoryPrecisionHigh_choiceOffload);
end

aa1 = mean(meta_trialLevel_x_overMismatch);
aa2 = mean(meta_trialLevel_x_underMismatch);


% match trials
trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh;
trialBoolIndex_memoryPrecisionLow_choiceOffloadLow;

a3 = sum(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh);
a4 = sum(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow);

if if_plot_fineTuningMismatch == 1
    meta_trialLevel_x_highMatch = meta_trialLevel_test(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh);
    meta_trialLevel_x_lowMatch = meta_trialLevel_test(trialBoolIndex_memoryPrecisionLow_choiceOffloadLow);
elseif if_plot_fineTuningMismatch == 0
    meta_trialLevel_x_highMatch = meta_trialLevel_test(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory);
    meta_trialLevel_x_lowMatch = meta_trialLevel_test(trialBoolIndex_memoryPrecisionLow_choiceOffload);
end

aa3 = mean(meta_trialLevel_x_highMatch);
aa4 = mean(meta_trialLevel_x_lowMatch);
aa34A = mean([meta_trialLevel_x_highMatch;meta_trialLevel_x_lowMatch]);
aa34B = (aa3+aa4)/2;

meta_trialLevel_x_match = [meta_trialLevel_x_highMatch;meta_trialLevel_x_lowMatch];

meta_trialLevel_x_mean = mean(meta_trialLevel_test);


%% ttest2
meta_trialLevel_x_choiceMemory;
meta_trialLevel_x_choiceOffload;
[~,p_choiceMemory_choiceOffload] = ttest2(meta_trialLevel_x_choiceMemory,meta_trialLevel_x_choiceOffload);


meta_trialLevel_x_lowMatch;
meta_trialLevel_x_highMatch;
meta_trialLevel_x_match;
[~,p_lowMatch_highMatch] = ttest2(meta_trialLevel_x_lowMatch,meta_trialLevel_x_highMatch);
[~,p_lowMatch_match] = ttest2(meta_trialLevel_x_lowMatch,meta_trialLevel_x_match);
[~,p_highMatch_match] = ttest2(meta_trialLevel_x_highMatch,meta_trialLevel_x_match);


meta_trialLevel_x_underMismatch;
meta_trialLevel_x_overMismatch;
meta_trialLevel_x_match;
[~,p_under_over] = ttest2(meta_trialLevel_x_underMismatch,meta_trialLevel_x_overMismatch);
[~,p_under_match] = ttest2(meta_trialLevel_x_underMismatch,meta_trialLevel_x_match);
[~,p_over_match] = ttest2(meta_trialLevel_x_overMismatch,meta_trialLevel_x_match);

[~,p_under_match_tail] = ttest2(meta_trialLevel_x_underMismatch,meta_trialLevel_x_match,'Tail','left');
[~,p_over_match_tail] = ttest2(meta_trialLevel_x_overMismatch,meta_trialLevel_x_match,'Tail','right');


%% Plot
if if_plot == 1
    close all
    
    if if_plot_A == 1
        %%% Pdf of memory precision in choiceMemory trials and choiceOffload trials
        fig = figure('Name','asd','NumberTitle','off');
        % set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 1600 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 720 240]);
        t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
        
        if if_memeoryPrecision_stimuli0_response1 == 0
            if if_plot_fineTuningMismatch == 0
                temp_str = 'stimuli, roughTuning';
            elseif if_plot_fineTuningMismatch == 1
                temp_str = 'stimuli, fineTuning';
            end
        elseif if_memeoryPrecision_stimuli0_response1 == 1
            if if_plot_fineTuningMismatch == 0
                temp_str = 'response, roughTuning';
            elseif if_plot_fineTuningMismatch == 1
                temp_str = 'response, fineTuning';
            end
        end
        
        if if_trainMeta_0baseline_1delay1 == 1
            if if_testMeta_0baselineBin_1delay1Bin == 0
                t.Title.String = sprintf('Meta-WM of single trial (train delay1, test baseline) \n %s, %s',...
                    FOVName_currentFOV2,temp_str);
            elseif if_testMeta_0baselineBin_1delay1Bin == 1
                t.Title.String = sprintf('Meta-WM of single trial (train delay1, test delay1) \n %s, %s',...
                    FOVName_currentFOV2,temp_str);
            end
        elseif if_trainMeta_0baseline_1delay1 == 0
            if if_testMeta_0baselineBin_1delay1Bin == 0
                t.Title.String = sprintf('Meta-WM of single trial (train baseline, test baseline) \n %s, %s',...
                    FOVName_currentFOV2,temp_str);
            elseif if_testMeta_0baselineBin_1delay1Bin == 1
                t.Title.String = sprintf('Meta-WM of single trial (train baseline, test delay1) \n %s, %s',...
                    FOVName_currentFOV2,temp_str);
            end
        end
        t.Title.FontSize = 12;
        t.Title.Interpreter = 'none';
        
        %% Plot all trials
        nexttile
        
        x2 = meta_trialLevel_x_choiceMemory;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        %[pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
        [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        temp_ratio2 = 1;
        %     temp_ratio2 = length(meta_trialLevel_x_choiceMemory)/...
        %         (length(meta_trialLevel_x_choiceMemory)+length(meta_trialLevel_x_choiceOffload));
        y2 = pdf2*temp_ratio2;
        h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',color_choiceMemory);
        hold on
        
        x3 = meta_trialLevel_x_choiceOffload;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        %[pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
        [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        temp_ratio3 = 1;
        %     temp_ratio3 = length(meta_trialLevel_x_choiceOffload)/...
        %         (length(meta_trialLevel_x_choiceMemory)+length(meta_trialLevel_x_choiceOffload));
        y3 = pdf3*temp_ratio3;
        h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',color_choiceOffload); %#ok<*NASGU>
        hold on
        
        [x2_min,x2_max] = bounds(xmesh2);
        [x3_min,x3_max] = bounds(xmesh3);
        x_min = min([x2_min,x3_min]);
        x_max = max([x2_max,x3_max]);
        
        [y2_min,y2_max] = bounds(y2);
        [y3_min,y3_max] = bounds(y3);
        y_min = min([y2_min,y3_min]);
        y_max = max([y2_max,y3_max]);
        
        %le = legend(sprintf('choiceMemory(%d)',length(x2)),...
        %    sprintf('choiceOffload(%d)',length(x3)),...
        %    'Location','northwest','fontsize',9);
        %le = legend(sprintf('memory'),...
        %   sprintf('offload'),...
        %   'Location','northwest','fontsize',9);
        le = legend(sprintf('Choice-memory'),...
            sprintf('Choice-offload'),...
            'Location','northwest','fontsize',9);
        le.ItemTokenSize = ones(1,2)*10;
        
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        xlim([0 1]);
        ylim([y_min y_max+(y_max-y_min)*0.4]);%0.1-->0.4-->1
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-WM', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Pdf', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Free-choice, p=%.4f',p_choiceMemory_choiceOffload),...
            'FontSize',12,'FontWeight','bold');
        temp_title.Interpreter = 'none';
        
        %% Plot match
        nexttile
        %     x1 = meta_trialLevel_x_match;
        %     n=100;
        %     n=2^ceil(log2(n)); % round up n to the next power of 2;
        %     [pdf1,xmesh1,bandwidth1] = ksdensity(x1','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection'); %#ok<*ASGLU>
        %     %temp_ratio1 = sum(choiceMemoryBoolIndex)/sum(choiceBoolIndex);
        %     temp_ratio1 = 1;
        %     y1 = pdf1*temp_ratio1;
        %     h1 = plot(xmesh1,y1,'LineWidth',1.5,'color',[0 0 0]); %#ok<*NASGU>
        %     hold on
        %     %h1 = histogram(x1,'Normalization','pdf');
        %     %hold on
        
        
        x2 = meta_trialLevel_x_highMatch;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');     %#ok<*ASGLU>
        %temp_ratio2 = sum(choiceOffloadBoolIndex)/sum(choiceBoolIndex);
        temp_ratio2 = 1;
        %     temp_ratio2 = length(meta_trialLevel_x_highMatch)/...
        %         (length(meta_trialLevel_x_highMatch)+length(meta_trialLevel_x_lowMatch));
        y2 = pdf2*temp_ratio2;
        h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',color_choiceMemory);
        hold on
        %h2 = histogram(x2,'Normalization','pdf');
        %hold on
        
        x3 = meta_trialLevel_x_lowMatch;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        %temp_ratio3 = sum(choiceOffloadBoolIndex)/sum(choiceBoolIndex);
        temp_ratio3 = 1;
        %     temp_ratio3 = length(meta_trialLevel_x_lowMatch)/...
        %         (length(meta_trialLevel_x_highMatch)+length(meta_trialLevel_x_lowMatch));
        y3 = pdf3*temp_ratio3;
        h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',color_choiceOffload);
        hold on
        %h3 = histogram(x3,'Normalization','pdf');
        %hold on
        
        
        %     [x1_min,x1_max] = bounds(xmesh1);
        [x2_min,x2_max] = bounds(xmesh2);
        [x3_min,x3_max] = bounds(xmesh3);
        %     x_min = min([x1_min,x2_min,x3_min]);
        %     x_max = max([x1_max,x2_max,x3_max]);
        x_min = min([x2_min,x3_min]);
        x_max = max([x2_max,x3_max]);
        
        
        %     [y1_min,y1_max] = bounds(y1);
        [y2_min,y2_max] = bounds(y2);
        [y3_min,y3_max] = bounds(y3);
        %     y_min = min([y1_min,y2_min,y3_min]);
        %     y_max = max([y1_max,y2_max,y3_max]);
        y_min = min([y2_min,y3_min]);
        y_max = max([y2_max,y3_max]);
        
        %     legend(sprintf('Match(%d)',length(x1)),...
        %         sprintf('HighMatch(%d)',length(x2)),...
        %         sprintf('LowMatch(%d)',length(x3)),...
        %         'Location','northwest','fontsize',11)
        le = legend(sprintf('High-match(%d)',length(x2)),...
            sprintf('Low-match(%d)',length(x3)),...
            'Location','northwest','fontsize',9);
        le.ItemTokenSize = ones(1,2)*10;
        
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        xlim([0 1]);
        ylim([y_min y_max+(y_max-y_min)*1]);%0.1
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-WM', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Pdf', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Match'),...
            'FontSize',12,'FontWeight','bold');
        temp_title.Interpreter = 'none';
        
        
        %% Plot mismatch
        nexttile
        
        %     x1 = meta_trialLevel_x_match;
        %     n=100;
        %     n=2^ceil(log2(n)); % round up n to the next power of 2;
        %     [pdf1,xmesh1,bandwidth1] = ksdensity(x1','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        %     %temp_ratio1 = sum(choiceMemoryBoolIndex)/sum(choiceBoolIndex);
        %     temp_ratio1 = 1;
        %     y1 = pdf1*temp_ratio1;
        %     h1 = plot(xmesh1,y1,'LineWidth',1.5,'color',[0 0 0]);
        %     hold on
        %     %h1 = histogram(x1,'Normalization','pdf');
        %     %hold on
        
        
        x2 = meta_trialLevel_x_overMismatch;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        %temp_ratio2 = sum(choiceOffloadBoolIndex)/sum(choiceBoolIndex);
        temp_ratio2 = 1;
        %     temp_ratio2 = length(meta_trialLevel_x_overMismatch)/...
        %         (length(meta_trialLevel_x_overMismatch)+length(meta_trialLevel_x_underMismatch));
        y2 = pdf2*temp_ratio2;
        h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',color_choiceMemory);
        hold on
        %h2 = histogram(x2,'Normalization','pdf');
        %hold on
        
        x3 = meta_trialLevel_x_underMismatch;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        %temp_ratio3 = sum(choiceOffloadBoolIndex)/sum(choiceBoolIndex);
        temp_ratio3 = 1;
        %     temp_ratio3 = length(meta_trialLevel_x_underMismatch)/...
        %         (length(meta_trialLevel_x_overMismatch)+length(meta_trialLevel_x_underMismatch));
        y3 = pdf3*temp_ratio3;
        h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',color_choiceOffload);
        hold on
        %h3 = histogram(x3,'Normalization','pdf');
        %hold on
        
        
        %     [x1_min,x1_max] = bounds(xmesh1);
        [x2_min,x2_max] = bounds(xmesh2);
        [x3_min,x3_max] = bounds(xmesh3);
        %     x_min = min([x1_min,x2_min,x3_min]);
        %     x_max = max([x1_max,x2_max,x3_max]);
        x_min = min([x2_min,x3_min]);
        x_max = max([x2_max,x3_max]);
        
        
        %     [y1_min,y1_max] = bounds(y1);
        [y2_min,y2_max] = bounds(y2);
        [y3_min,y3_max] = bounds(y3);
        %     y_min = min([y1_min,y2_min,y3_min]);
        %     y_max = max([y1_max,y2_max,y3_max]);
        y_min = min([y2_min,y3_min]);
        y_max = max([y2_max,y3_max]);
        
        % legend('Match','OverMismatch','UnderMismatch',...
        %     'Location','northwest','fontsize',11)
        %     legend(sprintf('Match(%d)',length(x1)),...
        %         sprintf('OverMismatch(%d)',length(x2)),...
        %         sprintf('UnderMismatch(%d)',length(x3)),...
        %         'Location','northwest','fontsize',11)
        le = legend(sprintf('Over-mismatch(%d)',length(x2)),...
            sprintf('Under-mismatch(%d)',length(x3)),...
            'Location','northwest','fontsize',9);
        le.ItemTokenSize = ones(1,2)*10;
        
        
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        xlim([0 1]);
        ylim([y_min y_max+(y_max-y_min)*1]);%0.1
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-WM', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Pdf', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Mismatch'),...
            'FontSize',12,'FontWeight','bold');
        temp_title.Interpreter = 'none';
        
    end
    
    if if_plot_B == 1
        %% Pdf of memory precision in choiceMemory trials and choiceOffload trials
        fig = figure('Name','asd','NumberTitle','off');
        %set(gcf,'Position',[10 50 240*3 240*1.1]);
        %set(gcf,'Position',[10 50 240*3 240*0.9]);        
        set(gcf,'Position',[10 50 240*3*0.73 240*0.9*0.73*0.95*0.95*0.89]);
        t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');

        %         if if_memeoryPrecision_stimuli0_response1 == 0
        %             if if_plot_fineTuningMismatch == 0
        %                 temp_str = 'stimuli, roughTuning';
        %             elseif if_plot_fineTuningMismatch == 1
        %                 temp_str = 'stimuli, fineTuning';
        %             end
        %         elseif if_memeoryPrecision_stimuli0_response1 == 1
        %             if if_plot_fineTuningMismatch == 0
        %                 temp_str = 'response, roughTuning';
        %             elseif if_plot_fineTuningMismatch == 1
        %                 temp_str = 'response, fineTuning';
        %             end
        %         end
        %
        %         if if_trainMeta_0baseline_1delay1 == 1
        %             if if_testMeta_0baselineBin_1delay1Bin == 0
        %                 t.Title.String = sprintf('Meta-memory of single trial (train delay1, test baseline) \n %s, %s',...
        %                     FOVName_currentFOV2,temp_str);
        %             elseif if_testMeta_0baselineBin_1delay1Bin == 1
        %                 %t.Title.String = sprintf('Meta-memory of single trial (train delay1, test delay1) \n %s, %s',...
        %                 %    FOVName_currentFOV2,temp_str);
        %                 t.Title.String = sprintf('%s',FOVName_currentFOV2);
        %             end
        %         elseif if_trainMeta_0baseline_1delay1 == 0
        %             if if_testMeta_0baselineBin_1delay1Bin == 0
        %                 t.Title.String = sprintf('Meta-memory of single trial (train baseline, test baseline) \n %s, %s',...
        %                     FOVName_currentFOV2,temp_str);
        %             elseif if_testMeta_0baselineBin_1delay1Bin == 1
        %                 t.Title.String = sprintf('Meta-memory of single trial (train baseline, test delay1) \n %s, %s',...
        %                     FOVName_currentFOV2,temp_str);
        %             end
        %         end
        %         t.Title.FontSize = 9;
        %         t.Title.Interpreter = 'none';
        
        %% Plot FreeChoice trials
        nexttile
        
        x2 = meta_trialLevel_x_choiceMemory;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        %[pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
        [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        temp_ratio2 = 1;
        %     temp_ratio2 = length(meta_trialLevel_x_choiceMemory)/...
        %         (length(meta_trialLevel_x_choiceMemory)+length(meta_trialLevel_x_choiceOffload));
        y2 = pdf2*temp_ratio2;
        h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',color_choiceMemory);
        hold on
        
        x3 = meta_trialLevel_x_choiceOffload;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        %[pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
        [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        temp_ratio3 = 1;
        %     temp_ratio3 = length(meta_trialLevel_x_choiceOffload)/...
        %         (length(meta_trialLevel_x_choiceMemory)+length(meta_trialLevel_x_choiceOffload));
        y3 = pdf3*temp_ratio3;
        h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',color_choiceOffload); %#ok<*NASGU>
        hold on
        
        [x2_min,x2_max] = bounds(xmesh2);
        [x3_min,x3_max] = bounds(xmesh3);
        x_min = min([x2_min,x3_min]);
        x_max = max([x2_max,x3_max]);
        
        [y2_min,y2_max] = bounds(y2);
        [y3_min,y3_max] = bounds(y3);
        y_min = min([y2_min,y3_min]);
        y_max = max([y2_max,y3_max]);
        
        plot([metaDecoderThreshold_delay1 metaDecoderThreshold_delay1],[y_min y_max+(y_max-y_min)*0.3],...
            'LineWidth',1,'color',[0 0 0]);%[0.9290 0.6940 0.1250]
        hold on
        
        
        %le = legend(sprintf('Choice-memory(%d)',length(x2)),...
        %   sprintf('Choice-offload(%d)',length(x3)),...
        %   'DecisionBoundary',...
        %   'Location','northwest','fontsize',8);
        %le.ItemTokenSize = ones(1,2)*10;
        
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        xlim([0 1]);
        ylim([y_min y_max+(y_max-y_min)*0.4]);%0.1-->0.4-->1
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-WM', 'FontSize', 9);
        ylabel('Pdf', 'FontSize', 9);
        temp_title = title(sprintf('Free-choice, %d trials',length(x2)+length(x3)),'FontSize',9);
        temp_title.Interpreter = 'none';
        
        
        
        %% Plot ChoiceMemory trials
        nexttile
        
        x2 = meta_trialLevel_x_choiceMemory;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        %[pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
        [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        temp_ratio2 = 1;
        %     temp_ratio2 = length(meta_trialLevel_x_choiceMemory)/...
        %         (length(meta_trialLevel_x_choiceMemory)+length(meta_trialLevel_x_choiceOffload));
        y2 = pdf2*temp_ratio2;
        h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',color_choiceMemory);
        hold on
        
        
        plot([metaDecoderThreshold_delay1 metaDecoderThreshold_delay1],[y_min y_max+(y_max-y_min)*0.3],...
            'LineWidth',1,'color',[0 0 0]);%[0.9290 0.6940 0.1250]
        hold on
        
        
        %le = legend(sprintf('ChoiceMemory(%d)',length(x2)),...
        %   'DecisionBoundary',...
        %   'Location','northwest','fontsize',8);
        %le.ItemTokenSize = ones(1,2)*10;
        
        tempTxt1 = sprintf('Hit\n%.2f',svm_choiceMemory_hit);
        tempTxt2 = sprintf('Miss\n%.2f',1-svm_choiceMemory_hit);
        
        if if_resample_meta == 0
            text(metaDecoderThreshold_delay1+0.15,y_min+(y_max-y_min)*1.1,tempTxt1,'Color','black','FontSize',6.5,'FontWeight','normal',...
                'HorizontalAlignment','center');
            text(metaDecoderThreshold_delay1-0.35,y_min+(y_max-y_min)*1.1,tempTxt2,'Color','black','FontSize',6.5,'FontWeight','normal',...
                'HorizontalAlignment','center');
        elseif if_resample_meta == 1
            text(metaDecoderThreshold_delay1+0.25,y_min+(y_max-y_min)*1.23,tempTxt1,'Color','black','FontSize',6.5,'FontWeight','normal',...
                'HorizontalAlignment','center');
            text(metaDecoderThreshold_delay1-0.25,y_min+(y_max-y_min)*1.23,tempTxt2,'Color','black','FontSize',6.5,'FontWeight','normal',...
                'HorizontalAlignment','center');
        end
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        xlim([0 1]);
        ylim([y_min y_max+(y_max-y_min)*0.4]);%0.1-->0.4-->1
        yticks('');
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-WM', 'FontSize', 9);
        %ylabel('Pdf', 'FontSize', 9);
        temp_title = title(sprintf('Choice-memory, %d trials',length(x2)),'FontSize',9);
        temp_title.Interpreter = 'none';

        
        %% Plot ChoiceOffload trials
        nexttile
        
        x3 = meta_trialLevel_x_choiceOffload;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        %[pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
        [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        temp_ratio3 = 1;
        %     temp_ratio3 = length(meta_trialLevel_x_choiceOffload)/...
        %         (length(meta_trialLevel_x_choiceMemory)+length(meta_trialLevel_x_choiceOffload));
        y3 = pdf3*temp_ratio3;
        h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',color_choiceOffload); %#ok<*NASGU>
        hold on
        
        
        plot([metaDecoderThreshold_delay1 metaDecoderThreshold_delay1],[y_min y_max+(y_max-y_min)*0.3],...
            'LineWidth',1,'color',[0 0 0]);%[0.9290 0.6940 0.1250]
        hold on
        
        
        %le = legend(sprintf('ChoiceOffload(%d)',length(x3)),...
        %   'DecisionBoundary',...
        %   'Location','northwest','fontsize',8);
        %le.ItemTokenSize = ones(1,2)*10;
        
        
        tempTxt1 = sprintf('False alarm\n%.2f',svm_choiceMemory_falseAlarm);
        tempTxt2 = sprintf('Correct rejection\n%.2f',1-svm_choiceMemory_falseAlarm);        
        
        
        if if_resample_meta == 0
            text(metaDecoderThreshold_delay1+0.18,y_min+(y_max-y_min)*1.1,tempTxt1,'Color','black','FontSize',6.5,'FontWeight','normal',...
                'HorizontalAlignment','center');
            text(metaDecoderThreshold_delay1-0.35,y_min+(y_max-y_min)*1.1,tempTxt2,'Color','black','FontSize',6.5,'FontWeight','normal',...
                'HorizontalAlignment','center');
        elseif if_resample_meta == 1
            text(metaDecoderThreshold_delay1+0.25,y_min+(y_max-y_min)*1.23,tempTxt1,'Color','black','FontSize',6.5,'FontWeight','normal',...
                'HorizontalAlignment','center');
            text(metaDecoderThreshold_delay1-0.25,y_min+(y_max-y_min)*1.23,tempTxt2,'Color','black','FontSize',6.5,'FontWeight','normal',...
                'HorizontalAlignment','center');
        end
        
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        xlim([0 1]);
        ylim([y_min y_max+(y_max-y_min)*0.4]);%0.1-->0.4-->1
        yticks('');
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-WM', 'FontSize', 9);
        %ylabel('Pdf', 'FontSize', 9);
        temp_title = title(sprintf('Choice-offload, %d trials',length(x3)),'FontSize',9);
        temp_title.Interpreter = 'none';
        
        
        
    end
    
    
    if if_plot_C == 1
        %% Pdf of meta-memory in choiceMemory trials and choiceOffload trials
        fig = figure('Name','asd','NumberTitle','off');
        %set(gcf,'Position',[10 260 240*3*0.73*0.33 240*0.9*0.73*0.95*0.95*0.89]);
        %set(gcf,'Position',[10 260 100 240*0.9*0.73*0.95*0.95*0.89*1.42*0.8*0.94*1.01]);        
        set(gcf,'Position',[10 260 100 240*0.9*0.73*0.95*0.95*0.89*1.42*0.8*0.94*1.01*0.88]);
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        %% Plot FreeChoice trials
        nexttile
        
        x2 = meta_trialLevel_x_choiceMemory;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        %[pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
        [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        temp_ratio2 = 1;
        %     temp_ratio2 = length(meta_trialLevel_x_choiceMemory)/...
        %         (length(meta_trialLevel_x_choiceMemory)+length(meta_trialLevel_x_choiceOffload));
        y2 = pdf2*temp_ratio2;
        h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',color_choiceMemory);
        hold on
        
        x3 = meta_trialLevel_x_choiceOffload;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        %[pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
        [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        temp_ratio3 = 1;
        %     temp_ratio3 = length(meta_trialLevel_x_choiceOffload)/...
        %         (length(meta_trialLevel_x_choiceMemory)+length(meta_trialLevel_x_choiceOffload));
        y3 = pdf3*temp_ratio3;
        h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',color_choiceOffload); %#ok<*NASGU>
        hold on
        
        [x2_min,x2_max] = bounds(xmesh2);
        [x3_min,x3_max] = bounds(xmesh3);
        x_min = min([x2_min,x3_min]);
        x_max = max([x2_max,x3_max]);
        
        [y2_min,y2_max] = bounds(y2);
        [y3_min,y3_max] = bounds(y3);
        y_min = min([y2_min,y3_min]);
        y_max = max([y2_max,y3_max]);
        
        %plot([metaDecoderThreshold_delay1 metaDecoderThreshold_delay1],[y_min y_max+(y_max-y_min)*0.3],...
        %    'LineWidth',1,'color',[0 0 0]);%[0.9290 0.6940 0.1250]
        plot([metaDecoderThreshold_delay1 metaDecoderThreshold_delay1],[y_min y_max+(y_max-y_min)*0.1],...
           'LineWidth',1,'color',[0 0 0]);%[0.9290 0.6940 0.1250]        
        hold on
        
        
        %le = legend(sprintf('Choice-memory(%d)',length(x2)),...
        %   sprintf('Choice-offload(%d)',length(x3)),...
        %   'DecisionBoundary',...
        %   'Location','northwest','fontsize',8);
        %le.ItemTokenSize = ones(1,2)*10;
        
        yticks([0 1 2]);
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        xlim([0 1]);
        ylim([y_min y_max+(y_max-y_min)*0.4]);%0.1-->0.4-->1
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-WM', 'FontSize', 9);
        ylabel('Pdf', 'FontSize', 9);
        %temp_title = title(sprintf('Free-choice, %d trials',length(x2)+length(x3)),'FontSize',9);
        %temp_title.Interpreter = 'none';
        
    end
    
    
    if if_plot_D == 1
        %% Pdf of memory precision in choiceMemoryCorrect trials and choiceMemoryError trials
        fig = figure('Name','asd','NumberTitle','off');
        set(gcf,'Position',[210 260 100 240*0.9*0.73*0.95*0.95*0.89*1.42*0.8*0.94*1.01*0.88]);
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        x2 = memoryPrecision_trialLevel_choiceMemoryCorrect;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        temp_ratio2 = 1;
        y2 = pdf2*temp_ratio2;
        h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',color_choiceMemory);
        hold on
        
        x3 = memoryPrecision_trialLevel_choiceMemoryError;
        n=100;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
        temp_ratio3 = 1;
        y3 = pdf3*temp_ratio3;
        h3 = plot(xmesh3,y3,'LineWidth',1.5,'color',color_choiceMemoryError); %#ok<*NASGU>
        hold on
        
        [x2_min,x2_max] = bounds(xmesh2);
        [x3_min,x3_max] = bounds(xmesh3);
        x_min = min([x2_min,x3_min]);
        x_max = max([x2_max,x3_max]);
        
        [y2_min,y2_max] = bounds(y2);
        [y3_min,y3_max] = bounds(y3);
        y_min = min([y2_min,y3_min]);
        y_max = max([y2_max,y3_max]);
        
        temp1 = y2;
        temp2 = y3;
        
        xmesh = xmesh3;
        
        if false
            temp_thresholdRange = 0.01:0.001:0.99; %#ok<*UNRCH>
            
            choiceMemoryCorrectError_hit_minus_falseAlarm_multi = zeros(1,length(temp_thresholdRange));
            choiceMemoryCorrectError_hit_multi = zeros(1,length(temp_thresholdRange));
            choiceMemoryCorrectError_falseAlarm_multi = zeros(1,length(temp_thresholdRange));
            
            for temptempi=1:length(temp_thresholdRange)
                
                temp_threshold = temp_thresholdRange(temptempi);
                
                %temp_hit = temp1 > temp_threshold
                
                temptempBoolIndex = xmesh > temp_threshold;
                temp_hit = sum(temp1(temptempBoolIndex))./sum(temp1);                
                temp_falseAlarm = sum(temp2(temptempBoolIndex))./sum(temp2);
                
                hit_minus_falseAlarm = temp_hit - temp_falseAlarm;
                
                choiceMemoryCorrectError_hit_multi(temptempi) = temp_hit;
                choiceMemoryCorrectError_falseAlarm_multi(temptempi) = temp_falseAlarm;
                
                choiceMemoryCorrectError_hit_minus_falseAlarm_multi(temptempi) = hit_minus_falseAlarm;
            end
            
            [M,I] = max(choiceMemoryCorrectError_hit_minus_falseAlarm_multi);
            
            tempThreshold = temp_thresholdRange(I);
            
            choiceMemoryCorrectError_hit = choiceMemoryCorrectError_hit_multi(I);
            choiceMemoryCorrectError_falseAlarm = choiceMemoryCorrectError_falseAlarm_multi(I);
            
            lowThreshold_memoryPrecision = tempThreshold;
        end
        
        
        
        plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],[y_min y_max+(y_max-y_min)*0.1],...
           'LineWidth',1,'color',[0 0 0]);      
        hold on
        
        yticks([0 1 2]);
        
        set(gca,'linewidth',1.5)
        xlim([0 1]);
        ylim([y_min y_max+(y_max-y_min)*0.4]);
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('WM strength', 'FontSize', 9);
        ylabel('Pdf', 'FontSize', 9);
        
        
    end
    
    
end
%% End