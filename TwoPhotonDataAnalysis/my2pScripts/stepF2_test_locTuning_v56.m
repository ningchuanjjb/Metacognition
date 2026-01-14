% Chuan's 2nd script (20251214)
% This script: To test location & offloading rate selectivity of all neurons from multi-FOVs.
% P.S. Please run stepF1 first.
%% Initialization


if_rsquare0_anova1_both2 = 2;

if_selective_0delay1_1delay2_2all = 0;

if_plot = 1;

if_plot_tuningStability = 1;
if_plot_maxPopTuning = 0;
if_plot_popTuningStd = 1;
if_plot_locPopTuning = 0;
if_plot_selectiveDistri_eachFOV = 0;
if_plot_rProbTuning = 0;
if_plot_tuningComparison = 0;

if_plot_length4 = 0;

if_maxPopTuning_onlySelective0_all1 = 1;

if_smoothHistogram = 1;

r2_threshold = 0;%0.05-->0.1


a = 1;
selevtivity_multiFOV;


if if_selective_0delay1_1delay2_2all == 0 || if_selective_0delay1_1delay2_2all == 2
    temp_r2 = glm_r2_lengthx_delay1Bin_multiFOV;
    temp_glm = glm_beta_lengthx_delay1Bin_multiFOV;
elseif if_selective_0delay1_1delay2_2all == 1
    temp_r2 = glm_r2_lengthx_delay2Bin_multiFOV;
    temp_glm = glm_beta_lengthx_delay2Bin_multiFOV;
end

roi_num = size(temp_glm,1);

seq_length_max = 4;
valid_length = 3;

if if_rsquare0_anova1_both2== 0
    tempBoolIndex = temp_r2 > r2_threshold;
    
    tempBoolIndex1 = sum(tempBoolIndex(:,1),2) == 1;
    tempBoolIndex2 = sum(tempBoolIndex(:,2),2) == 1;
    tempBoolIndex3 = sum(tempBoolIndex(:,3),2) == 1;
    tempBoolIndex4 = sum(tempBoolIndex(:,4),2) == 1;
else
    if if_rsquare0_anova1_both2 == 1
        if if_selective_0delay1_1delay2_2all == 2
            tempBoolIndex1 = selectiveCellBoolIndex_length1_all_multiFOV;
            tempBoolIndex2 = selectiveCellBoolIndex_length2_all_multiFOV;
            tempBoolIndex3 = selectiveCellBoolIndex_length3_all_multiFOV;
            tempBoolIndex4 = selectiveCellBoolIndex_length4_all_multiFOV;
        elseif if_selective_0delay1_1delay2_2all == 0
            tempBoolIndex1 = selectiveCellBoolIndex_length1_delay1Bin_multiFOV;
            tempBoolIndex2 = selectiveCellBoolIndex_length2_delay1Bin_multiFOV;
            tempBoolIndex3 = selectiveCellBoolIndex_length3_delay1Bin_multiFOV;
            tempBoolIndex4 = selectiveCellBoolIndex_length4_delay1Bin_multiFOV;
        elseif if_selective_0delay1_1delay2_2all == 1
            tempBoolIndex1 = selectiveCellBoolIndex_length1_delay2Bin_multiFOV;
            tempBoolIndex2 = selectiveCellBoolIndex_length2_delay2Bin_multiFOV;
            tempBoolIndex3 = selectiveCellBoolIndex_length3_delay2Bin_multiFOV;
            tempBoolIndex4 = selectiveCellBoolIndex_length4_delay2Bin_multiFOV;
        end
        
    elseif if_rsquare0_anova1_both2 == 2
        tempBoolIndex = temp_r2 > r2_threshold;
        if if_selective_0delay1_1delay2_2all == 2
            tempBoolIndex1 = selectiveCellBoolIndex_length1_all_multiFOV & sum(tempBoolIndex(:,1),2) == 1;
            tempBoolIndex2 = selectiveCellBoolIndex_length2_all_multiFOV & sum(tempBoolIndex(:,2),2) == 1;
            tempBoolIndex3 = selectiveCellBoolIndex_length3_all_multiFOV & sum(tempBoolIndex(:,3),2) == 1;
            tempBoolIndex4 = selectiveCellBoolIndex_length4_all_multiFOV & sum(tempBoolIndex(:,4),2) == 1;
        elseif if_selective_0delay1_1delay2_2all == 0
            tempBoolIndex1 = selectiveCellBoolIndex_length1_delay1Bin_multiFOV & sum(tempBoolIndex(:,1),2) == 1;
            tempBoolIndex2 = selectiveCellBoolIndex_length2_delay1Bin_multiFOV & sum(tempBoolIndex(:,2),2) == 1;
            tempBoolIndex3 = selectiveCellBoolIndex_length3_delay1Bin_multiFOV & sum(tempBoolIndex(:,3),2) == 1;
            tempBoolIndex4 = selectiveCellBoolIndex_length4_delay1Bin_multiFOV & sum(tempBoolIndex(:,4),2) == 1;
        elseif if_selective_0delay1_1delay2_2all == 1
            tempBoolIndex1 = selectiveCellBoolIndex_length1_delay2Bin_multiFOV & sum(tempBoolIndex(:,1),2) == 1;
            tempBoolIndex2 = selectiveCellBoolIndex_length2_delay2Bin_multiFOV & sum(tempBoolIndex(:,2),2) == 1;
            tempBoolIndex3 = selectiveCellBoolIndex_length3_delay2Bin_multiFOV & sum(tempBoolIndex(:,3),2) == 1;
            tempBoolIndex4 = selectiveCellBoolIndex_length4_delay2Bin_multiFOV & sum(tempBoolIndex(:,4),2) == 1;
        end
        
    end
    
end


tempBoolIndex_length = [tempBoolIndex1,tempBoolIndex2,tempBoolIndex3,tempBoolIndex4];

tempBoolIndex12_and = tempBoolIndex1 & tempBoolIndex2;
tempBoolIndex13_and = tempBoolIndex1 & tempBoolIndex3;
tempBoolIndex14_and = tempBoolIndex1 & tempBoolIndex4;
tempBoolIndex23_and = tempBoolIndex2 & tempBoolIndex3;
tempBoolIndex123_and = tempBoolIndex1 & tempBoolIndex2 & tempBoolIndex3;
tempBoolIndex1234_and = tempBoolIndex1 & tempBoolIndex2 & tempBoolIndex3 & tempBoolIndex4;

tempBoolIndex123_or = tempBoolIndex1 | tempBoolIndex2 | tempBoolIndex3;
tempBoolIndex1234_or = tempBoolIndex1 | tempBoolIndex2 | tempBoolIndex3 | tempBoolIndex4;



fprintf('length x r2 valid roi num is [%d, %d, %d, %d].\n',sum(tempBoolIndex1),sum(tempBoolIndex2),sum(tempBoolIndex3),sum(tempBoolIndex4));
fprintf('length [1&2, 1&3, 2&3, 1&4] r2 valid roi num is [%d, %d, %d, %d].\n',sum(tempBoolIndex12_and),sum(tempBoolIndex13_and),sum(tempBoolIndex23_and),sum(tempBoolIndex14_and));
fprintf('length 123 r2 valid roi num is [%d(and), %d(or)].\n',sum(tempBoolIndex123_and),sum(tempBoolIndex123_or));
fprintf('length 1234 r2 valid roi num is [%d(and), %d(or)].\n',sum(tempBoolIndex1234_and),sum(tempBoolIndex1234_or));


%% Get peak location, peak beta and beta std of each roi
% temp_glm_peakLoc = zeros(roi_num,valid_length);
temp_glm_peakLoc = zeros(roi_num,seq_length_max);
temp_glm_peakBeta = zeros(size(temp_glm_peakLoc));
temp_glm_std = zeros(size(temp_glm_peakLoc));
for tempi=1:seq_length_max
    temp_range = ((tempi-1)*numFrames+1):tempi*numFrames;
    for tempj=1:roi_num
        if tempBoolIndex_length(tempj,tempi) == false
            continue
        end
        temp_beta = temp_glm(tempj,temp_range);
        [M,I] = max(temp_beta);
        temp_glm_peakLoc(tempj,tempi) = I;
        
        temp_glm_peakBeta(tempj,tempi) = M;
        
        temp_glm_std(tempj,tempi) = std(temp_beta);
    end
end

%% Compare tuning property of across-length shared roi
temp_glm_peakLoc;
temp_glm_peakLoc_length12 = temp_glm_peakLoc(tempBoolIndex12_and,:);
temp_glm_peakLoc_length13 = temp_glm_peakLoc(tempBoolIndex13_and,:);
temp_glm_peakLoc_length123 = temp_glm_peakLoc(tempBoolIndex123_and,:);

temp_glm_peakLocDis_length12 = abs(temp_glm_peakLoc_length12(:,1)-temp_glm_peakLoc_length12(:,2));
temp_glm_peakLocDis_length12 = min([temp_glm_peakLocDis_length12,numFrames-temp_glm_peakLocDis_length12],[],2);

temp_glm_peakLocDis_length13 = abs(temp_glm_peakLoc_length13(:,1)-temp_glm_peakLoc_length13(:,3));
temp_glm_peakLocDis_length13 = min([temp_glm_peakLocDis_length13,numFrames-temp_glm_peakLocDis_length13],[],2);


temp_glm_peakBeta;
temp_glm_peakBeta_length12 = temp_glm_peakBeta(tempBoolIndex12_and,:);
temp_glm_peakBeta_length13 = temp_glm_peakBeta(tempBoolIndex13_and,:);
temp_glm_peakBeta_length123 = temp_glm_peakBeta(tempBoolIndex123_and,:);
temp_glm_peakBeta_length1234 = temp_glm_peakBeta(tempBoolIndex1234_and,:);
% temp_glm_peakBeta_length1234 = temp_glm_peakBeta_length123;

[~,p12,~,~] = ttest(temp_glm_peakBeta_length12(:,1),temp_glm_peakBeta_length12(:,2));

[~,p13,~,~] = ttest(temp_glm_peakBeta_length13(:,1),temp_glm_peakBeta_length13(:,3));

% cftool();

a = 1;
temp_glm_peakBeta_length123;
x = [];
y = [];
for tempi=1:valid_length
    x = [x;tempi*ones(size(temp_glm_peakBeta_length123,1),1)];    %#ok<*AGROW>
    y = [y;temp_glm_peakBeta_length123(:,tempi)];
end
%cftool(x,y);
temp_mdl = fitglm(x,y,'linear');
p123 = temp_mdl.Coefficients.pValue(2);

a = 1;
temp_glm_peakBeta_length1234;
x = [];
y = [];
for tempi=1:seq_length_max
    x = [x;tempi*ones(size(temp_glm_peakBeta_length1234,1),1)];    %#ok<*AGROW>
    y = [y;temp_glm_peakBeta_length1234(:,tempi)];
end
%cftool(x,y);
temp_mdl2 = fitglm(x,y,'linear');
p1234 = temp_mdl2.Coefficients.pValue(2);

temp_glm;

temp_glm_std;
temp_glm_std_length12 = temp_glm_std(tempBoolIndex12_and,:);
temp_glm_std_length13 = temp_glm_std(tempBoolIndex13_and,:);
temp_glm_std_length123 = temp_glm_std(tempBoolIndex123_and,:);
temp_glm_std_length1234 = temp_glm_std(tempBoolIndex1234_and,:);
% temp_glm_std_length1234 = temp_glm_std_length123;

[~,p12_std,~,~] = ttest(temp_glm_std_length12(:,1),temp_glm_std_length12(:,2));
[~,p13_std,~,~] = ttest(temp_glm_std_length13(:,1),temp_glm_std_length13(:,3));

x = [];
y = [];
for tempi=1:valid_length
    x = [x;tempi*ones(size(temp_glm_std_length123,1),1)];    %#ok<*AGROW>
    y = [y;temp_glm_std_length123(:,tempi)];
end
%cftool(x,y);
temp_mdl = fitglm(x,y,'linear');
p123_std = temp_mdl.Coefficients.pValue(2);

x = [];
y = [];
for tempi=1:seq_length_max
    x = [x;tempi*ones(size(temp_glm_std_length1234,1),1)];    %#ok<*AGROW>
    y = [y;temp_glm_std_length1234(:,tempi)];
end
%cftool(x,y);
temp_mdl2 = fitglm(x,y,'linear');
p1234_std = temp_mdl2.Coefficients.pValue(2);

a = 1;

%% Compare population tuning between length in each location
% valid_length = 3;

p_population_loc_12 = zeros(1,numFrames);
p_population_loc_13 = zeros(1,numFrames);
p_population_loc_23 = zeros(1,numFrames);

temp_beta_population_loc = cell(numFrames,valid_length);
for tempLoc=1:numFrames
    for tempi=1:valid_length
        temp_range = ((tempi-1)*numFrames+1):tempi*numFrames;
        temp_beta_population_loc{tempLoc,tempi} = temp_glm(tempBoolIndex_length(:,tempi),temp_range(tempLoc));
    end
    a = 1; %#ok<*NASGU>
    
    %[~,p_population_loc_12(tempLoc),~,~] = ttest2(temp_beta_population{1},temp_beta_population{2});
    [~,p_population_loc_12(tempLoc),~,~] = ttest2(temp_beta_population_loc{tempLoc,1},temp_beta_population_loc{tempLoc,2},'Tail','right');
    [~,p_population_loc_13(tempLoc),~,~] = ttest2(temp_beta_population_loc{tempLoc,1},temp_beta_population_loc{tempLoc,3},'Tail','right');
    [~,p_population_loc_23(tempLoc),~,~] = ttest2(temp_beta_population_loc{tempLoc,2},temp_beta_population_loc{tempLoc,3},'Tail','right');
end

p_population_loc_abs_12 = zeros(1,numFrames);
p_population_loc_abs_13 = zeros(1,numFrames);
p_population_loc_abs_23 = zeros(1,numFrames);

temp_beta_population_abs_loc = cell(numFrames,valid_length);
for tempLoc=1:numFrames
    for tempi=1:valid_length
        temp_range = ((tempi-1)*numFrames+1):tempi*numFrames;
        %temp_beta_population_positive_loc{tempLoc,tempi} = temp_glm(tempBoolIndex_length(:,tempi),temp_range(tempLoc));
        temp1 = temp_glm(tempBoolIndex_length(:,tempi),temp_range(tempLoc));
        %temp1 = temp1(temp1>0);
        %temp1 = temp1(temp1~=0);
        temp1 = abs(temp1);
        temp_beta_population_abs_loc{tempLoc,tempi} = temp1;
    end
    a = 1; %#ok<*NASGU>
    
    [~,p_population_loc_abs_12(tempLoc),~,~] = ttest2(temp_beta_population_abs_loc{tempLoc,1},temp_beta_population_abs_loc{tempLoc,2},'Tail','right');
    [~,p_population_loc_abs_13(tempLoc),~,~] = ttest2(temp_beta_population_abs_loc{tempLoc,1},temp_beta_population_abs_loc{tempLoc,3},'Tail','right');
    [~,p_population_loc_abs_23(tempLoc),~,~] = ttest2(temp_beta_population_abs_loc{tempLoc,2},temp_beta_population_abs_loc{tempLoc,3},'Tail','right');
end


a = 1;

%% Population max tuning
% temp_beta_population_peakBeta = cell(1,valid_length);
temp_beta_population_peakBeta = cell(1,seq_length_max);
temp_beta_population_peakLoc = cell(size(temp_beta_population_peakBeta));
temp_beta_population_std = cell(size(temp_beta_population_peakBeta));
% for tempi=1:valid_length
for tempi=1:seq_length_max
    temp_range = ((tempi-1)*numFrames+1):tempi*numFrames;
    if if_maxPopTuning_onlySelective0_all1 == 0
        temp_beta1 = temp_glm(tempBoolIndex_length(:,tempi),temp_range);
    elseif if_maxPopTuning_onlySelective0_all1 == 1
        temp_beta1 = temp_glm(tempBoolIndex123_or,temp_range);
    end
    [M,I] = max(temp_beta1,[],2);
    
    temp_beta_population_peakBeta{tempi} = M;
    temp_beta_population_peakLoc{tempi} = I;
    temp_beta_population_std{tempi} = std(temp_beta1,0,2);
end
a = 1;

[~,p_population_peakBeta_12,~,~] = ttest2(temp_beta_population_peakBeta{1},temp_beta_population_peakBeta{2},'Tail','right');
[~,p_population_peakBeta_13,~,~] = ttest2(temp_beta_population_peakBeta{1},temp_beta_population_peakBeta{3},'Tail','right');
[~,p_population_peakBeta_23,~,~] = ttest2(temp_beta_population_peakBeta{2},temp_beta_population_peakBeta{3},'Tail','right');

[~,p_population_peakBeta_14,~,~] = ttest2(temp_beta_population_peakBeta{1},temp_beta_population_peakBeta{4},'Tail','right');
[~,p_population_peakBeta_24,~,~] = ttest2(temp_beta_population_peakBeta{2},temp_beta_population_peakBeta{4},'Tail','right');
[~,p_population_peakBeta_34,~,~] = ttest2(temp_beta_population_peakBeta{3},temp_beta_population_peakBeta{4},'Tail','right');


temp_beta_population_peakBeta;
x = [];
y = [];
% for tempi=1:valid_length
for tempi=1:seq_length_max
    x = [x;tempi*ones(size(temp_beta_population_peakBeta{tempi},1),1)];    %#ok<*AGROW>
    y = [y;temp_beta_population_peakBeta{tempi}];
end
temp_mdl = fitglm(x,y,'linear');
p_population_peakBeta_123 = temp_mdl.Coefficients.pValue(2);
% p_population_peakBeta_1234 = temp_mdl.Coefficients.pValue(2);

a = 1;
% [mean(temp_beta_population_peakBeta{1}), mean(temp_beta_population_peakBeta{2}), mean(temp_beta_population_peakBeta{3}), mean(temp_beta_population_peakBeta{4})]

a = 1;

[~,p_population_std_12,~,~] = ttest2(temp_beta_population_std{1},temp_beta_population_std{2},'Tail','right');
[~,p_population_std_13,~,~] = ttest2(temp_beta_population_std{1},temp_beta_population_std{3},'Tail','right');
[~,p_population_std_23,~,~] = ttest2(temp_beta_population_std{2},temp_beta_population_std{3},'Tail','right');

[~,p_population_std_14,~,~] = ttest2(temp_beta_population_std{1},temp_beta_population_std{4},'Tail','right');
[~,p_population_std_24,~,~] = ttest2(temp_beta_population_std{2},temp_beta_population_std{4},'Tail','right');
[~,p_population_std_34,~,~] = ttest2(temp_beta_population_std{3},temp_beta_population_std{4},'Tail','right');


%% Example roi std
if false
    
    temp_FOVIndex = 20;%20    
    temp_cellIndex_suit2p_target = 19;%25,19
    
    %temp_FOVIndex = 10;%
    %temp_cellIndex_suit2p_target = 78;%
    
    temp_cellRange = selevtivity_multiFOV.FOVAllCellRange_multiFOV(temp_FOVIndex,:); %#ok<*UNRCH>
    temp_cellIndex_suite2p = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_cellRange(1):temp_cellRange(2));
    temp_cellIndex_suite2p;
    

    temp_cellIndex_target = find(temp_cellIndex_suite2p==temp_cellIndex_suit2p_target);
    
    temp_cellIndex_target
    temp_cellIndex_target_global = temp_cellRange(1)+temp_cellIndex_target-1;
    
    tempBoolIndex123_or;
    
    temptempIndex = find(tempBoolIndex123_or==true);    
    
    temptempIndex2 = find(temptempIndex==temp_cellIndex_target_global);
    
    temp_beta_population_std_target = ...
        [temp_beta_population_std{1}(temptempIndex2),...
        temp_beta_population_std{2}(temptempIndex2),...
        temp_beta_population_std{3}(temptempIndex2)];
    
    
    fig = figure('Name','Example neuron tuning strength (std)','NumberTitle','off');
    %set(gcf,'Position',[10 50 200 165]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 50 112*0.7 80*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    plot(1:3,temp_beta_population_std_target,'color',[0.25 0.25 0.25],'linewidth',1.5,'marker','o',...
        'MarkerFaceColor',[0.25 0.25 0.25],'MarkerSize',3);
    hold on
    
    temp_y_min = min(temp_beta_population_std_target);
    temp_y_max = max(temp_beta_population_std_target);
    
    xlim([0.5 3.5])
    
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.2 temp_y_max+(temp_y_max-temp_y_min)*0.2]);
    
    temptemp_stepNum = 3;
    
    temptemp_range = linspace(temp_y_min,temp_y_max,temptemp_stepNum);
    %temptemp_step = floor((temptemp_range(2)-temptemp_range(1))*100)/100;
    temptemp_step = ceil((temptemp_range(2)-temptemp_range(1))*10)/10;
    temp1 = temptemp_range(1);
    temp2 = floor(temp1*10)/10;
    
    set(gca,'ytick',temp2+(0:temptemp_step:temptemp_step*temptemp_stepNum))
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'XTick', [1 2 3]);
    set(gca,'XTickLabel', ['1';'2';'3'],'FontSize', 10);
    set(gca,'box','off');
    %ylabel('Tuning strength', 'FontSize', 15, 'FontWeight', 'bold');
    title(sprintf('Tuning strength'),'fontsize',10);
    
end

%% rProb_glm

% temp1 = sum(selectiveCellBoolIndex_rProb_glm_delay2_multiFOV);
% temp1_and = sum(tempBoolIndex123_or & selectiveCellBoolIndex_rProb_glm_delay2_multiFOV);

if if_selective_0delay1_1delay2_2all == 0
    selectiveCellBoolIndex_rProb_glm_delayx_multiFOV = selectiveCellBoolIndex_rProb_glm_delay1_multiFOV;
    r2_rProb_glm_delayx_multiFOV = r2_rProb_glm_delay1_multiFOV;
    p_rProb_glm_delayx_multiFOV = p_rProb_glm_delay1_multiFOV;
    
    selectiveCellBoolIndex_precision_glm_delayx_multiFOV = selectiveCellBoolIndex_precision_glm_delay1_multiFOV;
    r2_precision_glm_delayx_multiFOV = r2_precision_glm_delay1_multiFOV;
    p_precision_glm_delayx_multiFOV = p_precision_glm_delay1_multiFOV;
    
elseif if_selective_0delay1_1delay2_2all == 1
    selectiveCellBoolIndex_rProb_glm_delayx_multiFOV = selectiveCellBoolIndex_rProb_glm_delay2_multiFOV;
    r2_rProb_glm_delayx_multiFOV = r2_rProb_glm_delay2_multiFOV;
    p_rProb_glm_delayx_multiFOV = p_rProb_glm_delay2_multiFOV;
    
    
    selectiveCellBoolIndex_precision_glm_delayx_multiFOV = selectiveCellBoolIndex_precision_glm_delay2_multiFOV;
    r2_precision_glm_delayx_multiFOV = r2_precision_glm_delay2_multiFOV;
    p_precision_glm_delayx_multiFOV = p_precision_glm_delay2_multiFOV;
end

tempL = tempBoolIndex123_or;

tempA = selectiveCellBoolIndex_rProb_glm_delayx_multiFOV;
tempAL_and = tempL & selectiveCellBoolIndex_rProb_glm_delayx_multiFOV;

tempB = selectiveCellBoolIndex_precision_glm_delayx_multiFOV;
tempBL_and = tempL & selectiveCellBoolIndex_precision_glm_delayx_multiFOV;

tempAB_and = tempA & tempB;
tempAB_or = tempA | tempB;

tempABL_and = tempA & tempB & tempL;

fprintf('length 123(or) & precision roi num is %d (%d & %d).\n',sum(tempBL_and),sum(tempL),sum(tempB));
fprintf('length 123(or) & offloading roi num is %d (%d & %d).\n',sum(tempAL_and),sum(tempL),sum(tempA));
fprintf('precision & offloading roi num is %d (%d & %d).\n',sum(tempAB_and),sum(tempB),sum(tempA));

fprintf('length 123(or) & precision roi & offloading roi num is %d.\n',sum(tempABL_and));

tempAB_pureA = tempA & (~tempB);
tempAB_pureB = (~tempA) & tempB;


tempBL_pureB = tempB & (~tempL);
tempBL_pureL = (~tempB) & tempL;


tempAL_pureA = tempA & (~tempL);
tempAL_pureL = (~tempA) & tempL;

a = 1;


selectiveCellBoolIndex_rProb_glm_baseline_multiFOV;
selectiveCellBoolIndex_rProb_glm_delay1_multiFOV;
selectiveCellBoolIndex_rProb_glm_delay2_multiFOV;



%% selectiveCell num in each FOV
tempBoolIndex123_or_num_eachFOV = zeros(size(FOVAllCellRange_multiFOV,1),1);
selectiveCellBoolIndex_rProb_glm_delayx_num_eachFOV = zeros(size(FOVAllCellRange_multiFOV,1),1);
for tempi=1:size(FOVAllCellRange_multiFOV,1)
    temp_range = FOVAllCellRange_multiFOV(tempi,1):FOVAllCellRange_multiFOV(tempi,2);
    
    temptempBoolIndex = tempBoolIndex123_or(temp_range);
    tempBoolIndex123_or_num_eachFOV(tempi) = sum(temptempBoolIndex);
    
    temptempBoolIndex = selectiveCellBoolIndex_rProb_glm_delayx_multiFOV(temp_range);
    selectiveCellBoolIndex_rProb_glm_delayx_num_eachFOV(tempi) = sum(temptempBoolIndex);
end
tempBoolIndex123_or_num_eachFOV;
selectiveCellBoolIndex_rProb_glm_delayx_num_eachFOV;


%% Plot
if if_plot == 1
    close all
    if if_plot_tuningStability == 1
        %% Fig, single roi tuning stability
        fig = figure('Name','tuning stability','NumberTitle','off');
        set(gcf,'Position',[10 50 1100 820]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        histogram(temp_glm_peakLocDis_length12,'FaceColor',[1 1 1]*0.25);
        
        set(gca, 'FontSize', 14)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Tuning shift', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('Roi count', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('length 12 max tuning shift \n Num=%d',sum(tempBoolIndex12_and)),'fontsize',16);
        
        
        nexttile
        histogram(temp_glm_peakLocDis_length13,'FaceColor',[1 1 1]*0.25);
        
        set(gca, 'FontSize', 14)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Tuning shift', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('Roi count', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('length 13 max tuning shift \n Num=%d',sum(tempBoolIndex13_and)),'fontsize',16);
        
        %nexttile
        %set(gca,'visible','off');
        
        a = 1;
        nexttile
        
        %valid_length
        for tempi=1:valid_length
            %h = histogram(temp_glm_peakBeta_length123(:,tempi),'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.5);
            
            if if_smoothHistogram == 0
                h = histogram(temp_glm_peakBeta_length123(:,tempi),'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.5);
                
            elseif if_smoothHistogram == 1
                x = temp_glm_peakBeta_length123(:,tempi);
                
                temp_data = x';
                n=100;
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                [pdf2,xmesh2,bandwidth2] = ksdensity(temp_data,'NumPoints',n,'Function','pdf');
                
                plot(xmesh2,pdf2,'LineWidth',1.5);
                
            end
            
            hold on
        end
        
        legend('length 1','length 2','length 3',...
            'Location','northeast','fontsize',9)
        
        % set(gca,'YScale','log');
        set(gca, 'FontSize', 14)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Tuning strength', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('pdf', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('max tuning strength'),'fontsize',16);
        
        a = 1;
        
        nexttile
        
        temp_p = p12;
        temp_1 = temp_glm_peakBeta_length12(:,1);
        temp_2 = temp_glm_peakBeta_length12(:,2);
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        temp_y_min = 0;
        temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
        
        bar([0 1], [mean(temp_1) mean(temp_2)], ...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        errorbar([0 1], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(0.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca,'XTickLabel', ["length 1"; "length 2"],'FontSize', 12);%给坐标加标签
        set(gca, 'FontSize', 14)
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Tuning weight', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('Tuning strength', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('length 12 max tuning strength \n Num=%d, paired t test',sum(tempBoolIndex12_and)),'fontsize',16);
        
        
        nexttile
        
        temp_p = p13;
        temp_1 = temp_glm_peakBeta_length13(:,1);
        temp_2 = temp_glm_peakBeta_length13(:,3);
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        temp_y_min = 0;
        temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
        
        bar([0 1], [mean(temp_1) mean(temp_2)], ...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        errorbar([0 1], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(0.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca,'XTickLabel', ["length 1"; "length 3"],'FontSize', 12);%给坐标加标签
        set(gca, 'FontSize', 14)
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Tuning weight', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('Tuning strength', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('length 13 max tuning strength \n Num=%d, paired t test',sum(tempBoolIndex13_and)),'fontsize',16);
        
        
        nexttile
        
        temp_p = p123;
        temp_1 = temp_glm_peakBeta_length123(:,1);
        temp_2 = temp_glm_peakBeta_length123(:,2);
        temp_3 = temp_glm_peakBeta_length123(:,3);
        
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
        set(gca,'XTickLabel', ["length1"; "length2"; "length3"],'FontSize', 12);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Tuning weight', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('Tuning strength', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('length 123 max tuning strength \n Num=%d, linear regression',sum(tempBoolIndex123_and)),'fontsize',16);
        
        
        nexttile
        
        temp_p = p12_std;
        temp_1 = temp_glm_std_length12(:,1);
        temp_2 = temp_glm_std_length12(:,2);
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        temp_y_min = 0;
        temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
        
        bar([0 1], [mean(temp_1) mean(temp_2)], ...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        errorbar([0 1], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(0.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca,'XTickLabel', ["length 1"; "length 2"],'FontSize', 12);%给坐标加标签
        set(gca, 'FontSize', 14)
        set(gca,'box','off');% 取消右、上边框
        ylabel('Tuning variance (std)', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('length 12 tuning variance \n Num=%d, paired t test',sum(tempBoolIndex12_and)),'fontsize',16);
        
        
        nexttile
        
        temp_p = p13_std;
        temp_1 = temp_glm_std_length13(:,1);
        temp_2 = temp_glm_std_length13(:,3);
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        temp_y_min = 0;
        temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
        
        bar([0 1], [mean(temp_1) mean(temp_2)], ...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        errorbar([0 1], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(0.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca,'XTickLabel', ["length 1"; "length 3"],'FontSize', 12);%给坐标加标签
        set(gca, 'FontSize', 14)
        set(gca,'box','off');% 取消右、上边框
        ylabel('Tuning variance (std)', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('length 13 tuning variance \n Num=%d, paired t test',sum(tempBoolIndex13_and)),'fontsize',16);
        
        
        nexttile
        
        temp_p = p123_std;
        temp_1 = temp_glm_std_length123(:,1);
        temp_2 = temp_glm_std_length123(:,2);
        temp_3 = temp_glm_std_length123(:,3);
        
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
        set(gca,'XTickLabel', ["length1"; "length2"; "length3"],'FontSize', 12);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        ylabel('Tuning variance (std)', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('length 123 tuning variance \n Num=%d, linear regression',sum(tempBoolIndex123_and)),'fontsize',16);
        
    end
    
    if if_plot_maxPopTuning == 1
        %% Fig, population max tuning
        fig = figure('Name','population max tuning','NumberTitle','off');
        %set(gcf,'Position',[10 50 800 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 1200 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 1200 375]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
        
        
        nexttile
        if if_plot_length4 == 0
            temp_valid_length = 3;
        elseif if_plot_length4 == 1
            temp_valid_length = 4;
        end
        
        %valid_length
        %for tempi=1:temp_valid_length
        for tempi=temp_valid_length:-1:1
            if if_smoothHistogram == 0
                h = histogram(temp_beta_population_peakBeta{tempi},'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.5);
                
            elseif if_smoothHistogram == 1
                x = temp_beta_population_peakBeta{tempi};
                
                temp_data = x';
                n=100;
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                [pdf2,xmesh2,bandwidth2] = ksdensity(temp_data,'NumPoints',n,'Function','pdf');
                
                if tempi == 1
                    temp_color = [1 1 1]*0.3;
                elseif tempi == 2
                    temp_color = [1 1 1]*0.5;
                elseif tempi == 3
                    temp_color = [1 1 1]*0.8;
                end
                %plot(xmesh2,pdf2,'LineWidth',1.5);
                plot(xmesh2,pdf2,'LineWidth',1.5,'color',temp_color);
                
            end
            
            hold on
        end
        
        legend('length 3','length 2','length 1',...
            'Location','northeast','fontsize',9)
        
        % set(gca,'YScale','log');
        set(gca, 'FontSize', 14)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Tuning strength', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('pdf', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('max tuning strength'),'fontsize',16);
        
        a = 1;
        
        nexttile
        
        temp_p = p_population_peakBeta_123;
        temp_1 = temp_beta_population_peakBeta{1};
        temp_2 = temp_beta_population_peakBeta{2};
        temp_3 = temp_beta_population_peakBeta{3};
        temp_4 = temp_beta_population_peakBeta{4};
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        temp3_SEM = std(temp_3)/sqrt(length(temp_3));
        temp4_SEM = std(temp_4)/sqrt(length(temp_4));
        
        temp_y_min = 0;
        temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
        
        if if_plot_length4 == 0
            %             bar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)], ...
            %                 'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
            temp_bar = bar([0 1 2], [mean(temp_3) mean(temp_2) mean(temp_1)], ...
                'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
            hold on
            a = 1;
            temp_bar.CData(1,:) = [1 1 1]*0.8;
            temp_bar.CData(2,:) = [1 1 1]*0.5;
            temp_bar.CData(3,:) = [1 1 1]*0.3;
            
            %             errorbar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
            errorbar([0 1 2], [mean(temp_3) mean(temp_2) mean(temp_1)],[temp3_SEM temp2_SEM temp1_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
        elseif if_plot_length4 == 1
            bar([0 1 2 3], [mean(temp_4) mean(temp_3) mean(temp_2) mean(temp_1)], ...
                'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
            hold on
            errorbar([0 1 2 3], [mean(temp_4) mean(temp_3) mean(temp_2) mean(temp_1)],[temp4_SEM temp3_SEM temp2_SEM temp1_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
        end
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
        if if_plot_length4 == 0
            %set(gca,'XTickLabel', ["length1"; "length2"; "length3"],'FontSize', 12);%给坐标加标签
            set(gca,'XTickLabel', ["length3"; "length2"; "length1"],'FontSize', 12);%给坐标加标签
        elseif if_plot_length4 == 1
            set(gca,'XTickLabel', ["length4"; "length3"; "length2"; "length1"],'FontSize', 12);%给坐标加标签
        end
        set(gca,'box','off');% 取消右、上边框
        ylabel('Tuning strength', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('max tuning strength \nlinear regression'),'fontsize',16);
        
        
        nexttile
        
        %valid_length
        count_loc = zeros(valid_length,numFrames);
        for tempi=1:valid_length
            a = 1;
            temp_beta1 = temp_beta_population_peakLoc{tempi};
            
            temp_count_loc = zeros(1,numFrames);
            for tempj=1:numFrames
                temp_count_loc(tempj) = sum(temp_beta1==tempj);
            end
            count_loc(tempi,:) = temp_count_loc;
            
            plot(temp_count_loc./sum(temp_count_loc),'LineWidth',1);
            %h = histogram(temp_beta_population_peakLoc{tempi},'FaceAlpha',0.5,'Normalization','pdf','DisplayStyle','bar','LineWidth',0.5);
            hold on
        end
        
        mean_pdf_loc = sum(count_loc,1)./sum(count_loc,'all');
        
        x = 1:numFrames;
        y = mean_pdf_loc;
        %cftool(x,y);
        
        xq = 1:0.001:numFrames;
        yq = makima(x,y,xq);
        plot(xq,yq,'LineWidth',1.5,'Color',[0 0 0]);
        hold on
        
        %plot(sum(count_loc,1)./sum(count_loc,'all'),'LineWidth',1);
        %hold on
        
        legend('length 1','length 2','length 3','fit',...
            'Location','northeast','fontsize',9)
        
        
        xlim([1-0.5 numFrames+0.5]);
        
        % set(gca,'YScale','log');
        set(gca, 'FontSize', 14)
        set(gca,'box','off');% 取消右、上边框
        set(gca,'XTick',1:numFrames);%设置要显示坐标刻度的范围
        xlabel('Location', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('pdf', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('max tuning location'),'fontsize',16);
        
    end
    
    %if_plot_popTuningStd = 1;
    
    if if_plot_popTuningStd == 1
        %% Fig, population tuning std
        fig = figure('Name','population tuning std','NumberTitle','off');
        %set(gcf,'Position',[10 50 800 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 1200 355]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 800 375]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 700 190]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 700 190*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 933 190*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
        
        
        nexttile
        if if_plot_length4 == 0
            temp_valid_length = 3;
        elseif if_plot_length4 == 1
            temp_valid_length = 4;
        end
        
        temp_x_max = 0;
        
        %valid_length
        for tempi=1:temp_valid_length
            %for tempi=temp_valid_length:-1:1
            if if_smoothHistogram == 0
                h = histogram(temp_beta_population_std{tempi},'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.5);
                
            elseif if_smoothHistogram == 1
                x = temp_beta_population_std{tempi};
                
                temp_data = x';
                n=100;
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                %[pdf2,xmesh2,bandwidth2] = ksdensity(temp_data,'NumPoints',n,'Function','pdf');
                [pdf2,xmesh2,bandwidth2] = ksdensity(temp_data,'NumPoints',n,'Function','pdf','Support',[-0.01 5.01],'BoundaryCorrection', 'Reflection');
                
                if tempi == 1
                    temp_color = [1 1 1]*0.3;
                elseif tempi == 2
                    temp_color = [1 1 1]*0.5;
                elseif tempi == 3
                    temp_color = [1 1 1]*0.8;
                end
                %plot(xmesh2,pdf2,'LineWidth',1.5);
                plot(xmesh2,pdf2,'LineWidth',1.5,'color',temp_color);
                
                temp_x_max = max(temp_x_max,max(xmesh2));
                
            end
            
            hold on
        end
        
        le = legend('length 1','length 2','length 3',...
            'Location','northeast','fontsize',9);
        le.ItemTokenSize = ones(1,3)*10;
        
        xlim([0 temp_x_max]);
        
        set(gca,'linewidth',1.5)
        % set(gca,'YScale','log');
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Tuning strength (std)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('pdf', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Tuning strength (std)'),'fontsize',12);
        
        a = 1;
        
        nexttile
        
        temp_p = p_population_peakBeta_123;
        temp_1 = temp_beta_population_std{1};
        temp_2 = temp_beta_population_std{2};
        temp_3 = temp_beta_population_std{3};
        temp_4 = temp_beta_population_std{4};
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        temp3_SEM = std(temp_3)/sqrt(length(temp_3));
        temp4_SEM = std(temp_4)/sqrt(length(temp_4));
        
        temp_y_min = 0;
        temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
        
        if if_plot_length4 == 0
            %             bar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)], ...
            %                 'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
            temp_bar = bar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)], ...
                'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
            hold on
            a = 1;
            temp_bar.CData(1,:) = [1 1 1]*0.3;
            temp_bar.CData(2,:) = [1 1 1]*0.5;
            temp_bar.CData(3,:) = [1 1 1]*0.8;
            
            errorbar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
            %             errorbar([0 1 2], [mean(temp_3) mean(temp_2) mean(temp_1)],[temp3_SEM temp2_SEM temp1_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
        elseif if_plot_length4 == 1
            bar([0 1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3) mean(temp_4)], ...
                'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
            hold on
            errorbar([0 1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3) mean(temp_4)],[temp1_SEM temp2_SEM temp3_SEM temp4_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
        end
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca, 'FontSize', 12)
        if if_plot_length4 == 0
            set(gca,'XTickLabel', ["length1"; "length2"; "length3"],'FontSize', 12);%给坐标加标签
            %             set(gca,'XTickLabel', ["length3"; "length2"; "length1"],'FontSize', 12);%给坐标加标签
        elseif if_plot_length4 == 1
            set(gca,'XTickLabel', ["length1"; "length2"; "length3"; "length4"],'FontSize', 12);%给坐标加标签
        end
        set(gca,'box','off');% 取消右、上边框
        ylabel('Tuning strength', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Tuning strength (std) \nlinear regression'),'fontsize',12);
        
        
        nexttile
        
        temp_p = p_population_peakBeta_123;
        temp_1 = temp_beta_population_std{1};
        temp_2 = temp_beta_population_std{2};
        temp_3 = temp_beta_population_std{3};
        temp_4 = temp_beta_population_std{4};
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        temp3_SEM = std(temp_3)/sqrt(length(temp_3));
        temp4_SEM = std(temp_4)/sqrt(length(temp_4));
        
        temp_y_min = 0;
        %temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
        temp_y_max = max([temp_1;temp_2;temp_3]);
        
        if if_plot_length4 == 0
            %             temp_bar = bar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)], ...
            %                 'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1);
            %             hold on
            %             a = 1;
            %             temp_bar.CData(1,:) = [1 1 1]*0.3;
            %             temp_bar.CData(2,:) = [1 1 1]*0.5;
            %             temp_bar.CData(3,:) = [1 1 1]*0.8;
            %
            %             errorbar([0 1 2], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
            
            %boxplot([temp_1 temp_2 temp_3]);
            %violinplot([temp_1 temp_2 temp_3]);
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 3, 1);
            
            %             temptemp_color2 = [[1 1 1]*0.3;...
            %                 [1 1 1]*0.5;...
            %                 [1 1 1]*0.8];
            violinplot([temp_1 temp_2 temp_3],[],'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2);
            
            hold on
            
            
        elseif if_plot_length4 == 1
            %
        end
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(2,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.5 3.5])
        %ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca, 'FontSize', 12)
        if if_plot_length4 == 0
            set(gca,'XTickLabel', ["length1"; "length2"; "length3"],'FontSize', 11);%给坐标加标签
        elseif if_plot_length4 == 1
            set(gca,'XTickLabel', ["length1"; "length2"; "length3"; "length4"],'FontSize', 11);%给坐标加标签
        end
        xtickangle(0);
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Tuning strength', 'FontSize', 12, 'FontWeight', 'bold');
        %title(sprintf('Tuning strength (std) \nlinear regression'),'fontsize',12);
        title(sprintf('Tuning strength'),'fontsize',12);
        
        
        
        nexttile
        
        temp_1 = sum(tempBoolIndex1);
        temp_2 = sum(tempBoolIndex2);
        temp_3 = sum(tempBoolIndex3);
        temp_4 = sum(tempBoolIndex4);
        
        
        temp_y_min = 0;
        temp_y_max = max([temp_1;temp_2;temp_3]);
        
        if if_plot_length4 == 0
            %temptemp_color1 = [1 1 1]*0.5;
            %temptemp_color2 = repmat(temptemp_color1, 3, 1);
            %violinplot([temp_1 temp_2 temp_3],[],'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2);
            
            bar([0 1 2], [temp_1 temp_2 temp_3],0.5, ...
                'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
            hold on
            
            
        elseif if_plot_length4 == 1
            %
        end
        
        set(gca,'linewidth',1.5)
        xlim([-0.5 2.5])
        %ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca, 'FontSize', 12)
        if if_plot_length4 == 0
            set(gca,'XTickLabel', ["length1"; "length2"; "length3"],'FontSize', 11);%给坐标加标签
        elseif if_plot_length4 == 1
            set(gca,'XTickLabel', ["length1"; "length2"; "length3"; "length4"],'FontSize', 11);%给坐标加标签
        end
        xtickangle(0);
        
        temptemp_stepNum = 5;
        
        temptemp_range = linspace(0,temp_y_max,(temptemp_stepNum-1));
        temptemp_step = floor((temptemp_range(2)-temptemp_range(1))*0.01)/0.01;
        
        set(gca,'ytick',0:temptemp_step:temptemp_step*temptemp_stepNum)
        
        
        set(gca,'box','off');% 取消右、上边框
        %ylabel('Tuning strength', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Selective neuron number'),'fontsize',12);
        
    end
    
    
    if if_plot_locPopTuning == 1
        %% Fig, population tuning
        fig = figure('Name','population tuning','NumberTitle','off');
        set(gcf,'Position',[10 50 1100 700]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
        
        for tempLoc=1:numFrames
            nexttile
            
            p_population_loc_12;
            p_population_loc_13;
            p_population_loc_23;
            
            a = 1;
            
            for tempi=1:valid_length
                h = histogram(temp_beta_population_loc{tempLoc,tempi},'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1);
                hold on
            end
            if tempLoc == 1
                legend('length 1','length 2','length 3',...
                    'Location','northeast','fontsize',9)
                % legend
            end
            
            ylim([0 2]);
            %set(gca,'YScale','log');
            set(gca, 'FontSize', 14)
            set(gca,'box','off');% 取消右、上边框
            xlabel('Tuning weight', 'FontSize', 15, 'FontWeight', 'bold');
            ylabel('pdf', 'FontSize', 15, 'FontWeight', 'bold');
            title(sprintf('location %d',tempLoc),'fontsize',16);
        end
    end
    a = 1;
    
    if if_plot_selectiveDistri_eachFOV == 1
        %% Fig, selectiveDistri of each FOV
        fig = figure('Name','selectiveDistri of each FOV','NumberTitle','off');
        set(gcf,'Position',[10 50 1200 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        
        for tempi=1:2
            if tempi==1
                y = tempBoolIndex123_or_num_eachFOV;
                temp_color = [0.25 0.25 0.25];
            elseif tempi==2
                y = selectiveCellBoolIndex_rProb_glm_delayx_num_eachFOV;
                temp_color = [0.6350 0.0780 0.1840];
            end
            
            plot(y,'LineWidth',1.5,'color',temp_color);
            hold on
        end
        
        legend('Memory','Meta',...
            'Location','northeast','fontsize',12)
        
        
        xlim([1-0.5 FOV_num+0.5]);
        ylim([0 max(y)]);
        
        % set(gca,'YScale','log');
        set(gca, 'FontSize', 14)
        set(gca,'box','off');% 取消右、上边框
        set(gca,'XTick',1:FOV_num);%设置要显示坐标刻度的范围
        xlabel('FOV', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('Roi num', 'FontSize', 15, 'FontWeight', 'bold');
        if if_selective_0delay1_1delay2_2all == 0
            title(sprintf('Selective roi distribution of each FOV, delay 1'),'fontsize',16);
        elseif if_selective_0delay1_1delay2_2all == 1
            title(sprintf('Selective roi distribution of each FOV, delay 2'),'fontsize',16);
        end
    end
    
    
    
    
    if if_plot_rProbTuning == 1
        %% Fig, rProb tuning
        fig = figure('Name','rProb tuning','NumberTitle','off');
        %set(gcf,'Position',[10 50 235 190*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 235*1.2 162*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 235*0.95 194*1.3]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        
        if_smoothHistogram = 0;
        
        nexttile
        
        temp_x_min = 0;
        temp_x_max = 0;
        
        x = r2_rProb_glm_delay1_multiFOV;
        
        if if_smoothHistogram == 0
            %h = histogram(x,'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.5,'edgeColor',[0.25 0.25 0.25]);
            h = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5,'edgeColor',[1 1 1]*0.5);
            
            [temp_x_min,temp_x_max] = bounds(x);
            
            [temp_y_min,temp_y_max] = bounds( h.BinCounts);
            
        elseif if_smoothHistogram == 1
            temp_data = x';
            n=100;
            n=2^ceil(log2(n)); % round up n to the next power of 2;
            [pdf2,xmesh2,bandwidth2] = ksdensity(temp_data,'NumPoints',n,'Function','pdf');
            %[pdf2,xmesh2,bandwidth2] = ksdensity(temp_data,'NumPoints',n,'Function','pdf','Support',[-0.01 5.01],'BoundaryCorrection', 'Reflection');
            
            temp_color = [1 1 1]*0.3;
            plot(xmesh2,pdf2,'LineWidth',1.5,'color',temp_color);
            
            temp_x_max = max(temp_x_max,max(xmesh2));
            
        end
        
        hold on
        
                
        set(gca,'linewidth',1.5)
        
        xlim([temp_x_min-(temp_x_max-temp_x_min)*0.04 temp_x_max+(temp_x_max-temp_x_min)*0.04]);
        ylim([temp_y_min temp_y_max+(temp_y_max-temp_y_min)*0.04]);
        
        set(gca,'YScale','log');
        set(gca,'YTick',[1 10 100 1000]);
        set(gca,'YMinorTick','off')
        %set(gca,'TickLength',[0.02, 0.1])
        
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('r2', 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel('pdf', 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel('Count', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Neuron count', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Offloading rate\nTuning distribution'),'fontsize',12);
        
        a = 1;
        
    end
    
    
    if if_plot_tuningComparison == 1
        %% Fig, precision and rProb tuning
        fig = figure('Name','tuning comparison','NumberTitle','off');
        set(gcf,'Position',[10 350 223*2 252*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        
        %% Pure and mix comparison (precision)
        nexttile                                
        x1 = r2_precision_glm_delayx_multiFOV(tempAB_pureB);
        h1 = histogram(x1,'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1,'edgeColor',[1 1 1]*0.25);
        hold on
        [temp_x1_min,temp_x1_max] = bounds(x1);
        [temp_y1_min,temp_y1_max] = bounds(h1.Values);
        
        x2 = r2_precision_glm_delayx_multiFOV(tempAB_and);
        h2 = histogram(x2,'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1,'edgeColor',[1 1 1]*0.75);
        hold on
        [temp_x2_min,temp_x2_max] = bounds(x2);
        [temp_y2_min,temp_y2_max] = bounds(h2.Values);
        
        [~,p_AB_pureMixB,~,~] = ttest2(x1,x2);        
        temp_p = p_AB_pureMixB;
        
        
        temp_x_min = min(temp_x1_min,temp_x2_min);
        temp_x_max = max(temp_x1_max,temp_x2_max);
        temp_y_min = min(temp_y1_min,temp_y2_min);
        temp_y_max = max(temp_y1_max,temp_y2_max);
        
        hold on
        
        
        le = legend('Pure','Mix',...
            'Location','northeast','fontsize',9);
        le.ItemTokenSize = ones(1,2)*10;
               
        
        tempTxt = sprintf('ns');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(mean([temp_x_min temp_x_max]),temp_y_max+(temp_y_max-temp_y_min)*0.02,...
            tempTxt,'Color','black','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
        
        
        xlim([temp_x_min-(temp_x_max-temp_x_min)*0.04 temp_x_max+(temp_x_max-temp_x_min)*0.04]);
        ylim([temp_y_min temp_y_max+(temp_y_max-temp_y_min)*0.1]);
        
        set(gca,'linewidth',1.5)                
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('r2', 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel('Neuron count', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Pdf', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Precision\nTuning distribution'),'fontsize',12);        
        
        
        %% Pure and mix comparison (offloading rate)
        nexttile                                
        x1 = r2_rProb_glm_delayx_multiFOV(tempAB_pureA);
        h1 = histogram(x1,'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1,'edgeColor',[1 1 1]*0.25);
        hold on
        [temp_x1_min,temp_x1_max] = bounds(x1);
        [temp_y1_min,temp_y1_max] = bounds(h1.Values);
        
        x2 = r2_rProb_glm_delayx_multiFOV(tempAB_and);
        h2 = histogram(x2,'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1,'edgeColor',[1 1 1]*0.75);
        hold on
        [temp_x2_min,temp_x2_max] = bounds(x2);
        [temp_y2_min,temp_y2_max] = bounds(h2.Values);
        
        [~,p_AB_pureMixA,~,~] = ttest2(x1,x2);        
        temp_p = p_AB_pureMixA;
        
        temp_x_min = min(temp_x1_min,temp_x2_min);
        temp_x_max = max(temp_x1_max,temp_x2_max);
        temp_y_min = min(temp_y1_min,temp_y2_min);
        temp_y_max = max(temp_y1_max,temp_y2_max);
        
        hold on
        
        
        le = legend('Pure','Mix',...
            'Location','northeast','fontsize',9);
        le.ItemTokenSize = ones(1,2)*10;
        
        
        tempTxt = sprintf('ns');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(mean([temp_x_min temp_x_max]),temp_y_max+(temp_y_max-temp_y_min)*0.02,...
            tempTxt,'Color','black','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
        
        
        xlim([temp_x_min-(temp_x_max-temp_x_min)*0.04 temp_x_max+(temp_x_max-temp_x_min)*0.04]);
        ylim([temp_y_min temp_y_max+(temp_y_max-temp_y_min)*0.1]);
        
        
        set(gca,'linewidth',1.5)                
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('r2', 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel('Neuron count', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Pdf', 'FontSize', 12, 'FontWeight', 'bold');        
        title(sprintf('Offloading rate\nTuning distribution'),'fontsize',12);
        
        
        
        %% Fig, location precision tuning
        fig = figure('Name','tuning comparison','NumberTitle','off');
        set(gcf,'Position',[10 650 223*2 252*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,5,'TileSpacing','Compact','Padding','Compact');
        
        %% Pure and mix comparison (location)
        nexttile([1 3])   
        
        tempL_index = find(tempL==true);                
        tempBL_pureL_index = find(tempBL_pureL==true);
        tempBL_and_index = find(tempBL_and==true);
        
        tempBL_pureL_index2 = ismember(tempL_index,tempBL_pureL_index);
        tempBL_and_index2 = ismember(tempL_index,tempBL_and_index);
        
        x1 = tempBL_pureL_index2;
        x2 = tempBL_and_index2;
        
        temp_1_pure = temp_beta_population_std{1}(x1);
        temp_2_pure = temp_beta_population_std{2}(x1);
        temp_3_pure = temp_beta_population_std{3}(x1);
                
        temp_1_mix = temp_beta_population_std{1}(x2);
        temp_2_mix = temp_beta_population_std{2}(x2);
        temp_3_mix = temp_beta_population_std{3}(x2);
                
        [~,p_length1_pureMix,~,~] = ttest2(temp_1_pure,temp_1_mix);
        [~,p_length2_pureMix,~,~] = ttest2(temp_2_pure,temp_2_mix);
        [~,p_length3_pureMix,~,~] = ttest2(temp_3_pure,temp_3_mix);
        
        temp_p1 = p_length1_pureMix;
        temp_p2 = p_length2_pureMix;
        temp_p3 = p_length3_pureMix;
        
        
        temp_y_min = 0;
        temp_y_max = max([temp_1_pure;temp_2_pure;temp_3_pure;temp_1_mix;temp_2_mix;temp_3_mix]);
        
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 6, 1);
        
        
        temp_data = [temp_1_pure;temp_1_mix;temp_2_pure;temp_2_mix;temp_3_pure;temp_3_mix];
        
        g1 = repmat({'Pure1'},length(temp_1_pure),1);
        g2 = repmat({'Mix1'},length(temp_1_mix),1);
        g3 = repmat({'Pure2'},length(temp_2_pure),1);
        g4 = repmat({'Mix2'},length(temp_2_mix),1);
        g5 = repmat({'Pure3'},length(temp_3_pure),1);
        g6 = repmat({'Mix3'},length(temp_3_mix),1);
        
        temp_label = [g1; g2; g3; g4; g5; g6];

        violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'Pure1'};{'Mix1'};{'Pure2'};{'Mix2'};{'Pure3'};{'Mix3'}]);
        
        hold on
        
        
        tempTxt = sprintf('ns');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center'); 
        
        tempTxt = sprintf('ns');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(3.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center'); 
        
        tempTxt = sprintf('ns');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(5.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');         
        
        set(gca,'linewidth',1.5)
        xlim([0.5 6.5])
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca, 'FontSize', 12)
        %set(gca,'XTickLabel', ["pure1"; "mix1";"pure2"; "mix2";"pure2"; "mix2";],'FontSize', 10);%给坐标加标签
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',10);
        xtickangle(0);
        set(gca,'box','off');% 取消右、上边框
        title(sprintf('Location\nTuning strength'),'fontsize',12);
                        
        
        %% Pure and mix comparison (precision)
        nexttile([1 2]) 

        x1 = r2_precision_glm_delayx_multiFOV(tempBL_pureB);
        x2 = r2_precision_glm_delayx_multiFOV(tempBL_and);
                        
        [~,p_BL_pureMix,~,~] = ttest2(x1,x2);
        
        temp_p = p_BL_pureMix;
        
        temp_y_min = 0;
        temp_y_max = max([x1;x2]);
        
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 6, 1);
        
        
        temp_data = [x1;x2];
        
        g1 = repmat({'Pure'},length(x1),1);
        g2 = repmat({'Mix'},length(x2),1);
        
        temp_label = [g1; g2];

        violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'Pure'};{'Mix'}]);
        
        hold on
        
        
        tempTxt = sprintf('ns');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        set(gca,'linewidth',1.5)
        xlim([0.5 2.5])
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca, 'FontSize', 12)
        %set(gca,'XTickLabel', ["pure"; "mix"],'FontSize', 10);%给坐标加标签
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',10);        
        xtickangle(0);
        set(gca,'box','off');% 取消右、上边框
        title(sprintf('Precision\nTuning distribution'),'fontsize',12);
        
        
        
        
        
        
        %% Fig, location offloading tuning
        fig = figure('Name','tuning comparison','NumberTitle','off');
        set(gcf,'Position',[500 650 223*2 252*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,5,'TileSpacing','Compact','Padding','Compact');
        
        %% Pure and mix comparison (location)
        nexttile([1 3])   
        
        tempL_index = find(tempL==true);                
        tempAL_pureL_index = find(tempAL_pureL==true);
        tempAL_and_index = find(tempAL_and==true);
        
        tempAL_pureL_index2 = ismember(tempL_index,tempAL_pureL_index);
        tempAL_and_index2 = ismember(tempL_index,tempAL_and_index);
        
        x1 = tempAL_pureL_index2;
        x2 = tempAL_and_index2;
        
        temp_1_pure = temp_beta_population_std{1}(x1);
        temp_2_pure = temp_beta_population_std{2}(x1);
        temp_3_pure = temp_beta_population_std{3}(x1);
                
        temp_1_mix = temp_beta_population_std{1}(x2);
        temp_2_mix = temp_beta_population_std{2}(x2);
        temp_3_mix = temp_beta_population_std{3}(x2);
                
        [~,p_length1_pureMix,~,~] = ttest2(temp_1_pure,temp_1_mix);
        [~,p_length2_pureMix,~,~] = ttest2(temp_2_pure,temp_2_mix);
        [~,p_length3_pureMix,~,~] = ttest2(temp_3_pure,temp_3_mix);
        
        temp_p1 = p_length1_pureMix;
        temp_p2 = p_length2_pureMix;
        temp_p3 = p_length3_pureMix;
        
        
        temp_y_min = 0;
        temp_y_max = max([temp_1_pure;temp_2_pure;temp_3_pure;temp_1_mix;temp_2_mix;temp_3_mix]);
        
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 6, 1);
        
        
        temp_data = [temp_1_pure;temp_1_mix;temp_2_pure;temp_2_mix;temp_3_pure;temp_3_mix];
        
        g1 = repmat({'Pure1'},length(temp_1_pure),1);
        g2 = repmat({'Mix1'},length(temp_1_mix),1);
        g3 = repmat({'Pure2'},length(temp_2_pure),1);
        g4 = repmat({'Mix2'},length(temp_2_mix),1);
        g5 = repmat({'Pure3'},length(temp_3_pure),1);
        g6 = repmat({'Mix3'},length(temp_3_mix),1);
        
        temp_label = [g1; g2; g3; g4; g5; g6];

        violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'Pure1'};{'Mix1'};{'Pure2'};{'Mix2'};{'Pure3'};{'Mix3'}]);
        
        hold on
        
        
        tempTxt = sprintf('ns');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center'); 
        
        tempTxt = sprintf('ns');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(3.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center'); 
        
        tempTxt = sprintf('ns');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(5.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');         
        
        set(gca,'linewidth',1.5)
        xlim([0.5 6.5])
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca, 'FontSize', 12)
        %set(gca,'XTickLabel', ["pure1"; "mix1";"pure2"; "mix2";"pure2"; "mix2";],'FontSize', 10);%给坐标加标签
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',10);
        xtickangle(0);
        set(gca,'box','off');% 取消右、上边框
        title(sprintf('Location\nTuning strength'),'fontsize',12);
                        
        
        %% Pure and mix comparison (offloading)
        nexttile([1 2]) 

        x1 = r2_rProb_glm_delayx_multiFOV(tempAL_pureA);
        x2 = r2_rProb_glm_delayx_multiFOV(tempAL_and);
                        
        [~,p_AL_pureMix,~,~] = ttest2(x1,x2);
        
        temp_p = p_AL_pureMix;
        
        temp_y_min = 0;
        temp_y_max = max([x1;x2]);
        
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 6, 1);
        
        
        temp_data = [x1;x2];
        
        g1 = repmat({'Pure'},length(x1),1);
        g2 = repmat({'Mix'},length(x2),1);
        
        temp_label = [g1; g2];

        violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'Pure'};{'Mix'}]);
        
        hold on
        
        
        tempTxt = sprintf('ns');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,(temp_y_max-temp_y_min)*1.03,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        set(gca,'linewidth',1.5)
        xlim([0.5 2.5])
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.1]);
        set(gca, 'FontSize', 12)
        %set(gca,'XTickLabel', ["pure"; "mix"],'FontSize', 10);%给坐标加标签
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',10);        
        xtickangle(0);
        set(gca,'box','off');% 取消右、上边框
        title(sprintf('Offloading rate\nTuning distribution'),'fontsize',12);
        
    end
    
    
end


if false
   %[M,I]=max(r2_rProb_glm_delayx_multiFOV(tempAB_pureA)); 
   [M,I]=sort(r2_rProb_glm_delayx_multiFOV(tempAB_pureA),'descend');

   temp1 = 19;%range from 1 to 35. 6,10(bad),14(bad),17,19,20(bad)
   
   tempIndex = find(tempAB_pureA==true);
   tempIndexB = tempIndex(I(temp1));
   a1 = r2_rProb_glm_delayx_multiFOV(tempIndexB);
   a2 = r2_precision_glm_delayx_multiFOV(tempIndexB);
   fprintf('rProb=%.4f, precision=%.4f.\n',a1,a2);  
end


if false
   [M,I]=sort(r2_precision_glm_delayx_multiFOV(tempBL_and),'descend');

   temp1 = 16;%1(miss),2(miss),3(conflict),7,12
   
   tempIndex = find(tempBL_and==true);
   tempIndexB = tempIndex(I(temp1));
   a2 = r2_precision_glm_delayx_multiFOV(tempIndexB);
   fprintf('precision=%.4f.\n',a2);  
   tempIndexB
end


if false
   [M,I]=sort(r2_precision_glm_delayx_multiFOV(tempBL_pureB),'descend');

   temp1 = 7;%1(miss)
   
   tempIndex = find(tempBL_pureB==true);
   tempIndexB = tempIndex(I(temp1));
   a2 = r2_precision_glm_delayx_multiFOV(tempIndexB);
   fprintf('precision=%.4f.\n',a2);  
   tempIndexB
end


% if false
if true
   glm_beta_lengthx_delay1Bin_multiFOV;
   glm_r2_lengthx_delay1Bin_multiFOV;
   
   temp1 = glm_beta_lengthx_delay1Bin_multiFOV(:,1:6);
   temp2 = glm_beta_lengthx_delay1Bin_multiFOV(:,7:12);
   temp3 = glm_beta_lengthx_delay1Bin_multiFOV(:,13:18);   
   glm_beta_6locMean_delay1Bin_multiFOV = (temp1+temp2+temp3)/3;
    
   temp1 = glm_r2_lengthx_delay1Bin_multiFOV(:,1);
   temp2 = glm_r2_lengthx_delay1Bin_multiFOV(:,2);
   temp3 = glm_r2_lengthx_delay1Bin_multiFOV(:,3);   
   glm_r2_6locMean_delay1Bin_multiFOV = (temp1+temp2+temp3)/3;
   

end