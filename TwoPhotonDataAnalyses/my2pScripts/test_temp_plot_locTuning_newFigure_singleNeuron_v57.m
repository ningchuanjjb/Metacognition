% Chuan's 23th script (20251214)
% This script: To plot single-neuron's location memory selectivity, related to figure 2.

% please run 'stepF4A_memoryMetaMismatch_singleNeuron.m' first

%% Initialization for demo data

% if_loadDemoData = if_loadDemoData;
if_loadDemoData = 1;

loadExampleNeuronID_demo = 1;%1 or 2;1 for Figure 1G, 2 for Figure 1H.


if if_loadDemoData == 1
    if loadExampleNeuronID_demo == 1
        temp_F_dff_exampleNeuron_memoryX = decodingDataSimplified.F_dff_exampleNeuron_memory1;
    elseif loadExampleNeuronID_demo == 2
        temp_F_dff_exampleNeuron_memoryX = decodingDataSimplified.F_dff_exampleNeuron_memory2;
    end
end

%% Initialization
if_plot = 1;

if_plot_multiPanel = 0;
if_plot_singlePanel = 0;
if_plot_multiPanelB = 0;
if_plot_multiPanelB2 = 0;% plot length 1
if_plot_multiPanelB3 = 0;% plot mean of length 1 2 3, from lastT (forward)
if_plot_multiPanelB3_2 = 1;% plot mean of length 1 2 3, from lastT (forward), but trial sorted as seq-A VS. non-seqA
if_plot_multiPanelB4 = 0;% plot mean of length 1 2 3 (for Fig. 4)
if_plot_multiPanelC = 0;% plot activity of each seq 


if if_loadDemoData == 1
    if loadExampleNeuronID_demo == 1
        if_plot_multiPanelB3 = 1;
        if_plot_multiPanelB3_2 = 0;
    elseif loadExampleNeuronID_demo == 2
        if_plot_multiPanelB3 = 0;
        if_plot_multiPanelB3_2 = 1;
    end
end


if_delay1_forward0_backward1 = 0; % new

if if_plot_multiPanelB3 == 1
    if_delay1_forward0_backward1 = 0;
end

if if_plot_multiPanelB4 == 1
    if_delay1_forward0_backward1 = 1;
    %if_delay1_forward0_backward1 = 0;
end


if_colormap_loadEnhanced = 0;

if_seqIdentity_default0_tunedLoc1 = 1;
% if_seqIdentity_once0_multi1 = 0;

% high location, high meta
% cellIndex_suite2p_temptempRaw = 92;%
% cellIndex_suite2p_temptempRaw = 25;%
% cellIndex_suite2p_temptempRaw = 20;% good
% cellIndex_suite2p_temptempRaw = 124;% (0510A) good
% cellIndex_suite2p_temptempRaw = 1444;% (0510A)



% high location, low meta
% cellIndex_suite2p_temptempRaw = 2;% good
% cellIndex_suite2p_temptempRaw = 11;%
% cellIndex_suite2p_temptempRaw = 197;% (0510A)
% cellIndex_suite2p_temptempRaw = 50;% (0614A)
% cellIndex_suite2p_temptempRaw = 17;% (0312A_Z)
% cellIndex_suite2p_temptempRaw = 19;% (0312A_Z)
% cellIndex_suite2p_temptempRaw = 101;% (0312A_Z)
% cellIndex_suite2p_temptempRaw = 29;% (0124A_Z)


% low location, high meta
% cellIndex_suite2p_temptempRaw = 724;%
% cellIndex_suite2p_temptempRaw = 695;% good
% cellIndex_suite2p_temptempRaw = 816;%


% suite2p cellIndex in sessionB
% cellIndex_suite2p_temptempRaw = 19;% loc5
% cellIndex_suite2p_temptempRaw = 11;% loc3

cellIndex_suite2p_temptempRaw = 25;% all loc
% cellIndex_suite2p_temptempRaw = 41;% loc6

% cellIndex_suite2p_temptempRaw = 186;% loc6, (0528A)
% cellIndex_suite2p_temptempRaw = 121;% loc3, (0528A)
% cellIndex_suite2p_temptempRaw = 886;% loc3, (0528A)

% cellIndex_suite2p_temptempRaw = 41;% loc4, (0522A)
% cellIndex_suite2p_temptempRaw = 46;% loc3, (0522A)
% cellIndex_suite2p_temptempRaw = 346;% loc3, (0522A)

% cellIndex_suite2p_temptempRaw = 9;% loc1, (0527A)
% cellIndex_suite2p_temptempRaw = 138;% loc3, (0527A)

% cellIndex_suite2p_temptempRaw = 8;% loc2, (0504A)
% cellIndex_suite2p_temptempRaw = 6;% loc2, (0526A)
% cellIndex_suite2p_temptempRaw = 1;% loc2, (0502A)
% cellIndex_suite2p_temptempRaw = 25;% loc3, (0605A)

% pure precision tuning neuron
% cellIndex_suite2p_temptempRaw = 656;%(0513A)
% cellIndex_suite2p_temptempRaw = 21;%(0509A)

% mix tuning for location and precision neuron
% cellIndex_suite2p_temptempRaw = 461;%(0512A)
% cellIndex_suite2p_temptempRaw = 428;%(0513A)
% cellIndex_suite2p_temptempRaw = 429;%(0513A)
% cellIndex_suite2p_temptempRaw = 0;%(0513A)
% cellIndex_suite2p_temptempRaw = 27;%(0526A)
% cellIndex_suite2p_temptempRaw = 5;%(0513A)

if if_loadDemoData == 1
    if loadExampleNeuronID_demo == 1
        cellIndex_suite2p_temptempRaw = 19;% loc5
    elseif loadExampleNeuronID_demo == 2
        cellIndex_suite2p_temptempRaw = 25;% all loc
    end
end




TunedLoc = [];
if cellIndex_suite2p_temptempRaw == 19
    TunedLoc = 5;
elseif cellIndex_suite2p_temptempRaw == 25
    TunedLoc = [1 2 6];
elseif cellIndex_suite2p_temptempRaw == 9
    TunedLoc = 1;
elseif cellIndex_suite2p_temptempRaw == 121
    TunedLoc = 3;
elseif cellIndex_suite2p_temptempRaw == 41
    TunedLoc = 4;
elseif cellIndex_suite2p_temptempRaw == 8
    TunedLoc = 2;
elseif cellIndex_suite2p_temptempRaw == 6
    TunedLoc = 2;        
end

if isempty(TunedLoc)
    if_seqIdentity_default0_tunedLoc1 = 0;
end
    
% multi_rgbColor = ...
%     [228,26,28;
%     55,126,184;
%     77,175,74;
%     152,78,163;
%     255,127,0;
%     255,255,51]/255;
multi_rgbColor = ...
    [228,26,28;
    55,126,184;
    77,175,74;
    152,78,163;
    255,127,0;
    [255,255,51].*0.8]/255;

% temp_EdgeColor = [1 1 1]*0.725;%0.62,0.725
temp_EdgeColor = 'none';%[0.62,0.62,0.62]-->'none'
temp_LineWidth = 1;%1
temp_FaceAlpha = 0.35;%0.1-->0.05-->0.35

decodingDataSimplified = decodingDataSimplified; %#ok<*ASGSL>

if if_dff_singleSession1_twoSession2_allSession3 == 2
    temp_cellIndex_target = find(decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end)==cellIndex_suite2p_temptempRaw);
    cellIndex_suite2p_temptemp = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(temp_cellIndex_target,1);
elseif if_dff_singleSession1_twoSession2_allSession3 == 1
    cellIndex_suite2p_temptemp = cellIndex_suite2p_temptempRaw;
end

xlim_padSize1 = 10;%2,5.5
% periodSkipInterval = -0.5;%3-->2-->0,-0.5
periodSkipInterval = 0;%3-->2-->0,-0.5
% temp_ticks_baselinePeriod = baselinePeriod_interval;
% temp_ticks_decisionPeriodA = temp_ticks_baselinePeriod(end)+periodSkipInterval+decisionPeriodA_interval;
% temp_ticks_decisionPeriod = temp_ticks_decisionPeriodA(end)+periodSkipInterval+decisionPeriod_interval;

% backgounrdColor = [1 1 1]*0.825;%0.875
backgounrdColor = [1 1 1];

temp_id2 = find(cellIndex_suite2p == cellIndex_suite2p_temptemp);

r2_6loc(temp_id2);
p_6loc;
fprintf('r2_6loc = %.3f, p_6loc = %.3f.\n',r2_6loc(temp_id2),p_6loc(temp_id2));


%% cellID_F_dff_decisionBin1_seqMerged
if if_loadDemoData == 0
    F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
elseif if_loadDemoData == 1
    F_dff_decisionBin1 = decodingDataSimplified.F_dff_decisionBin_demo;
end
cellID_F_dff_decisionBin1 = F_dff_decisionBin1(temp_id2,:);
cellID_F_dff_decisionBin1_seqMerged = zeros(1,sum(numSeq));
for tempi=1:sum(numSeq)
    temptempBoolIndex = seqIndex==tempi & trialIndex_bool_memoryCorrect;    
    
    temp_dff = cellID_F_dff_decisionBin1(temptempBoolIndex);
    cellID_F_dff_decisionBin1_seqMerged(tempi) = mean(temp_dff);
    %if tempi == 41
    %   a = 1; 
    %end
end

cellID_F_dff_decisionBin1_seqMerged;

cellID_F_dff_decisionBin1_seqMerged_eachLoc = zeros(numFrames,16);
temp_seqIndex_eachLoc = zeros(numFrames,16);
for tempi=1:numFrames
    tempCount = 0;
    for tempj=1:sum(numSeq(1:3))
        if boolIndex_location_seq(tempi,tempj) == true
            tempCount = tempCount+1;
            cellID_F_dff_decisionBin1_seqMerged_eachLoc(tempi,tempCount) = cellID_F_dff_decisionBin1_seqMerged(tempj);
            temp_seqIndex_eachLoc(tempi,tempCount) = tempj;
        end
    end
end


%% cellID_F_dff_length3_sample_location
if if_loadDemoData == 0
    cellID_F_dff_length1_sample = squeeze(F_dff_length1_sample(temp_id2,:,:));
    cellID_F_dff_length2_sample = squeeze(F_dff_length2_sample(temp_id2,:,:));
    cellID_F_dff_length3_sample = squeeze(F_dff_length3_sample(temp_id2,:,:));
elseif if_loadDemoData == 1
    cellID_F_dff_length1_sample = temp_F_dff_exampleNeuron_memoryX.F_dff_length1_sample_exampleNeuron_memory;
    cellID_F_dff_length2_sample = temp_F_dff_exampleNeuron_memoryX.F_dff_length2_sample_exampleNeuron_memory;
    cellID_F_dff_length3_sample = temp_F_dff_exampleNeuron_memoryX.F_dff_length3_sample_exampleNeuron_memory;
end
% cellID_F_dff_lastT_sample = squeeze(F_dff_length1_sample(temp_id2,:,:));

boolIndex_location_allTrial;
for tempi=1:3
    temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
    
    if tempi==1
        cellID_F_dff_lengthx_sample = cellID_F_dff_length1_sample;
    elseif tempi==2
        cellID_F_dff_lengthx_sample = cellID_F_dff_length2_sample;
    elseif tempi==3
        cellID_F_dff_lengthx_sample = cellID_F_dff_length3_sample;
    end
    
    cellID_F_dff_lengthx_sample_location = cell(numFrames,1);
    
    % temptempBoolIndex1 is length label
    % temptempBoolIndex2 is location label
    % temptempBoolIndex3 is length & location & memory correct label
    temptempBoolIndex1 = ismember(seqIndex,temp_range);
    for tempj=1:numFrames
        temptempBoolIndex2 = boolIndex_location_allTrial(tempj,:)==true;
        temptempBoolIndex3 = temptempBoolIndex1 & temptempBoolIndex2 & trialIndex_bool_memoryCorrect;
        
        %temptempBoolIndex3(1:1035) = false;
        
        cellID_F_dff_lengthx_sample_location{tempj} = cellID_F_dff_lengthx_sample(temptempBoolIndex3,:);        
    end
    
    if tempi==1
        cellID_F_dff_length1_sample_location = cellID_F_dff_lengthx_sample_location;
    elseif tempi==2
        cellID_F_dff_length2_sample_location = cellID_F_dff_lengthx_sample_location;
    elseif tempi==3
        cellID_F_dff_length3_sample_location = cellID_F_dff_lengthx_sample_location;
    end
end

cellID_F_dff_lastT_sample_location = cell(numFrames,1);
for tempi=1:numFrames
    temp1 = cellID_F_dff_length1_sample_location{tempi};
    temp2 = cellID_F_dff_length2_sample_location{tempi};
    temp3 = cellID_F_dff_length3_sample_location{tempi};
    
    tempLastT = [temp1;temp2(:,end-17:end);temp3(:,end-17:end)];
    cellID_F_dff_lastT_sample_location{tempi} = tempLastT;
end

%% cellID_F_dff_lastT_sample
cellID_F_dff_lastT_sample = nan(size(cellID_F_dff_length1_sample));

boolIndex_location_allTrial;
for tempi=1:3
    temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
    
    if tempi==1
        cellID_F_dff_lengthx_sample = cellID_F_dff_length1_sample;
    elseif tempi==2
        cellID_F_dff_lengthx_sample = cellID_F_dff_length2_sample;
    elseif tempi==3
        cellID_F_dff_lengthx_sample = cellID_F_dff_length3_sample;
    end

    temptempBoolIndex = ismember(seqIndex,temp_range);
    
    temp1 = cellID_F_dff_lengthx_sample(temptempBoolIndex,end-17:end);
    
    cellID_F_dff_lastT_sample(temptempBoolIndex,:) = temp1;    
end





%% cellID_F_dff_lengthx_decisionPeriodA_location
if if_loadDemoData == 0
    if if_delay1_forward0_backward1 == 1
        cellID_F_dff_length1_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,:));
        cellID_F_dff_length2_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,:));
        cellID_F_dff_length3_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,:));
    elseif if_delay1_forward0_backward1 == 0
        cellID_F_dff_length1_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,:));
        cellID_F_dff_length2_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,:));
        cellID_F_dff_length3_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,:));
        
        cellID_F_dff_length1_decisionPeriodA(:,1:decisionPeriodB_interval(2)) = squeeze(F_dff_decisionPeriodB(temp_id2,:,:));
        cellID_F_dff_length2_decisionPeriodA(:,1:decisionPeriodB_interval(2)) = squeeze(F_dff_decisionPeriodB(temp_id2,:,:));
        cellID_F_dff_length3_decisionPeriodA(:,1:decisionPeriodB_interval(2)) = squeeze(F_dff_decisionPeriodB(temp_id2,:,:));
    end
elseif if_loadDemoData == 1
    cellID_F_dff_length1_decisionPeriodA = temp_F_dff_exampleNeuron_memoryX.F_dff_decisionPeriodA_exampleNeuron_memory;
    cellID_F_dff_length2_decisionPeriodA = temp_F_dff_exampleNeuron_memoryX.F_dff_decisionPeriodA_exampleNeuron_memory;
    cellID_F_dff_length3_decisionPeriodA = temp_F_dff_exampleNeuron_memoryX.F_dff_decisionPeriodA_exampleNeuron_memory;
end

decisionPeriodA_interval = decisionPeriodA_interval;

% cellID_F_dff_length1_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,1:decisionPeriodA_interval(2)));
% cellID_F_dff_length2_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,1:decisionPeriodA_interval(2)));
% cellID_F_dff_length3_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,1:decisionPeriodA_interval(2)));

boolIndex_location_allTrial;
for tempi=1:3
    temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
    
    if tempi==1
        cellID_F_dff_lengthx_decisionPeriodA = cellID_F_dff_length1_decisionPeriodA;
    elseif tempi==2
        cellID_F_dff_lengthx_decisionPeriodA = cellID_F_dff_length2_decisionPeriodA;
    elseif tempi==3
        cellID_F_dff_lengthx_decisionPeriodA = cellID_F_dff_length3_decisionPeriodA;
    end
    
    cellID_F_dff_lengthx_decisionPeriodA_location = cell(numFrames,1);
    
    % temptempBoolIndex1 is length label
    % temptempBoolIndex2 is location label
    % temptempBoolIndex3 is length & location & memory correct label
    temptempBoolIndex1 = ismember(seqIndex,temp_range);
    for tempj=1:numFrames
        temptempBoolIndex2 = boolIndex_location_allTrial(tempj,:)==true;
        temptempBoolIndex3 = temptempBoolIndex1 & temptempBoolIndex2 & trialIndex_bool_memoryCorrect;
        
        %temptempBoolIndex3(1:1035) = false;
        
        cellID_F_dff_lengthx_decisionPeriodA_location{tempj} = cellID_F_dff_lengthx_decisionPeriodA(temptempBoolIndex3,:);
    end
    
    if tempi==1
        cellID_F_dff_length1_decisionPeriodA_location = cellID_F_dff_lengthx_decisionPeriodA_location;
    elseif tempi==2
        cellID_F_dff_length2_decisionPeriodA_location = cellID_F_dff_lengthx_decisionPeriodA_location;
    elseif tempi==3
        cellID_F_dff_length3_decisionPeriodA_location = cellID_F_dff_lengthx_decisionPeriodA_location;
    end
end

cellID_F_dff_lastT_decisionPeriodA_location = cell(numFrames,1);
for tempi=1:numFrames
    temp1 = cellID_F_dff_length1_decisionPeriodA_location{tempi};
    temp2 = cellID_F_dff_length2_decisionPeriodA_location{tempi};
    temp3 = cellID_F_dff_length3_decisionPeriodA_location{tempi};
    
    tempLastT = [temp1;temp2;temp3];
    cellID_F_dff_lastT_decisionPeriodA_location{tempi} = tempLastT;
end

%%
for plot_lengthFlag=1:4
    if plot_lengthFlag == 1
        cellID_F_dff_lengthx_sample_location = cellID_F_dff_length1_sample_location;
        cellID_F_dff_lengthx_decisionPeriodA_location = cellID_F_dff_length1_decisionPeriodA_location;        
        lengthx_sample_interval = length1_sample_interval;
    elseif plot_lengthFlag == 2
        cellID_F_dff_lengthx_sample_location = cellID_F_dff_length2_sample_location;
        cellID_F_dff_lengthx_decisionPeriodA_location = cellID_F_dff_length2_decisionPeriodA_location;                
        lengthx_sample_interval = length2_sample_interval;
    elseif plot_lengthFlag == 3
        cellID_F_dff_lengthx_sample_location = cellID_F_dff_length3_sample_location;
        cellID_F_dff_lengthx_decisionPeriodA_location = cellID_F_dff_length3_decisionPeriodA_location;                
        lengthx_sample_interval = length3_sample_interval;
    elseif plot_lengthFlag == 4 % lastT
        cellID_F_dff_lengthx_sample_location = cellID_F_dff_lastT_sample_location;
        cellID_F_dff_lengthx_decisionPeriodA_location = cellID_F_dff_lastT_decisionPeriodA_location;                
        lengthx_sample_interval = length1_sample_interval;        
    end
    
    cellID_F_dff_lengthx_sample_location_mean = nan(numFrames,lengthx_sample_interval(end));
    cellID_F_dff_lengthx_sample_location_sem = nan(numFrames,lengthx_sample_interval(end));
    for tempi=1:numFrames
        temp_dff = cellID_F_dff_lengthx_sample_location{tempi};
        cellID_F_dff_lengthx_sample_location_mean(tempi,:) = mean(temp_dff,1,'omitnan');
        cellID_F_dff_lengthx_sample_location_sem(tempi,:) = std(temp_dff,1,1,'omitnan')./sqrt(size(temp_dff,1));
    end
    
    cellID_F_dff_lengthx_decisionPeriodA_location_mean = nan(numFrames,decisionPeriodA_interval(end));
    cellID_F_dff_lengthx_decisionPeriodA_location_sem = nan(numFrames,decisionPeriodA_interval(end));
    %cellID_F_dff_lengthx_decisionPeriodA_location_mean = nan(numFrames,decisionPeriodA_interval(2));
    %cellID_F_dff_lengthx_decisionPeriodA_location_sem = nan(numFrames,decisionPeriodA_interval(2));    
    for tempi=1:numFrames
        temp_dff = cellID_F_dff_lengthx_decisionPeriodA_location{tempi};
        cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,:) = mean(temp_dff,1,'omitnan');
        cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,:) = std(temp_dff,1,1,'omitnan')./sqrt(size(temp_dff,1));
    end    
    
    if plot_lengthFlag == 1
        cellID_F_dff_length1_sample_location_mean = cellID_F_dff_lengthx_sample_location_mean;
        cellID_F_dff_length1_sample_location_sem = cellID_F_dff_lengthx_sample_location_sem;
        cellID_F_dff_length1_decisionPeriodA_location_mean = cellID_F_dff_lengthx_decisionPeriodA_location_mean;
        cellID_F_dff_length1_decisionPeriodA_location_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem;        
    elseif plot_lengthFlag == 2
        cellID_F_dff_length2_sample_location_mean = cellID_F_dff_lengthx_sample_location_mean;
        cellID_F_dff_length2_sample_location_sem = cellID_F_dff_lengthx_sample_location_sem;
        cellID_F_dff_length2_decisionPeriodA_location_mean = cellID_F_dff_lengthx_decisionPeriodA_location_mean;
        cellID_F_dff_length2_decisionPeriodA_location_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem;                
    elseif plot_lengthFlag == 3
        cellID_F_dff_length3_sample_location_mean = cellID_F_dff_lengthx_sample_location_mean;
        cellID_F_dff_length3_sample_location_sem = cellID_F_dff_lengthx_sample_location_sem;
        cellID_F_dff_length3_decisionPeriodA_location_mean = cellID_F_dff_lengthx_decisionPeriodA_location_mean;
        cellID_F_dff_length3_decisionPeriodA_location_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem;                
    elseif plot_lengthFlag == 4
        cellID_F_dff_lastT_sample_location_mean = cellID_F_dff_lengthx_sample_location_mean;
        cellID_F_dff_lastT_sample_location_sem = cellID_F_dff_lengthx_sample_location_sem;
        cellID_F_dff_lastT_decisionPeriodA_location_mean = cellID_F_dff_lengthx_decisionPeriodA_location_mean;
        cellID_F_dff_lastT_decisionPeriodA_location_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem;                        
    end
    
end
a = 1;

% %% Merge sample and decisionPeriodA
temp1_sample = min(cellID_F_dff_length1_sample_location_mean-cellID_F_dff_length1_sample_location_sem,[],'all');
temp2_sample = min(cellID_F_dff_length2_sample_location_mean-cellID_F_dff_length2_sample_location_sem,[],'all');
temp3_sample = min(cellID_F_dff_length3_sample_location_mean-cellID_F_dff_length3_sample_location_sem,[],'all');
min_sample = min([temp1_sample,temp2_sample,temp3_sample]);

temp1_decisionPeriodA = min(cellID_F_dff_length1_decisionPeriodA_location_mean-cellID_F_dff_length1_decisionPeriodA_location_sem,[],'all');
temp2_decisionPeriodA = min(cellID_F_dff_length2_decisionPeriodA_location_mean-cellID_F_dff_length2_decisionPeriodA_location_sem,[],'all');
temp3_decisionPeriodA = min(cellID_F_dff_length3_decisionPeriodA_location_mean-cellID_F_dff_length3_decisionPeriodA_location_sem,[],'all');
min_decisionPeriodA = min([temp1_decisionPeriodA,temp2_decisionPeriodA,temp3_decisionPeriodA]);

min_sample_decisionPeriodA = min([min_sample,min_decisionPeriodA]);

temp1_sample = max(cellID_F_dff_length1_sample_location_mean+cellID_F_dff_length1_sample_location_sem,[],'all');
temp2_sample = max(cellID_F_dff_length2_sample_location_mean+cellID_F_dff_length2_sample_location_sem,[],'all');
temp3_sample = max(cellID_F_dff_length3_sample_location_mean+cellID_F_dff_length3_sample_location_sem,[],'all');
max_sample = max([temp1_sample,temp2_sample,temp3_sample]);

temp1_decisionPeriodA = max(cellID_F_dff_length1_decisionPeriodA_location_mean+cellID_F_dff_length1_decisionPeriodA_location_sem,[],'all');
temp2_decisionPeriodA = max(cellID_F_dff_length2_decisionPeriodA_location_mean+cellID_F_dff_length2_decisionPeriodA_location_sem,[],'all');
temp3_decisionPeriodA = max(cellID_F_dff_length3_decisionPeriodA_location_mean+cellID_F_dff_length3_decisionPeriodA_location_sem,[],'all');
max_decisionPeriodA = max([temp1_decisionPeriodA,temp2_decisionPeriodA,temp3_decisionPeriodA]);

max_sample_decisionPeriodA = max([max_sample,max_decisionPeriodA]);

if cellIndex_suite2p_temptempRaw == 11 && if_plot_multiPanelB == 1
    max_sample_decisionPeriodA = 1.6;
end



%% cellID_F_dff_decisionPeriodA_seqSubset

seqSubset_locBoolIndex = [1,1,0,0,0,1]==1;

if if_loadDemoData == 0
    if if_delay1_forward0_backward1 == 1
        cellID_F_dff_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,:));
    elseif if_delay1_forward0_backward1 == 0
        cellID_F_dff_decisionPeriodA = squeeze(F_dff_decisionPeriodA(temp_id2,:,:));
        cellID_F_dff_decisionPeriodA(:,1:decisionPeriodB_interval(2)) = squeeze(F_dff_decisionPeriodB(temp_id2,:,:));
    end
elseif if_loadDemoData == 1
    cellID_F_dff_decisionPeriodA = temp_F_dff_exampleNeuron_memoryX.F_dff_decisionPeriodA_exampleNeuron_memory;
end
% cellID_F_dff_decisionPeriodA_seqSubset

boolIndex_location_allTrial;
boolIndex_location_allTrial_T = boolIndex_location_allTrial';

seqSubset_trialBoolIndex = nan(size(boolIndex_location_allTrial_T,1),1);

% seqSubset_trialBoolIndex =
temp1 = boolIndex_location_allTrial_T & seqSubset_locBoolIndex;
temp2 = sum(temp1,2)>0;% hit
% seqSubset_trialBoolIndex = temp2;

temp3 = boolIndex_location_allTrial_T & ~seqSubset_locBoolIndex;
temp4 = sum(temp3,2)>0;% false alarm

temp24 = temp2 & ~temp4;
seqSubset_trialBoolIndex = temp24;


temptempBoolIndex = seqSubset_trialBoolIndex & (seqIndex' <= 41) & trialIndex_bool_memoryCorrect';
cellID_F_dff_decisionPeriodA_seqSubset = cellID_F_dff_decisionPeriodA(temptempBoolIndex,:);

temptempBoolIndex = (~seqSubset_trialBoolIndex) & (seqIndex' <= 41) & trialIndex_bool_memoryCorrect';
cellID_F_dff_decisionPeriodA_nonSeqSubset = cellID_F_dff_decisionPeriodA(temptempBoolIndex,:);


%% cellID_F_dff_lastT_sample_seqSubset

temptempBoolIndex = seqSubset_trialBoolIndex & (seqIndex' <= 41) & trialIndex_bool_memoryCorrect';
cellID_F_dff_lastT_sample_seqSubset = cellID_F_dff_lastT_sample(temptempBoolIndex,:);

temptempBoolIndex = (~seqSubset_trialBoolIndex) & (seqIndex' <= 41) & trialIndex_bool_memoryCorrect';
cellID_F_dff_lastT_sample_nonSeqSubset = cellID_F_dff_lastT_sample(temptempBoolIndex,:);


%% Plot
if if_plot == 1
    close all
    
    if if_plot_multiPanel == 1
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 520*0.97 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[35+0 42+0 297*1.05 229]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point        
        %set(gcf,'Position',[35+0 42+0 260*1.05 216*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                
        %set(gcf,'Position',[35+0 42+0 512 380*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+0 512*1.28 342*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 512*1.28*0.6*0.8*1.1*0.96 342*0.94*0.85*0.95*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        
        
        %t = tiledlayout(3,1,'TileSpacing','tight','Padding','compact'); %#ok<*NASGU>
        t = tiledlayout(3,1,'TileSpacing','tight','Padding','tight'); %#ok<*NASGU>        
        
        %t.Title.String = sprintf('suite2pID%d,%s',cellIndex_suite2p_temptemp,FOVName_currentFOV2);
        %t.Title.String = sprintf('suite2pBID%d,%s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2);        
        t.Title.String = sprintf(' suite2pBID%d,%s \n Location tuning in length 1 2 3',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2);                
        t.Title.FontSize = 10;
        t.Title.Interpreter = 'none';
        
        
        for plot_lengthFlag=1:3
        %if true
        %    plot_lengthFlag = 1;
            
            nexttile
            
            if plot_lengthFlag == 1
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length1_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length1_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length1_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length1_decisionPeriodA_location_sem;                
                lengthx_sample_interval = length1_sample_interval;
            elseif plot_lengthFlag == 2
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length2_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length2_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length2_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length2_decisionPeriodA_location_sem;                                
                lengthx_sample_interval = length2_sample_interval;
            elseif plot_lengthFlag == 3
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length3_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length3_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length3_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length3_decisionPeriodA_location_sem;                                
                lengthx_sample_interval = length3_sample_interval;
            end
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            %temp_ticks_decisionPeriodA = temp_ticks_sample(end)-1+decisionPeriodA_interval;
            
            
            %% Plot sample
            h_line = [];
            for tempi=1:numFrames
                x = temp_ticks_sample(1):temp_ticks_sample(end);
                y = cellID_F_dff_lengthx_sample_location_mean(tempi,:);
                %h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                h_line = [h_line plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>                
                hold on
                y_sem = cellID_F_dff_lengthx_sample_location_sem(tempi,:);
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
                                    
            %for tempi=1:length(lengthx_sample_interval)
            for tempi=1:length(lengthx_sample_interval(1:end-1))
                plot(temp_ticks_sample(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
                    '-','LineWidth',1,'Color',[0 0 0]);
                hold on
            end
            
            for tempi=1:plot_lengthFlag
                x = [temp_ticks_sample(1+tempi-1) min_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1) max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 min_sample_decisionPeriodA];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            %% Plot dicisionPeriodA
            for tempi=1:numFrames
                %x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);                
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));                
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,1:length(x));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end

            
            %for tempi=1:length(decisionPeriodA_interval)
            %for tempi=1:(length(decisionPeriodA_interval)-1)
            for tempi=2:length(decisionPeriodA_interval(1:end-1))
                plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
                    '-','LineWidth',1,'Color',[0 0 0]);
                hold on
            end            
            
%             if plot_lengthFlag == 2
%                 %le = legend(h_line,'location1','location2','location3',...
%                 %    'location4','location5','location6','Location','northeast','fontsize',7);%10
%                 le = legend(h_line,'1','2','3',...
%                    '4','5','6','Location','southeast','fontsize',9,'NumColumns',2);%10--.8                
%                 le.ItemTokenSize = ones(1,6)*12;%10
%                 le.Color = backgounrdColor;                                
%             end
            %uistack(le,'top');
            a = 1;
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 temp_ticks_sample(end)+xlim_padSize2]);
            %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+7]);
            %xlim([1-xlim_padSize1 length3_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(end)+7]);                                    
            %xlim([1-xlim_padSize1 length3_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1]);                                                
            xlim([1 length3_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1/2]);                                                            
            ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            %ylim([min_sample_decisionPeriodA,temptemp_step*2]);
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])            
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            if plot_lengthFlag == 1
                %xticklabels({'PreFixation','Fixation','T1','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                xticklabels({'T1','Delay1',''});
            elseif plot_lengthFlag == 2
                xticklabels({'T1','T2','Delay1',''});
            elseif plot_lengthFlag == 3
                xticklabels({'T1','T2','T3','Delay1',''});
            end
            xtickangle(0);
            
            %set(gca,'ytick',[min_sample_decisionPeriodA:0.4:max_sample_decisionPeriodA])
            %set(gca,'ytick',[0:0.4:max_sample_decisionPeriodA]) %#ok<*NBRAK>

            temptemp_range = linspace(0,max_sample_decisionPeriodA,3);
            %temptemp_step = round((temptemp_range(2)-temptemp_range(1))*10)/10;
            %temptemp_step = floor((temptemp_range(2)-temptemp_range(1))*10)/10;
            %set(gca,'ytick',0:temptemp_step:temptemp_step*2)
            
            temptemp_step = 2*floor((temptemp_range(2)-temptemp_range(1))*1)/1;            
            %set(gca,'ytick',0:temptemp_step:temptemp_step*1)
            
            %yticks([0 1]);
            
            
            %text(temp_ticks_sample(end),min_sample_decisionPeriodA,'/','fontsize',12);
            %text(temp_ticks_decisionPeriodA(1),min_sample_decisionPeriodA,'/','fontsize',12);
            %text(temp_ticks_decisionPeriodA(end),min_sample_decisionPeriodA,'/','fontsize',12);
            text(temp_ticks_decisionPeriodA(2),min_sample_decisionPeriodA,'/','fontsize',12);            
            
            set(gca, 'FontSize', 9)
            set(gca,'box','off');% 取消右、上边框
            %ylabel(sprintf('dF/F'), 'FontSize', 12, 'FontWeight', 'bold');
            ylabel(sprintf('dF/F'), 'FontSize', 9);            
        end
    end
    
    
    
    
    if if_plot_singlePanel == 1
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 512*1.07 120+7]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 512*1.07*0.7-2 120+7]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        t = tiledlayout(1,1,'TileSpacing','tight','Padding','loose'); %#ok<*NASGU>        
        
        %t.Title.String = sprintf(' suite2pBID%d,%s \n Location tuning in length 1 2 3',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2);                
        %t.Title.FontSize = 11;
        %t.Title.Interpreter = 'none';        
        
        for plot_lengthFlag=1:1
        %if true
        %    plot_lengthFlag = 1;
            
            nexttile
            
            if plot_lengthFlag == 1
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length1_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length1_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length1_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length1_decisionPeriodA_location_sem;                
                lengthx_sample_interval = length1_sample_interval;
            elseif plot_lengthFlag == 2
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length2_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length2_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length2_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length2_decisionPeriodA_location_sem;                                
                lengthx_sample_interval = length2_sample_interval;
            elseif plot_lengthFlag == 3
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length3_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length3_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length3_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length3_decisionPeriodA_location_sem;                                
                lengthx_sample_interval = length3_sample_interval;
            end
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            %temp_ticks_decisionPeriodA = temp_ticks_sample(end)-1+decisionPeriodA_interval;
            
            
            %% Plot sample
            h_line = [];
            for tempi=1:numFrames
                x = temp_ticks_sample(1):temp_ticks_sample(end);
                y = cellID_F_dff_lengthx_sample_location_mean(tempi,:);
                %h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                h_line = [h_line plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>                
                hold on
                y_sem = cellID_F_dff_lengthx_sample_location_sem(tempi,:);
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
                                    
            %for tempi=1:length(lengthx_sample_interval)
            for tempi=1:length(lengthx_sample_interval(1:end-1))
                plot(temp_ticks_sample(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
                    '-','LineWidth',1,'Color',[0 0 0]);
                hold on
            end
            
            for tempi=1:plot_lengthFlag
                x = [temp_ticks_sample(1+tempi-1) min_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1) max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 min_sample_decisionPeriodA];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            %% Plot dicisionPeriodA
            for tempi=1:numFrames
                %x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);                
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));                
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,1:length(x));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end

            
            %for tempi=1:length(decisionPeriodA_interval)
            %for tempi=1:(length(decisionPeriodA_interval)-1)
            for tempi=2:length(decisionPeriodA_interval(1:end-1))
                plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
                    '-','LineWidth',1,'Color',[0 0 0]);
                hold on
            end            
            
%             if plot_lengthFlag == 2
%                 %le = legend(h_line,'location1','location2','location3',...
%                 %    'location4','location5','location6','Location','northeast','fontsize',7);%10
%                 le = legend(h_line,'1','2','3',...
%                    '4','5','6','Location','southeast','fontsize',9,'NumColumns',2);%10--.8                
%                 le.ItemTokenSize = ones(1,6)*12;%10
%                 le.Color = backgounrdColor;                                
%             end
            %uistack(le,'top');
            a = 1;
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 temp_ticks_sample(end)+xlim_padSize2]);
            %xlim([1-xlim_padSize1 temp_ticks_decisionPeriodA(end)+7]);
            %xlim([1-xlim_padSize1 length3_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(end)+7]);                                    
            %xlim([1-xlim_padSize1 length3_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1]);                                                
            xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1]);                                                            
            ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            %ylim([min_sample_decisionPeriodA,temptemp_step*2]);
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])            
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            if plot_lengthFlag == 1
                %xticklabels({'PreFixation','Fixation','T1','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                xticklabels({'T1','Delay1',''});
            elseif plot_lengthFlag == 2
                xticklabels({'T1','T2','Delay1',''});
            elseif plot_lengthFlag == 3
                xticklabels({'T1','T2','T3','Delay1',''});
            end
            xtickangle(0);
            
            %set(gca,'ytick',[min_sample_decisionPeriodA:0.4:max_sample_decisionPeriodA])
            %set(gca,'ytick',[0:0.4:max_sample_decisionPeriodA]) %#ok<*NBRAK>

            temptemp_range = linspace(0,max_sample_decisionPeriodA,3);
            %temptemp_step = round((temptemp_range(2)-temptemp_range(1))*10)/10;
            %temptemp_step = floor((temptemp_range(2)-temptemp_range(1))*10)/10;
            %set(gca,'ytick',0:temptemp_step:temptemp_step*2)
            
            temptemp_step = 2*floor((temptemp_range(2)-temptemp_range(1))*1)/1;            
            %set(gca,'ytick',0:temptemp_step:temptemp_step*1)
            
            %yticks([0 1]);
            
            
            %text(temp_ticks_sample(end),min_sample_decisionPeriodA,'/','fontsize',12);
            %text(temp_ticks_decisionPeriodA(1),min_sample_decisionPeriodA,'/','fontsize',12);
            %text(temp_ticks_decisionPeriodA(end),min_sample_decisionPeriodA,'/','fontsize',12);
            text(temp_ticks_decisionPeriodA(2),min_sample_decisionPeriodA,'/','fontsize',12);            
            
            set(gca, 'FontSize', 12)
            set(gca,'box','off');% 取消右、上边框
            %ylabel(sprintf('dF/F'), 'FontSize', 12, 'FontWeight', 'bold');
            ylabel(sprintf('dF/F'), 'FontSize', 12);            
        end
    end
    
    
    
    
    if if_plot_multiPanelB == 1
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 381*0.9 127*1.05*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977 127*1.05*0.95*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67 127*1.05*0.95*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        %t = tiledlayout(2,1,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>                    
        t = tiledlayout(2,10,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>       
        
        %nexttile
        nexttile([1 9])
        if true
            plot_lengthFlag = 1;
            
            cellID_F_dff_lengthx_sample_location = cellID_F_dff_length1_sample_location;
            cellID_F_dff_lengthx_decisionPeriodA_location = cellID_F_dff_length1_decisionPeriodA_location;
            
            cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length1_sample_location_mean;
            cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length1_sample_location_sem;
            cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length1_decisionPeriodA_location_mean;
            cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length1_decisionPeriodA_location_sem;
            lengthx_sample_interval = length1_sample_interval;
            
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            
            
            cellID_F_dff_lengthx_sample_location;
            cellID_F_dff_lengthx_sample_location_collapsed = [];
            cellID_F_dff_lengthx_decisionPeriodA_location;
            cellID_F_dff_lengthx_decisionPeriodA_location_collapsed = [];
            for tempi=1:numFrames
                cellID_F_dff_lengthx_sample_location_collapsed = [cellID_F_dff_lengthx_sample_location_collapsed; cellID_F_dff_lengthx_sample_location{tempi}];
                cellID_F_dff_lengthx_decisionPeriodA_location_collapsed = [cellID_F_dff_lengthx_decisionPeriodA_location_collapsed; cellID_F_dff_lengthx_decisionPeriodA_location{tempi}];
            end
            
            C1 = cellID_F_dff_lengthx_sample_location_collapsed;
            [C1_min,C1_max] = bounds(C1(:));
            
            C2 = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            [C2_min,C2_max] = bounds(C2(:));
            
            C_min = min([C1_min,C2_min]);
            C_min_valid = 0;
            
            C_max = max([C1_max,C2_max]);
            %C_max_valid = C_max*0.65;%0.5,0.65,0.55
            C_max_valid = max_sample_decisionPeriodA*2;
            if cellIndex_suite2p_temptempRaw==11 || cellIndex_suite2p_temptempRaw==19
                C_max_valid = 3.3;%5.5,2.5,3.3,2.3
            end
            
            temp_trialNum = nan(1,numFrames);
            for tempi=1:numFrames
                temp_trialNum(tempi) =  size(cellID_F_dff_lengthx_sample_location{tempi},1);
            end
            temp_lengthx_locationIndex = nan(1,numFrames+1);
            for tempi=1:(numFrames+1)
                temp_lengthx_locationIndex(tempi) = sum(temp_trialNum(1:tempi-1))+1;
            end
            
            
            x = [temp_ticks_sample(1) temp_ticks_sample(end)];
            y = [1,size(cellID_F_dff_lengthx_sample_location_collapsed,1)];
            C = cellID_F_dff_lengthx_sample_location_collapsed;
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(2)];
            y = [1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)];
            %C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed(:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %plot((temp_ticks_decisionPeriodA(1)-0.5)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
            %    '-','LineWidth',2,'Color',[1 1 1]);%3
            plot((temp_ticks_decisionPeriodA(1)-0.25)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
               '-','LineWidth',2,'Color',[1 1 1]);%3            
            hold on
            
            %             for tempi=1:length(lengthx_sample_interval(1:end-1))
            %                 plot(temp_ticks_sample(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
            %                     '-','LineWidth',1,'Color',[0 0 0]);
            %                 hold on
            %             end
            %
            %             for tempi=2:length(decisionPeriodA_interval(1:end-1))
            %                 plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
            %                     '-','LineWidth',1,'Color',[0 0 0]);
            %                 hold on
            %             end
            
            
            
            
            
            %             for tempi=1:numFrames
            %                 %plot(temp_ticks_decisionPeriodA(2)+[0.5 5],(temp_lengthx_locationIndex(tempi)-0.5)*[1 1],...
            %                 %    '-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
            %                 %hold on
            %
            %                     x = [temp_ticks_decisionPeriodA(2)+0.5 temp_lengthx_locationIndex(tempi)-0.5;...
            %                         temp_ticks_decisionPeriodA(2)+0.5 temp_lengthx_locationIndex(tempi+1)-1.5;...
            %                         temp_ticks_decisionPeriodA(2)+10 temp_lengthx_locationIndex(tempi+1)-1.5;...
            %                         temp_ticks_decisionPeriodA(2)+10 temp_lengthx_locationIndex(tempi)-0.5];
            %                     y = [1 2 3 4];
            %                     patch('Faces',y,'Vertices',x,'FaceColor',multi_rgbColor(tempi,:),'FaceAlpha',1,'EdgeColor','none');%0.1
            %                     hold on
            %
            %             end
            
            for tempi=1:numFrames                
            %if true
                
                temp_barLength = 3.5;%4
%                 
%                 tempi = 5;
%                 if cellIndex_suite2p_temptempRaw == 11
%                     tempi = 3;
%                 end                
                x = [temp_ticks_decisionPeriodA(2)+0.5 temp_lengthx_locationIndex(tempi)-0.5;...
                    temp_ticks_decisionPeriodA(2)+0.5 temp_lengthx_locationIndex(tempi+1)-1.5;...
                    temp_ticks_decisionPeriodA(2)+temp_barLength temp_lengthx_locationIndex(tempi+1)-1.5;...
                    temp_ticks_decisionPeriodA(2)+temp_barLength temp_lengthx_locationIndex(tempi)-0.5];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',multi_rgbColor(tempi,:),'FaceAlpha',1,'EdgeColor','none');%0.1
                hold on     
                
                %text(temp_ticks_decisionPeriodA(2)+temp_barLength+0.25,...
                %    -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
                %    sprintf('Loc%d',tempi),'fontsize',6.5,'color',multi_rgbColor(tempi,:));
                
                temp_colorStr = ['Loc','\color[rgb]',sprintf('{%f %f %f}',multi_rgbColor(tempi,:)),sprintf('%d',tempi)];
                
                text(temp_ticks_decisionPeriodA(2)+temp_barLength+0.25,...
                   -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
                   temp_colorStr,'fontsize',6.5);

            end
            
%             for tempi=1:numFrames
%                 %plot([0.5 temp_ticks_decisionPeriodA(2)+0.5],(temp_lengthx_locationIndex(tempi)-0.5)*[1 1],...
%                 %    ':','LineWidth',0.5,'Color',[1 1 1]);
%                 h = plot([0.5 temp_ticks_decisionPeriodA(2)+0.5],(temp_lengthx_locationIndex(tempi)-0.5)*[1 1],...
%                    '--','LineWidth',0.5,'Color',[1 1 1 0.3]);                
%                 hold on
%             end
            
            
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);
            xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);            
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            ylim([1-1 sum(temp_trialNum)+1]);
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticklabels('');
            
            %yticks([0 400]);
            yticks([]);
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('Trials'), 'FontSize', 8);
            %colorbar;
            
            %c = colorbar('westoutside','FontSize',7);
            c = colorbar('FontSize',6.5);%9
            %c.Position = c.Position+[0 0 -0.0375 -0.23];
            %c.Position = c.Position+[-0.01 0.095 -0.02 -0.1];
            c.Position = c.Position+[-0.045 0.095 -0.038 -0.1];
            
            c.Ticks = [0 floor(C_max_valid)];
            
            if if_colormap_loadEnhanced == 1
                load('parula_enhanced');
                colormap(parula_enhanced); 
                %load('gray_enhanced');
                %colormap(gray_enhanced);
            elseif if_colormap_loadEnhanced == 0
                colormap parula
                %colormap gray
            end            
            
            
        end
        
        
        
        %nexttile
        nexttile([1 9])
        
        temp_LineWidth = 0.5;%1
        
        if true
            plot_lengthFlag = 1;
            
            if plot_lengthFlag == 1
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length1_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length1_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length1_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length1_decisionPeriodA_location_sem;
                lengthx_sample_interval = length1_sample_interval;
            elseif plot_lengthFlag == 2
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length2_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length2_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length2_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length2_decisionPeriodA_location_sem;
                lengthx_sample_interval = length2_sample_interval;
            elseif plot_lengthFlag == 3
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length3_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length3_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length3_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length3_decisionPeriodA_location_sem;
                lengthx_sample_interval = length3_sample_interval;
            end
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            %temp_ticks_decisionPeriodA = temp_ticks_sample(end)-1+decisionPeriodA_interval;
            
            
            %% Plot sample
            h_line = [];
            for tempi=1:numFrames
                x = temp_ticks_sample(1):temp_ticks_sample(end);
                y = cellID_F_dff_lengthx_sample_location_mean(tempi,:);
                %h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                h_line = [h_line plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                hold on
                y_sem = cellID_F_dff_lengthx_sample_location_sem(tempi,:);
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            
            %             %for tempi=1:length(lengthx_sample_interval)
            %             for tempi=1:length(lengthx_sample_interval(1:end-1))
            %                 plot(temp_ticks_sample(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
            %                     '-','LineWidth',1,'Color',[0 0 0]);
            %                 hold on
            %             end
            
            for tempi=1:plot_lengthFlag
                x = [temp_ticks_sample(1+tempi-1) min_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1) max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 min_sample_decisionPeriodA];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            %% Plot dicisionPeriodA
            for tempi=1:numFrames
                %x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,1:length(x));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            
            
            %for tempi=1:length(decisionPeriodA_interval)
            %for tempi=1:(length(decisionPeriodA_interval)-1)
            %for tempi=2:length(decisionPeriodA_interval(1:end-1))
            %    plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
            %        '-','LineWidth',1,'Color',[0 0 0]);
            %    hold on
            %end
            
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                            
            xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                                        
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);            
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])            
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            if plot_lengthFlag == 1
                xticklabels({'T1','Delay1',''});
            elseif plot_lengthFlag == 2
                xticklabels({'T1','T2','Delay1',''});
            elseif plot_lengthFlag == 3
                xticklabels({'T1','T2','T3','Delay1',''});
            end
            xtickangle(0);
            
            temptemp_range = linspace(0,max_sample_decisionPeriodA,3);
            
            temptemp_step = 2*floor((temptemp_range(2)-temptemp_range(1))*1)/1;            
            text(temp_ticks_decisionPeriodA(2),min_sample_decisionPeriodA,'/','fontsize',12);            
            
            if cellIndex_suite2p_temptempRaw == 19
                yticks([0 2]);
            end
            if cellIndex_suite2p_temptempRaw == 11
                yticks([0 1]);
            end
            if cellIndex_suite2p_temptempRaw == 121
                yticks([0 3]);
            end
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('dF/F'), 'FontSize', 8);            
        end
        
        
    end    
    
    if if_plot_multiPanelB2 == 1
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 381*0.9 127*1.05*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977 127*1.05*0.95*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04 127*1.05*0.95*0.93*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        %t = tiledlayout(2,1,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>                    
        t = tiledlayout(2,10,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>       
        
        %nexttile
        nexttile([1 9])
        if true
            plot_lengthFlag = 1;
            
            cellID_F_dff_lengthx_sample_location = cellID_F_dff_length1_sample_location;
            cellID_F_dff_lengthx_decisionPeriodA_location = cellID_F_dff_length1_decisionPeriodA_location;
            
            cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length1_sample_location_mean;
            cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length1_sample_location_sem;
            cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length1_decisionPeriodA_location_mean;
            cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length1_decisionPeriodA_location_sem;
            lengthx_sample_interval = length1_sample_interval;
            
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            
            
            cellID_F_dff_lengthx_sample_location;
            cellID_F_dff_lengthx_sample_location_collapsed = [];
            cellID_F_dff_lengthx_decisionPeriodA_location;
            cellID_F_dff_lengthx_decisionPeriodA_location_collapsed = [];
            for tempi=1:numFrames
                cellID_F_dff_lengthx_sample_location_collapsed = [cellID_F_dff_lengthx_sample_location_collapsed; cellID_F_dff_lengthx_sample_location{tempi}];
                cellID_F_dff_lengthx_decisionPeriodA_location_collapsed = [cellID_F_dff_lengthx_decisionPeriodA_location_collapsed; cellID_F_dff_lengthx_decisionPeriodA_location{tempi}];
            end
            
            C1 = cellID_F_dff_lengthx_sample_location_collapsed;
            [C1_min,C1_max] = bounds(C1(:));
            
            C2 = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            [C2_min,C2_max] = bounds(C2(:));
            
            C_min = min([C1_min,C2_min]);
            C_min_valid = 0;
            
            C_max = max([C1_max,C2_max]);
            %C_max_valid = C_max*0.65;%0.5,0.65,0.55
            C_max_valid = max_sample_decisionPeriodA*2;
            if cellIndex_suite2p_temptempRaw==11 || cellIndex_suite2p_temptempRaw==19
                C_max_valid = 3.3;%5.5,2.5,3.3,2.3
            end
            
            temp_trialNum = nan(1,numFrames);
            for tempi=1:numFrames
                temp_trialNum(tempi) =  size(cellID_F_dff_lengthx_sample_location{tempi},1);
            end
            temp_lengthx_locationIndex = nan(1,numFrames+1);
            for tempi=1:(numFrames+1)
                temp_lengthx_locationIndex(tempi) = sum(temp_trialNum(1:tempi-1))+1;
            end
            
            
            x = [temp_ticks_sample(1) temp_ticks_sample(end)];
            y = [1,size(cellID_F_dff_lengthx_sample_location_collapsed,1)];
            C = cellID_F_dff_lengthx_sample_location_collapsed;
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(2)];
            x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];            
            y = [1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)];
            %C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed(:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %plot((temp_ticks_decisionPeriodA(1)-0.5)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
            %    '-','LineWidth',2,'Color',[1 1 1]);%3
            plot((temp_ticks_decisionPeriodA(1)-0.25)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
               '-','LineWidth',2,'Color',[1 1 1]);%3            
            hold on
            
            %             for tempi=1:length(lengthx_sample_interval(1:end-1))
            %                 plot(temp_ticks_sample(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
            %                     '-','LineWidth',1,'Color',[0 0 0]);
            %                 hold on
            %             end
            %
            %             for tempi=2:length(decisionPeriodA_interval(1:end-1))
            %                 plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
            %                     '-','LineWidth',1,'Color',[0 0 0]);
            %                 hold on
            %             end
            
            
            
            
            
            %             for tempi=1:numFrames
            %                 %plot(temp_ticks_decisionPeriodA(2)+[0.5 5],(temp_lengthx_locationIndex(tempi)-0.5)*[1 1],...
            %                 %    '-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
            %                 %hold on
            %
            %                     x = [temp_ticks_decisionPeriodA(2)+0.5 temp_lengthx_locationIndex(tempi)-0.5;...
            %                         temp_ticks_decisionPeriodA(2)+0.5 temp_lengthx_locationIndex(tempi+1)-1.5;...
            %                         temp_ticks_decisionPeriodA(2)+10 temp_lengthx_locationIndex(tempi+1)-1.5;...
            %                         temp_ticks_decisionPeriodA(2)+10 temp_lengthx_locationIndex(tempi)-0.5];
            %                     y = [1 2 3 4];
            %                     patch('Faces',y,'Vertices',x,'FaceColor',multi_rgbColor(tempi,:),'FaceAlpha',1,'EdgeColor','none');%0.1
            %                     hold on
            %
            %             end
            
            for tempi=1:numFrames                
            %if true
                
                temp_barLength = 3.5;%4
%                 
%                 tempi = 5;
%                 if cellIndex_suite2p_temptempRaw == 11
%                     tempi = 3;
%                 end                
                x = [temp_ticks_decisionPeriodA(3)+0.5 temp_lengthx_locationIndex(tempi)-0.5;...
                    temp_ticks_decisionPeriodA(3)+0.5 temp_lengthx_locationIndex(tempi+1)-1.5;...
                    temp_ticks_decisionPeriodA(3)+temp_barLength temp_lengthx_locationIndex(tempi+1)-1.5;...
                    temp_ticks_decisionPeriodA(3)+temp_barLength temp_lengthx_locationIndex(tempi)-0.5];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',multi_rgbColor(tempi,:),'FaceAlpha',1,'EdgeColor','none');%0.1
                hold on     
                
                %text(temp_ticks_decisionPeriodA(2)+temp_barLength+0.25,...
                %    -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
                %    sprintf('Loc%d',tempi),'fontsize',6.5,'color',multi_rgbColor(tempi,:));
                
                temp_colorStr = ['Loc','\color[rgb]',sprintf('{%f %f %f}',multi_rgbColor(tempi,:)),sprintf('%d',tempi)];
                
                text(temp_ticks_decisionPeriodA(3)+temp_barLength+0.25,...
                   -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
                   temp_colorStr,'fontsize',6.5);

            end
            
%             for tempi=1:numFrames
%                 %plot([0.5 temp_ticks_decisionPeriodA(2)+0.5],(temp_lengthx_locationIndex(tempi)-0.5)*[1 1],...
%                 %    ':','LineWidth',0.5,'Color',[1 1 1]);
%                 h = plot([0.5 temp_ticks_decisionPeriodA(2)+0.5],(temp_lengthx_locationIndex(tempi)-0.5)*[1 1],...
%                    '--','LineWidth',0.5,'Color',[1 1 1 0.3]);                
%                 hold on
%             end
            
            
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);            
            xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);            
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            ylim([1-1 sum(temp_trialNum)+1]);
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticklabels('');
            
            %yticks([0 400]);
            yticks([]);
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('Trials'), 'FontSize', 8);
            %colorbar;
            
            %c = colorbar('westoutside','FontSize',7);
            c = colorbar('FontSize',6.5);%9
            %c.Position = c.Position+[0 0 -0.0375 -0.23];
            %c.Position = c.Position+[-0.01 0.095 -0.02 -0.1];
            %c.Position = c.Position+[-0.045 0.095 -0.038 -0.1];
            c.Position = c.Position+[-0.045 0.080 -0.038 -0.1];
            
            c.Ticks = [0 floor(C_max_valid)];
            
            if if_colormap_loadEnhanced == 1
                load('parula_enhanced');
                colormap(parula_enhanced); 
                %load('gray_enhanced');
                %colormap(gray_enhanced);
            elseif if_colormap_loadEnhanced == 0
                colormap parula
                %colormap gray
            end            
            
            
        end
        
        
        
        %nexttile
        nexttile([1 9])
        
        temp_LineWidth = 0.5;%1
        
        if true
            plot_lengthFlag = 1;
            
            if plot_lengthFlag == 1
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length1_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length1_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length1_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length1_decisionPeriodA_location_sem;
                lengthx_sample_interval = length1_sample_interval;
            elseif plot_lengthFlag == 2
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length2_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length2_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length2_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length2_decisionPeriodA_location_sem;
                lengthx_sample_interval = length2_sample_interval;
            elseif plot_lengthFlag == 3
                cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_length3_sample_location_mean;
                cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_length3_sample_location_sem;
                cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_length3_decisionPeriodA_location_mean;
                cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_length3_decisionPeriodA_location_sem;
                lengthx_sample_interval = length3_sample_interval;
            end
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            %temp_ticks_decisionPeriodA = temp_ticks_sample(end)-1+decisionPeriodA_interval;
            
            
            %% Plot sample
            h_line = [];
            for tempi=1:numFrames
                x = temp_ticks_sample(1):temp_ticks_sample(end);
                y = cellID_F_dff_lengthx_sample_location_mean(tempi,:);
                %h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                h_line = [h_line plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                hold on
                y_sem = cellID_F_dff_lengthx_sample_location_sem(tempi,:);
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            
            %             %for tempi=1:length(lengthx_sample_interval)
            %             for tempi=1:length(lengthx_sample_interval(1:end-1))
            %                 plot(temp_ticks_sample(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
            %                     '-','LineWidth',1,'Color',[0 0 0]);
            %                 hold on
            %             end
            
            for tempi=1:plot_lengthFlag
                x = [temp_ticks_sample(1+tempi-1) min_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1) max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 min_sample_decisionPeriodA];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            %% Plot dicisionPeriodA
            for tempi=1:numFrames
                x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                %x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,1:length(x));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            
            
            %for tempi=1:length(decisionPeriodA_interval)
            %for tempi=1:(length(decisionPeriodA_interval)-1)
            %for tempi=2:length(decisionPeriodA_interval(1:end-1))
            %    plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_sample_decisionPeriodA max_sample_decisionPeriodA],...
            %        '-','LineWidth',1,'Color',[0 0 0]);
            %    hold on
            %end
            
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                            
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                                        
            xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);            
            ylim([min_sample_decisionPeriodA-(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.05,...
                max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);                        
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])            
            %xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end)]);            
            if plot_lengthFlag == 1
                xticklabels({'T1','Delay1','ChoiceCue'});
            elseif plot_lengthFlag == 2
                xticklabels({'T1','T2','Delay1','ChoiceCue'});
            elseif plot_lengthFlag == 3
                xticklabels({'T1','T2','T3','Delay1','ChoiceCue'});
            end
            xtickangle(0);
            
            temptemp_range = linspace(0,max_sample_decisionPeriodA,3);
            
            temptemp_step = 2*floor((temptemp_range(2)-temptemp_range(1))*1)/1;            
            %text(temp_ticks_decisionPeriodA(2),min_sample_decisionPeriodA,'/','fontsize',12);            
            text(temp_ticks_decisionPeriodA(3),min_sample_decisionPeriodA,'/','fontsize',12);            
            
            if cellIndex_suite2p_temptempRaw == 19
                yticks([0 2]);
            end
            if cellIndex_suite2p_temptempRaw == 11
                yticks([0 1]);
            end
            if cellIndex_suite2p_temptempRaw == 121
                yticks([0 3]);
            end
            if cellIndex_suite2p_temptempRaw == 25
                yticks([0 1]);
            end
                
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('dF/F'), 'FontSize', 8);            
        end
        
        
    end   
    
    if if_plot_multiPanelB3 == 1
                        
        if cellIndex_suite2p_temptempRaw==19
            max_sample_decisionPeriodA = 2;
        elseif cellIndex_suite2p_temptempRaw==25
            max_sample_decisionPeriodA = 0.8;
        end
        
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04 127*1.05*0.95*0.93*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 233.4*0.93 176.7*0.86]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        t = tiledlayout(2,10,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>       
        
        nexttile([1 9])
        if true
            plot_lengthFlag = 1;
            
            cellID_F_dff_lengthx_sample_location = cellID_F_dff_lastT_sample_location;
            cellID_F_dff_lengthx_decisionPeriodA_location = cellID_F_dff_lastT_decisionPeriodA_location;
            
            cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_lastT_sample_location_mean;
            cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_lastT_sample_location_sem;
            cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_lastT_decisionPeriodA_location_mean;
            cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_lastT_decisionPeriodA_location_sem;
            lengthx_sample_interval = length1_sample_interval;
            
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            
            
            cellID_F_dff_lengthx_sample_location;
            cellID_F_dff_lengthx_sample_location_collapsed = [];
            cellID_F_dff_lengthx_decisionPeriodA_location;
            cellID_F_dff_lengthx_decisionPeriodA_location_collapsed = [];
            for tempi=1:numFrames
                cellID_F_dff_lengthx_sample_location_collapsed = [cellID_F_dff_lengthx_sample_location_collapsed; cellID_F_dff_lengthx_sample_location{tempi}];
                cellID_F_dff_lengthx_decisionPeriodA_location_collapsed = [cellID_F_dff_lengthx_decisionPeriodA_location_collapsed; cellID_F_dff_lengthx_decisionPeriodA_location{tempi}];
            end
            
            C1 = cellID_F_dff_lengthx_sample_location_collapsed;
            [C1_min,C1_max] = bounds(C1(:));
            
            C2 = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            [C2_min,C2_max] = bounds(C2(:));
            
            C_min = min([C1_min,C2_min]);
            C_min_valid = 0;
            
            C_max = max([C1_max,C2_max]);
            %C_max_valid = C_max*0.65;%0.5,0.65,0.55
            C_max_valid = max_sample_decisionPeriodA*2;
            if cellIndex_suite2p_temptempRaw==11 || cellIndex_suite2p_temptempRaw==19
                %C_max_valid = 3.3;%5.5,2.5,3.3,2.3
                C_max_valid = 2.1;%2.5
            end
            
            temp_trialNum = nan(1,numFrames);
            for tempi=1:numFrames
                temp_trialNum(tempi) =  size(cellID_F_dff_lengthx_sample_location{tempi},1);
            end
            temp_lengthx_locationIndex = nan(1,numFrames+1);
            for tempi=1:(numFrames+1)
                temp_lengthx_locationIndex(tempi) = sum(temp_trialNum(1:tempi-1))+1;
            end
            
            
            x = [temp_ticks_sample(1) temp_ticks_sample(end)];
            y = [1,size(cellID_F_dff_lengthx_sample_location_collapsed,1)];
            C = cellID_F_dff_lengthx_sample_location_collapsed;
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(2)];
            x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];            
            %x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]-0.5;
            y = [1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)];
            %C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed(:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %plot((temp_ticks_decisionPeriodA(1)-0.5)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
            %    '-','LineWidth',2,'Color',[1 1 1]);%3
            plot((temp_ticks_decisionPeriodA(1)-0.25)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
              '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            hold on
            
            plot((temp_ticks_decisionPeriodA(2)+0.25)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
              '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            hold on
            

            
            for tempi=1:numFrames                
            %if true
                
                temp_barLength = 3.5;%4
                
                x = [temp_ticks_decisionPeriodA(3)+0.5 temp_lengthx_locationIndex(tempi)-0.5;...
                    temp_ticks_decisionPeriodA(3)+0.5 temp_lengthx_locationIndex(tempi+1)-1.5;...
                    temp_ticks_decisionPeriodA(3)+temp_barLength temp_lengthx_locationIndex(tempi+1)-1.5;...
                    temp_ticks_decisionPeriodA(3)+temp_barLength temp_lengthx_locationIndex(tempi)-0.5];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',multi_rgbColor(tempi,:),'FaceAlpha',1,'EdgeColor','none');%0.1
                hold on     
                
                %text(temp_ticks_decisionPeriodA(2)+temp_barLength+0.25,...
                %    -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
                %    sprintf('Loc%d',tempi),'fontsize',6.5,'color',multi_rgbColor(tempi,:));
                
                temp_colorStr = ['Loc','\color[rgb]',sprintf('{%f %f %f}',multi_rgbColor(tempi,:)),sprintf('%d',tempi)];
                
                text(temp_ticks_decisionPeriodA(3)+temp_barLength+0.25,...
                   -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
                   temp_colorStr,'fontsize',6.5);

            end
                        
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);            
            xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);            
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            ylim([1-1 sum(temp_trialNum)+1]);
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticklabels('');
            
            %yticks([0 400]);
            yticks([]);
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('Trials'), 'FontSize', 8);
            %colorbar;
            
            %c = colorbar('westoutside','FontSize',7);
            c = colorbar('FontSize',6.5);%9
            c.Position = c.Position+[-0.045 0.080 -0.038 -0.1];
            
            c.Ticks = [0 floor(C_max_valid)];
            
            if if_colormap_loadEnhanced == 1
                %load('parula_enhanced');
                %colormap(parula_enhanced);
                %load('gray_enhanced');
                %colormap(gray_enhanced);
                
                %load('gray_reversed_enhanced');
                %colormap(gray_reversed_enhanced);
                
                load('coolwarm_enhanced');
                colormap(coolwarm_enhanced);
                
                %load('coolwarm_n11n_enhanced');
                %colormap(coolwarm_n11n_enhanced);

            elseif if_colormap_loadEnhanced == 0
                %colormap parula
                %colormap gray
                
                %temp1 = gray;
                %temp1 = temp1(end:-1:1,:);
                %colormap(temp1);
                
                colormap(coolwarm());
                
                %temp1 = coolwarm(300);
                %temp2 = ((300-256)/2)+1;
                %temp3 = temp1(temp2:temp2+255,:);
                %colormap(temp3);
            end
            
            
        end
        
        
        
        nexttile([1 9])
        
        temp_LineWidth = 0.5;%1
        
        if true
            plot_lengthFlag = 1;
            
            cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_lastT_sample_location_mean;
            cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_lastT_sample_location_sem;
            cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_lastT_decisionPeriodA_location_mean;
            cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_lastT_decisionPeriodA_location_sem;
            lengthx_sample_interval = length1_sample_interval;
            
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            
            
            %% Plot sample
            h_line = [];
            for tempi=1:numFrames
                x = temp_ticks_sample(1):temp_ticks_sample(end);
                y = cellID_F_dff_lengthx_sample_location_mean(tempi,:);
                %h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                h_line = [h_line plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                hold on
                y_sem = cellID_F_dff_lengthx_sample_location_sem(tempi,:);
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            
            for tempi=1:plot_lengthFlag
                x = [temp_ticks_sample(1+tempi-1) min_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1) max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 min_sample_decisionPeriodA];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            %% Plot dicisionPeriodA
            for tempi=1:numFrames
                %x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);
                %y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,decisionPeriodA_interval(1):decisionPeriodA_interval(2));                
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            for tempi=1:numFrames
                x = temp_ticks_decisionPeriodA(2)+1:temp_ticks_decisionPeriodA(end);
                %y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end            
                                   
            %plot((temp_ticks_decisionPeriodA(1))*[1 1],[min_sample_decisionPeriodA,max_sample_decisionPeriodA],...
            %   '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            %hold on
            %plot((temp_ticks_decisionPeriodA(2)-0.25)*[1 1],[min_sample_decisionPeriodA,max_sample_decisionPeriodA],...
            %   '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            %hold on            
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                            
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                                        
            xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);            
            ylim([min_sample_decisionPeriodA-(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.05,...
                max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);                        
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])            
            %xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end)]);            
            %xticklabels({'T1','Delay1','ChoiceCue'});
            %xticklabels({'LastT','Delay1','Decision'});
            xticklabels({'LastT','Delay-on','Decision'});            
            xtickangle(0);
            
            temptemp_range = linspace(0,max_sample_decisionPeriodA,3);
            
            temptemp_step = 2*floor((temptemp_range(2)-temptemp_range(1))*1)/1;            
            %text(temp_ticks_decisionPeriodA(2),min_sample_decisionPeriodA,'/','fontsize',12);            
            text(temp_ticks_decisionPeriodA(3),min_sample_decisionPeriodA,'/','fontsize',12);            
            
            if cellIndex_suite2p_temptempRaw == 19
                yticks([0 2]);
            end
            if cellIndex_suite2p_temptempRaw == 11
                yticks([0 1]);
            end
            if cellIndex_suite2p_temptempRaw == 121
                yticks([0 3]);
            end
            if cellIndex_suite2p_temptempRaw == 25
                yticks([0 0.7]);
            end
                
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('dF/F'), 'FontSize', 8);            
        end
        
        
    end        
    
    
    
    if if_plot_multiPanelB3_2 == 1
                        
        if cellIndex_suite2p_temptempRaw==19
            max_sample_decisionPeriodA = 2;
        elseif cellIndex_suite2p_temptempRaw==25
            max_sample_decisionPeriodA = 0.8;
            min_sample_decisionPeriodA = 0;
        end
        
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04 127*1.05*0.95*0.93*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 233.4*0.93 176.7*0.86]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        t = tiledlayout(2,10,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>       
        
        nexttile([1 9])
        if true
            plot_lengthFlag = 1;
            
            
            cellID_F_dff_lastT_sample_seqSubset;
            cellID_F_dff_lastT_sample_nonSeqSubset;
            
            cellID_F_dff_decisionPeriodA_seqSubset;
            cellID_F_dff_decisionPeriodA_nonSeqSubset;
                        
            cellID_F_dff_lastT_sample_x_collapsed = [cellID_F_dff_lastT_sample_seqSubset;cellID_F_dff_lastT_sample_nonSeqSubset];
            cellID_F_dff_decisionPeriodA_x_collapsed = [cellID_F_dff_decisionPeriodA_seqSubset;cellID_F_dff_decisionPeriodA_nonSeqSubset];
                                    
            tempTrialBoundary = size(cellID_F_dff_lastT_sample_seqSubset,1);
            
            temp_trialNum = [size(cellID_F_dff_lastT_sample_seqSubset,1),size(cellID_F_dff_lastT_sample_nonSeqSubset,1)];
            temp_trialIndex = nan(1,length(temp_trialNum)+1);
            for tempi=1:length(temp_trialIndex)
                temp_trialIndex(tempi) = sum(temp_trialNum(1:tempi-1))+1;
            end
            
            lengthx_sample_interval = length1_sample_interval;            
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            
            
            
            
            C1 = cellID_F_dff_lastT_sample_x_collapsed;
            [C1_min,C1_max] = bounds(C1(:));
            
            C2 = cellID_F_dff_decisionPeriodA_x_collapsed;
            [C2_min,C2_max] = bounds(C2(:));
            
            C_min = min([C1_min,C2_min]);
            C_min_valid = 0;
            
            C_max = max([C1_max,C2_max]);
            %C_max_valid = C_max*0.65;%0.5,0.65,0.55
            C_max_valid = max_sample_decisionPeriodA*2;
            %C_max_valid = 1.4;
            if cellIndex_suite2p_temptempRaw==11 || cellIndex_suite2p_temptempRaw==19
                %C_max_valid = 3.3;%5.5,2.5,3.3,2.3
                C_max_valid = 2.5;
            end
            
            if cellIndex_suite2p_temptempRaw == 25
                C_max_valid = 1.4;
                C_min_valid = 0.4;
            end
            
            %             temp_trialNum = nan(1,numFrames);
            %             for tempi=1:numFrames
            %                 temp_trialNum(tempi) =  size(cellID_F_dff_lengthx_sample_location{tempi},1);
            %             end
            %             temp_lengthx_locationIndex = nan(1,numFrames+1);
            %             for tempi=1:(numFrames+1)
            %                 temp_lengthx_locationIndex(tempi) = sum(temp_trialNum(1:tempi-1))+1;
            %             end
            
            
            x = [temp_ticks_sample(1) temp_ticks_sample(end)];
            y = [1,size(cellID_F_dff_lastT_sample_x_collapsed,1)];
            C = cellID_F_dff_lastT_sample_x_collapsed;
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(2)];
            x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];            
            %x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]-0.5;
            y = [1,size(cellID_F_dff_decisionPeriodA_x_collapsed,1)];
            %C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            C = cellID_F_dff_decisionPeriodA_x_collapsed(:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %plot((temp_ticks_decisionPeriodA(1)-0.5)*[1 1],[1,size(cellID_F_dff_decisionPeriodA_x_collapsed,1)],...
            %    '-','LineWidth',2,'Color',[1 1 1]);%3
            plot((temp_ticks_decisionPeriodA(1)-0.25)*[1 1],[1,size(cellID_F_dff_decisionPeriodA_x_collapsed,1)],...
              '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            hold on
            
            plot((temp_ticks_decisionPeriodA(2)+0.25)*[1 1],[1,size(cellID_F_dff_decisionPeriodA_x_collapsed,1)],...
              '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            hold on
            
            
            plot([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2],...
                [1 1].*(temp_trialIndex(2)+1),...
              '--','LineWidth',1,'Color',[0 0 0]);%3,2           
            hold on
            
            for tempi=1:(length(temp_trialIndex)-1)
                if tempi == 1
                    %temp_str = sprintf('Subset\nof 126');
                    temp_str = [sprintf('Subset\nof '),...
                        '\color[rgb]',sprintf('{%f %f %f}',multi_rgbColor(1,:)),sprintf('1'),...
                        '\color[rgb]',sprintf('{%f %f %f}',multi_rgbColor(2,:)),sprintf('2'),...
                        '\color[rgb]',sprintf('{%f %f %f}',multi_rgbColor(6,:)),sprintf('6')];
                elseif tempi == 2
                    temp_str = sprintf('Others');
                end
                text(temp_ticks_decisionPeriodA(3)+0.7,...
                    -20+(temp_trialIndex(tempi)+temp_trialIndex(tempi+1))/2,...
                    temp_str,'fontsize',6);%6.5
            end

%             for tempi=1:numFrames                
%             %if true
%                 
%                 temp_barLength = 3.5;%4
%                 
%                 x = [temp_ticks_decisionPeriodA(3)+0.5 temp_lengthx_locationIndex(tempi)-0.5;...
%                     temp_ticks_decisionPeriodA(3)+0.5 temp_lengthx_locationIndex(tempi+1)-1.5;...
%                     temp_ticks_decisionPeriodA(3)+temp_barLength temp_lengthx_locationIndex(tempi+1)-1.5;...
%                     temp_ticks_decisionPeriodA(3)+temp_barLength temp_lengthx_locationIndex(tempi)-0.5];
%                 y = [1 2 3 4];
%                 patch('Faces',y,'Vertices',x,'FaceColor',multi_rgbColor(tempi,:),'FaceAlpha',1,'EdgeColor','none');%0.1
%                 hold on     
%                 
%                 %text(temp_ticks_decisionPeriodA(2)+temp_barLength+0.25,...
%                 %    -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
%                 %    sprintf('Loc%d',tempi),'fontsize',6.5,'color',multi_rgbColor(tempi,:));
%                 
%                 temp_colorStr = ['Loc','\color[rgb]',sprintf('{%f %f %f}',multi_rgbColor(tempi,:)),sprintf('%d',tempi)];
%                 
%                 text(temp_ticks_decisionPeriodA(3)+temp_barLength+0.25,...
%                    -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
%                    temp_colorStr,'fontsize',6.5);
% 
%             end
                        
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);            
            xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);            
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            ylim([1-1 sum(temp_trialNum)+1]);
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticklabels('');
            
            %yticks([0 400]);
            yticks([]);
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('Trials'), 'FontSize', 8);
            %colorbar;
            
            %c = colorbar('westoutside','FontSize',7);
            c = colorbar('FontSize',6.5);%9
            c.Position = c.Position+[-0.045 0.080 -0.038 -0.1];
            
            %c.Ticks = [0 floor(C_max_valid)];
            c.Ticks = [C_min_valid floor(C_max_valid)];
           
            if if_colormap_loadEnhanced == 1
                %load('parula_enhanced');
                %colormap(parula_enhanced);
                %load('gray_enhanced');
                %colormap(gray_enhanced);
                
                %load('gray_reversed_enhanced');
                %colormap(gray_reversed_enhanced);
                
                load('coolwarm_enhanced');
                colormap(coolwarm_enhanced);
                
                %load('coolwarm_n11n_enhanced');
                %colormap(coolwarm_n11n_enhanced);

            elseif if_colormap_loadEnhanced == 0
                %colormap parula
                %colormap gray
                
                %temp1 = gray;
                %temp1 = temp1(end:-1:1,:);
                %colormap(temp1);
                
                colormap(coolwarm());
                
                %temp1 = coolwarm(300);
                %temp2 = ((300-256)/2)+1;
                %temp3 = temp1(temp2:temp2+255,:);
                %colormap(temp3);
            end
        
            
            
        end
        
        
        
        nexttile([1 9])
        
        temp_LineWidth = 0.5;%1
        
        if true
            plot_lengthFlag = 1;
            
            cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_lastT_sample_location_mean;
            cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_lastT_sample_location_sem;
            cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_lastT_decisionPeriodA_location_mean;
            cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_lastT_decisionPeriodA_location_sem;
            lengthx_sample_interval = length1_sample_interval;
            
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            
            
            %% Plot sample
            h_line = [];
            for tempi=1:numFrames
                x = temp_ticks_sample(1):temp_ticks_sample(end);
                y = cellID_F_dff_lengthx_sample_location_mean(tempi,:);
                %h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                h_line = [h_line plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                hold on
                y_sem = cellID_F_dff_lengthx_sample_location_sem(tempi,:);
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            
            for tempi=1:plot_lengthFlag
                x = [temp_ticks_sample(1+tempi-1) min_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1) max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 min_sample_decisionPeriodA];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            %% Plot dicisionPeriodA
            for tempi=1:numFrames
                %x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);
                %y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,decisionPeriodA_interval(1):decisionPeriodA_interval(2));                
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            for tempi=1:numFrames
                x = temp_ticks_decisionPeriodA(2)+1:temp_ticks_decisionPeriodA(end);
                %y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end            
                                   
            %plot((temp_ticks_decisionPeriodA(1))*[1 1],[min_sample_decisionPeriodA,max_sample_decisionPeriodA],...
            %   '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            %hold on
            %plot((temp_ticks_decisionPeriodA(2)-0.25)*[1 1],[min_sample_decisionPeriodA,max_sample_decisionPeriodA],...
            %   '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            %hold on            
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                            
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                                        
            xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);            
            ylim([min_sample_decisionPeriodA-(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.05,...
                max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);                        
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])            
            %xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end)]);            
            %xticklabels({'T1','Delay1','ChoiceCue'});
            %xticklabels({'LastT','Delay1','Decision'});
            xticklabels({'LastT','Delay-on','Decision'});            
            xtickangle(0);
            
            temptemp_range = linspace(0,max_sample_decisionPeriodA,3);
            
            temptemp_step = 2*floor((temptemp_range(2)-temptemp_range(1))*1)/1;            
            %text(temp_ticks_decisionPeriodA(2),min_sample_decisionPeriodA,'/','fontsize',12);            
            text(temp_ticks_decisionPeriodA(3),min_sample_decisionPeriodA,'/','fontsize',12);            
            
            if cellIndex_suite2p_temptempRaw == 19
                yticks([0 2]);
            end
            if cellIndex_suite2p_temptempRaw == 11
                yticks([0 1]);
            end
            if cellIndex_suite2p_temptempRaw == 121
                yticks([0 3]);
            end
            if cellIndex_suite2p_temptempRaw == 25
                %yticks([0 0.7]);
                yticks([0 0.7]);
            end
                
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('dF/F'), 'FontSize', 8);            
        end
        
        
    end        
    
    
    %if if_plot_multiPanelB4 == 1
    if false
                        
        if cellIndex_suite2p_temptempRaw==19
            max_sample_decisionPeriodA = 2;
        end
        if cellIndex_suite2p_temptempRaw==25
            max_sample_decisionPeriodA = 0.8;
        end
        
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04 127*1.05*0.95*0.93*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+0 233.4*0.93 176.7*0.86]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 233.4*0.93*0.87*1.06*1.05*1.05 176.7*0.86*1.04]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        t = tiledlayout(2,10,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>       
        
        nexttile([1 9])
        if true
            plot_lengthFlag = 1;
            
            cellID_F_dff_lengthx_sample_location = cellID_F_dff_lastT_sample_location;
            cellID_F_dff_lengthx_decisionPeriodA_location = cellID_F_dff_lastT_decisionPeriodA_location;
            
            cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_lastT_sample_location_mean;
            cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_lastT_sample_location_sem;
            cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_lastT_decisionPeriodA_location_mean;
            cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_lastT_decisionPeriodA_location_sem;
            lengthx_sample_interval = length1_sample_interval;
            
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            
            
            cellID_F_dff_lengthx_sample_location;
            cellID_F_dff_lengthx_sample_location_collapsed = [];
            cellID_F_dff_lengthx_decisionPeriodA_location;
            cellID_F_dff_lengthx_decisionPeriodA_location_collapsed = [];
            for tempi=1:numFrames
                cellID_F_dff_lengthx_sample_location_collapsed = [cellID_F_dff_lengthx_sample_location_collapsed; cellID_F_dff_lengthx_sample_location{tempi}];
                cellID_F_dff_lengthx_decisionPeriodA_location_collapsed = [cellID_F_dff_lengthx_decisionPeriodA_location_collapsed; cellID_F_dff_lengthx_decisionPeriodA_location{tempi}];
            end
            
            C1 = cellID_F_dff_lengthx_sample_location_collapsed;
            [C1_min,C1_max] = bounds(C1(:));
            
            C2 = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            [C2_min,C2_max] = bounds(C2(:));
            
            C_min = min([C1_min,C2_min]);
            C_min_valid = 0;
            
            C_max = max([C1_max,C2_max]);
            %C_max_valid = C_max*0.65;%0.5,0.65,0.55
            C_max_valid = max_sample_decisionPeriodA*2;
            if cellIndex_suite2p_temptempRaw==11 || cellIndex_suite2p_temptempRaw==19
                %C_max_valid = 3.3;%5.5,2.5,3.3,2.3
                C_max_valid = 2.5;
            end
            
            temp_trialNum = nan(1,numFrames);
            for tempi=1:numFrames
                temp_trialNum(tempi) =  size(cellID_F_dff_lengthx_sample_location{tempi},1);
            end
            temp_lengthx_locationIndex = nan(1,numFrames+1);
            for tempi=1:(numFrames+1)
                temp_lengthx_locationIndex(tempi) = sum(temp_trialNum(1:tempi-1))+1;
            end
            
            
            x = [temp_ticks_sample(1) temp_ticks_sample(end)];
            y = [1,size(cellID_F_dff_lengthx_sample_location_collapsed,1)];
            C = cellID_F_dff_lengthx_sample_location_collapsed;
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(2)];
            x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];            
            %x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]-0.5;
            y = [1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)];
            %C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed(:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %plot((temp_ticks_decisionPeriodA(1)-0.5)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
            %    '-','LineWidth',2,'Color',[1 1 1]);%3
            plot((temp_ticks_decisionPeriodA(1)-0.25)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
              '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            hold on
            
            plot((temp_ticks_decisionPeriodA(2)+0.25)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
              '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            hold on
            

            
            for tempi=1:numFrames                
            %if true
                
                temp_barLength = 3.5;%4
                
                x = [temp_ticks_decisionPeriodA(3)+0.5 temp_lengthx_locationIndex(tempi)-0.5;...
                    temp_ticks_decisionPeriodA(3)+0.5 temp_lengthx_locationIndex(tempi+1)-1.5;...
                    temp_ticks_decisionPeriodA(3)+temp_barLength temp_lengthx_locationIndex(tempi+1)-1.5;...
                    temp_ticks_decisionPeriodA(3)+temp_barLength temp_lengthx_locationIndex(tempi)-0.5];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',multi_rgbColor(tempi,:),'FaceAlpha',1,'EdgeColor','none');%0.1
                hold on     
                
                %text(temp_ticks_decisionPeriodA(2)+temp_barLength+0.25,...
                %    -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
                %    sprintf('Loc%d',tempi),'fontsize',6.5,'color',multi_rgbColor(tempi,:));
                
                temp_colorStr = ['Loc','\color[rgb]',sprintf('{%f %f %f}',multi_rgbColor(tempi,:)),sprintf('%d',tempi)];
                
                text(temp_ticks_decisionPeriodA(3)+temp_barLength+0.25,...
                   -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
                   temp_colorStr,'fontsize',6.5);

            end
                        
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);            
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);            
            xlim([temp_ticks_decisionPeriodA(1) length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);                        
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            ylim([1-1 sum(temp_trialNum)+1]);
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticklabels('');
            
            %yticks([0 400]);
            yticks([]);
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('Trials'), 'FontSize', 8);
            %colorbar;
            
            %c = colorbar('westoutside','FontSize',7);
            c = colorbar('FontSize',6.5);%9
            c.Position = c.Position+[-0.065 0.080 -0.038 -0.1];
            
            if C_max_valid > 1
                c.Ticks = [0 floor(C_max_valid)];
            else
                c.Ticks = [0 C_max_valid];
            end
            
            if if_colormap_loadEnhanced == 1
                load('parula_enhanced');
                colormap(parula_enhanced); 
                %load('gray_enhanced');
                %colormap(gray_enhanced);
            elseif if_colormap_loadEnhanced == 0
                colormap parula
                %colormap gray
            end            
            
            
        end
        
        
        
        nexttile([1 9])
        
        temp_LineWidth = 0.5;%1
        
        if true
            plot_lengthFlag = 1;
            
            cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_lastT_sample_location_mean;
            cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_lastT_sample_location_sem;
            cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_lastT_decisionPeriodA_location_mean;
            cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_lastT_decisionPeriodA_location_sem;
            lengthx_sample_interval = length1_sample_interval;
            
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            
            
            %% Plot sample
            h_line = [];
            for tempi=1:numFrames
                x = temp_ticks_sample(1):temp_ticks_sample(end);
                y = cellID_F_dff_lengthx_sample_location_mean(tempi,:);
                %h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                h_line = [h_line plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                hold on
                y_sem = cellID_F_dff_lengthx_sample_location_sem(tempi,:);
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            
            for tempi=1:plot_lengthFlag
                x = [temp_ticks_sample(1+tempi-1) min_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1) max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 min_sample_decisionPeriodA];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            %% Plot dicisionPeriodA
            for tempi=1:numFrames
                %x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);
                %y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,decisionPeriodA_interval(1):decisionPeriodA_interval(2));                
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            for tempi=1:numFrames
                x = temp_ticks_decisionPeriodA(2)+1:temp_ticks_decisionPeriodA(end);
                %y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end            
                                   
            %plot((temp_ticks_decisionPeriodA(1))*[1 1],[min_sample_decisionPeriodA,max_sample_decisionPeriodA],...
            %   '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            %hold on
            %plot((temp_ticks_decisionPeriodA(2)-0.25)*[1 1],[min_sample_decisionPeriodA,max_sample_decisionPeriodA],...
            %   '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            %hold on            
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                            
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                                        
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);
            xlim([temp_ticks_decisionPeriodA(1) length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);                        
            
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);            
            ylim([min_sample_decisionPeriodA-(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.05,...
                max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);                        
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])            
            %xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end)]);            
            %xticklabels({'T1','Delay1','ChoiceCue'});
            xticklabels({'LastT','Delay1','Decision'});
            xtickangle(0);
            
            temptemp_range = linspace(0,max_sample_decisionPeriodA,3);
            
            temptemp_step = 2*floor((temptemp_range(2)-temptemp_range(1))*1)/1;            
            %text(temp_ticks_decisionPeriodA(2),min_sample_decisionPeriodA,'/','fontsize',12);            
            text(temp_ticks_decisionPeriodA(3),min_sample_decisionPeriodA,'/','fontsize',12);            
            
            if cellIndex_suite2p_temptempRaw == 19
                yticks([0 2]);
            end
            if cellIndex_suite2p_temptempRaw == 11
                yticks([0 1]);
            end
            if cellIndex_suite2p_temptempRaw == 121
                yticks([0 3]);
            end
            if cellIndex_suite2p_temptempRaw == 25
                yticks([0 0.7]);
            end
            if cellIndex_suite2p_temptempRaw == 92
                yticks([0 0.3]);
            end            
                
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('dF/F'), 'FontSize', 8);            
        end
        
        
    end    
    if if_plot_multiPanelB4 == 1
                        
        if cellIndex_suite2p_temptempRaw==19
            max_sample_decisionPeriodA = 2;
        end
        if cellIndex_suite2p_temptempRaw==25
            max_sample_decisionPeriodA = 0.8;
        end
        if cellIndex_suite2p_temptempRaw == 2
            max_sample_decisionPeriodA = 2;
        end
        if cellIndex_suite2p_temptempRaw == 124 % 0510A
            max_sample_decisionPeriodA = 1.1;
            min_sample_decisionPeriodA = 0.05;
        end
        if cellIndex_suite2p_temptempRaw == 11
            max_sample_decisionPeriodA = 1.3;
        end       
        if cellIndex_suite2p_temptempRaw == 197 % 0510A
            max_sample_decisionPeriodA = 0.45;
            min_sample_decisionPeriodA = 0.03;
        end         
        if cellIndex_suite2p_temptempRaw == 50 % 0614A
            max_sample_decisionPeriodA = 0.55;
            min_sample_decisionPeriodA = 0.05;
        end    
        if cellIndex_suite2p_temptempRaw == 17 % 0312A_Z
            max_sample_decisionPeriodA = 1.3;
            min_sample_decisionPeriodA = 0.05;
        end            
        if cellIndex_suite2p_temptempRaw == 19 && if_monkey_D0_Z1 == 1 % 0312A_Z
            max_sample_decisionPeriodA = 1.2;
            min_sample_decisionPeriodA = 0.05;
        end   
        if cellIndex_suite2p_temptempRaw == 101 && if_monkey_D0_Z1 == 1  % 0312A_Z
            max_sample_decisionPeriodA = 0.6;
            min_sample_decisionPeriodA = 0.05;
        end   
        if cellIndex_suite2p_temptempRaw == 29 && if_monkey_D0_Z1 == 1  % 0124A_Z
            max_sample_decisionPeriodA = 1.6;
            min_sample_decisionPeriodA = 0.1;
        end        
        
        
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.04 127*1.05*0.95*0.93*1.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+0 233.4*0.93 176.7*0.86]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+300 233.4*0.93*0.87*1.06*1.05*1.05 176.7*0.86*1.04]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+300 220.7*0.80*0.96 158.0]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+300 188.6 159.0]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+300 188.6*0.92 159.0]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        %set(gcf,'Position',[35+0 42+300 188.6*0.92 159.0*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+300 188.6*0.92 159.0*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        t = tiledlayout(2,5,'TileSpacing','none','Padding','tight'); %#ok<*NASGU>       
        
        nexttile([1 4])
        if true
            plot_lengthFlag = 1;
            
            cellID_F_dff_lengthx_sample_location = cellID_F_dff_lastT_sample_location;
            cellID_F_dff_lengthx_decisionPeriodA_location = cellID_F_dff_lastT_decisionPeriodA_location;
            
            cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_lastT_sample_location_mean;
            cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_lastT_sample_location_sem;
            cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_lastT_decisionPeriodA_location_mean;
            cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_lastT_decisionPeriodA_location_sem;
            lengthx_sample_interval = length1_sample_interval;
            
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            
            
            cellID_F_dff_lengthx_sample_location;
            cellID_F_dff_lengthx_sample_location_collapsed = [];
            cellID_F_dff_lengthx_decisionPeriodA_location;
            cellID_F_dff_lengthx_decisionPeriodA_location_collapsed = [];
            for tempi=1:numFrames
                cellID_F_dff_lengthx_sample_location_collapsed = [cellID_F_dff_lengthx_sample_location_collapsed; cellID_F_dff_lengthx_sample_location{tempi}];
                cellID_F_dff_lengthx_decisionPeriodA_location_collapsed = [cellID_F_dff_lengthx_decisionPeriodA_location_collapsed; cellID_F_dff_lengthx_decisionPeriodA_location{tempi}];
            end
            
            C1 = cellID_F_dff_lengthx_sample_location_collapsed;
            [C1_min,C1_max] = bounds(C1(:));
            
            C2 = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            [C2_min,C2_max] = bounds(C2(:));
            
            C_min = min([C1_min,C2_min]);
            C_min_valid = 0;
            
            C_max = max([C1_max,C2_max]);
            %C_max_valid = C_max*0.65;%0.5,0.65,0.55
            C_max_valid = max_sample_decisionPeriodA*2;
            if cellIndex_suite2p_temptempRaw==11 || cellIndex_suite2p_temptempRaw==19
                %C_max_valid = 3.3;%5.5,2.5,3.3,2.3
                C_max_valid = 2.5;
            end
            if cellIndex_suite2p_temptempRaw == 2
                C_max_valid = 2.9;%3
            end
            if cellIndex_suite2p_temptempRaw == 20
                C_max_valid = 2.05;%2
            end      
            if cellIndex_suite2p_temptempRaw == 124 % 0510A
                C_max_valid = 0.69;
                C_min_valid = 0;
            end           
            if cellIndex_suite2p_temptempRaw == 29 && if_monkey_D0_Z1 == 1  % 0124A_Z
                C_max_valid = 1.7;
                C_min_valid = 0;
            end
            
            temp_trialNum = nan(1,numFrames);
            for tempi=1:numFrames
                temp_trialNum(tempi) =  size(cellID_F_dff_lengthx_sample_location{tempi},1);
            end
            temp_lengthx_locationIndex = nan(1,numFrames+1);
            for tempi=1:(numFrames+1)
                temp_lengthx_locationIndex(tempi) = sum(temp_trialNum(1:tempi-1))+1;
            end
            
            
            x = [temp_ticks_sample(1) temp_ticks_sample(end)];
            y = [1,size(cellID_F_dff_lengthx_sample_location_collapsed,1)];
            C = cellID_F_dff_lengthx_sample_location_collapsed;
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(2)];
            x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];            
            %x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]-0.5;
            y = [1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)];
            %C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed;
            C = cellID_F_dff_lengthx_decisionPeriodA_location_collapsed(:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
            %imagesc(x,y,C);
            imagesc(x,y,C,[C_min_valid,C_max_valid]);
            hold on
            
            %plot((temp_ticks_decisionPeriodA(1)-0.5)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
            %    '-','LineWidth',2,'Color',[1 1 1]);%3
            plot((temp_ticks_decisionPeriodA(1)-0.25)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
              '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            hold on
            
            plot((temp_ticks_decisionPeriodA(2)+0.25)*[1 1],[1,size(cellID_F_dff_lengthx_decisionPeriodA_location_collapsed,1)],...
              '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            hold on
            

            if true
                for tempi=1:numFrames
                    %if true
                    
                    temp_barLength = 3.5;%4
                    
                    x = [temp_ticks_decisionPeriodA(3)+0.5-3.5 temp_lengthx_locationIndex(tempi)-0.5;...
                        temp_ticks_decisionPeriodA(3)+0.5-3.5 temp_lengthx_locationIndex(tempi+1)-1.5;...
                        temp_ticks_decisionPeriodA(3)+temp_barLength-3.5 temp_lengthx_locationIndex(tempi+1)-1.5;...
                        temp_ticks_decisionPeriodA(3)+temp_barLength-3.5 temp_lengthx_locationIndex(tempi)-0.5];
                    y = [1 2 3 4];
                    patch('Faces',y,'Vertices',x,'FaceColor',multi_rgbColor(tempi,:),'FaceAlpha',1,'EdgeColor','none');%0.1
                    hold on
                                        
                    %temp_colorStr = ['Loc','\color[rgb]',sprintf('{%f %f %f}',multi_rgbColor(tempi,:)),sprintf('%d',tempi)];
                    %temp_colorStr = ['\color[rgb]',sprintf('{%f %f %f}',multi_rgbColor(tempi,:)),sprintf('%d',tempi)];
                    temp_colorStr = ['\color[rgb]',sprintf('{%f %f %f}',[1 1 1]*0),sprintf('%d',tempi)];
                    
                    text(temp_ticks_decisionPeriodA(3)+temp_barLength+0.25-6,...
                        -20+(temp_lengthx_locationIndex(tempi)+temp_lengthx_locationIndex(tempi+1))/2,...
                        temp_colorStr,'fontsize',6.5);
                    
                end
            end
                        
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);            
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);            
            %xlim([temp_ticks_decisionPeriodA(1) length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);                        
            xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(3)]);                                    
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            ylim([1-1 sum(temp_trialNum)+1]);
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticklabels('');
            
            %yticks([0 400]);
            yticks([]);
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('Trials'), 'FontSize', 8);
            %colorbar;
            
            %c = colorbar('westoutside','FontSize',7);
            c = colorbar('FontSize',6.5);%9
            %c.Position = c.Position+[-0.025 0.080 -0.038 -0.1];
            %c.Position = c.Position+[-0.01 0.06 -0.038 -0.075];
            c.Position = c.Position+[-0.024 0.06 -0.046 -0.075];
            
            
            if C_max_valid > 1
                c.Ticks = [0 floor(C_max_valid)];
            else
                c.Ticks = [0 floor(C_max_valid*10)/10];
            end
            
            if if_colormap_loadEnhanced == 1
                %load('parula_enhanced');
                %colormap(parula_enhanced);
                %load('gray_enhanced');
                %colormap(gray_enhanced);
                
                %load('gray_reversed_enhanced');
                %colormap(gray_reversed_enhanced);
                
                load('coolwarm_enhanced');
                colormap(coolwarm_enhanced);
                
                %load('coolwarm_n11n_enhanced');
                %colormap(coolwarm_n11n_enhanced);

            elseif if_colormap_loadEnhanced == 0
                %colormap parula
                %colormap gray
                
                %temp1 = gray;
                %temp1 = temp1(end:-1:1,:);
                %colormap(temp1);
                
                colormap(coolwarm());
                
                %temp1 = coolwarm(300);
                %temp2 = ((300-256)/2)+1;
                %temp3 = temp1(temp2:temp2+255,:);
                %colormap(temp3);
            end

            
            
        end
        
        
        
        nexttile([1 4])
        
        temp_LineWidth = 0.5;%1
        
        if true
            plot_lengthFlag = 1;
            
            cellID_F_dff_lengthx_sample_location_mean = cellID_F_dff_lastT_sample_location_mean;
            cellID_F_dff_lengthx_sample_location_sem = cellID_F_dff_lastT_sample_location_sem;
            cellID_F_dff_lengthx_decisionPeriodA_location_mean = cellID_F_dff_lastT_decisionPeriodA_location_mean;
            cellID_F_dff_lengthx_decisionPeriodA_location_sem = cellID_F_dff_lastT_decisionPeriodA_location_sem;
            lengthx_sample_interval = length1_sample_interval;
            
            temp_ticks_sample = lengthx_sample_interval;
            temp_ticks_decisionPeriodA = temp_ticks_sample(end)+periodSkipInterval+decisionPeriodA_interval;
            
            
            %% Plot sample
            h_line = [];
            for tempi=1:numFrames
                x = temp_ticks_sample(1):temp_ticks_sample(end);
                y = cellID_F_dff_lengthx_sample_location_mean(tempi,:);
                %h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                h_line = [h_line plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:))]; %#ok<*AGROW>
                hold on
                y_sem = cellID_F_dff_lengthx_sample_location_sem(tempi,:);
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            
            for tempi=1:plot_lengthFlag
                x = [temp_ticks_sample(1+tempi-1) min_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1) max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 max_sample_decisionPeriodA;...
                    temp_ticks_sample(1+tempi-1)+6 min_sample_decisionPeriodA];
                y = [1 2 3 4];
                patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.25,'EdgeColor','none');%0.1
                hold on
            end
            
            %% Plot dicisionPeriodA
            for tempi=1:numFrames
                %x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(2);
                %y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,decisionPeriodA_interval(1):decisionPeriodA_interval(2));                
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end
            for tempi=1:numFrames
                x = temp_ticks_decisionPeriodA(2)+1:temp_ticks_decisionPeriodA(end);
                %y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,1:length(x));
                y = cellID_F_dff_lengthx_decisionPeriodA_location_mean(tempi,decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
                %plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                plot(x,y,'-','LineWidth',temp_LineWidth,'Color',multi_rgbColor(tempi,:));
                hold on
                y_sem = cellID_F_dff_lengthx_decisionPeriodA_location_sem(tempi,decisionPeriodA_interval(2)+1:decisionPeriodA_interval(end));
                patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',temp_FaceAlpha,'EdgeColor',temp_EdgeColor);
                hold on
            end            
                                   
            %plot((temp_ticks_decisionPeriodA(1))*[1 1],[min_sample_decisionPeriodA,max_sample_decisionPeriodA],...
            %   '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            %hold on
            %plot((temp_ticks_decisionPeriodA(2)-0.25)*[1 1],[min_sample_decisionPeriodA,max_sample_decisionPeriodA],...
            %   '-','LineWidth',1,'Color',[1 1 1]);%3,2           
            %hold on            
            
            set(gca,'linewidth',1.5);
            set(gca,'color',backgounrdColor);
            %xlim([1-xlim_padSize1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                            
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(2)+xlim_padSize1+2]);                                                                        
            %xlim([1 length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);
            %xlim([temp_ticks_decisionPeriodA(1) length1_sample_interval(end)+periodSkipInterval+decisionPeriodA_interval(3)+xlim_padSize1+2+2]);                        
            xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(3)]);                                    
            
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA]);
            %ylim([min_sample_decisionPeriodA,max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);            
            ylim([min_sample_decisionPeriodA-(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.05,...
                max_sample_decisionPeriodA+(max_sample_decisionPeriodA-min_sample_decisionPeriodA)*0.125]);                        
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])            
            %xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end-1)]);
            xticks([temp_ticks_sample(1:end-1) temp_ticks_decisionPeriodA(1:end)]);            
            %xticklabels({'T1','Delay1','ChoiceCue'});
            %xticklabels({'LastT','Delay1','Decision'});
            xticklabels({'LastT','Delay-on','Decision'});            
            xtickangle(0);
            
            temptemp_range = linspace(0,max_sample_decisionPeriodA,3);
            
            temptemp_step = 2*floor((temptemp_range(2)-temptemp_range(1))*1)/1;            
            %text(temp_ticks_decisionPeriodA(2),min_sample_decisionPeriodA,'/','fontsize',12);            
            %text(temp_ticks_decisionPeriodA(3),min_sample_decisionPeriodA,'/','fontsize',12);            
            
            if cellIndex_suite2p_temptempRaw == 19 && if_monkey_D0_Z1 == 0
                yticks([0 2]);
            end
            if cellIndex_suite2p_temptempRaw == 11
                yticks([0 1]);
            end
            if cellIndex_suite2p_temptempRaw == 121
                yticks([0 3]);
            end
            if cellIndex_suite2p_temptempRaw == 25
                yticks([0 0.7]);
            end
            if cellIndex_suite2p_temptempRaw == 92
                yticks([0 0.3]);
            end    
            
            if cellIndex_suite2p_temptempRaw == 2
                yticks([0 1.9]);
            end            
                
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            ylabel(sprintf('dF/F'), 'FontSize', 8);            
        end
        
        
    end        
    
    
    
    
    if if_plot_multiPanelC == 1
        fig1 = figure('Name',sprintf('Cell id %d, FOV id %s',cellIndex_suite2p_temptempRaw,FOVName_currentFOV2),'NumberTitle','off');
        %set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67 127*1.05*0.95*0.93*0.6]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
        set(gcf,'Position',[35+0 42+0 381*0.9*0.977*0.67*1.2 127*1.05*0.95*0.93*0.6*1.2*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                        
                
        t = tiledlayout(4,10,'TileSpacing','none','Padding','loose'); %#ok<*NASGU>

        %temp1 = locDistri_confusionMatrix_B(:,:,1)';                        
        %temp1 = temp1(:,1:16);
        
        temp1 = cellID_F_dff_decisionBin1_seqMerged_eachLoc;                
        
        %nexttile
        %nexttile([5 9])
        nexttile([3 10])
        
        C = temp1;
        
        
        temp_repeatBoolIndex = false(size(C));
        if if_seqIdentity_default0_tunedLoc1 == 1
            TunedLoc;
            
            for tempi=1:size(C,1) % loc
                for tempj=1:size(C,2) % seq
                    tempBoolIndex_loc = boolIndex_location_seq_T(temp_seqIndex_eachLoc(tempi,tempj),:);
                    
                    if ~ismember(tempi,TunedLoc)
                        if sum(tempBoolIndex_loc(TunedLoc)) > 0
                            
                            %C(tempi,tempj) = nan;
                            C(tempi,tempj) = 0;
                            
                            temp_repeatBoolIndex(tempi,tempj) = true;                            
                        end
                    end
                    
                end
            end
        end
        
        
        
        
        %C_max = 4;
        C_max = max(C,[],'all');
        
        C_max = max([C_max,1]);
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[0 C_max]);
        hold on
        
        % Plot length bound in response
        for temptempi=1:size(C,2)
            
            if temptempi==1 || temptempi==6
                plot([temptempi temptempi]+0.5,[0.5 size(C,1)+0.5], 'color', [1 1 1]);
                hold on
            end
        end

        %         for tempi=1:size(C,1)
        %             for tempj=1:size(C,2)
        %                 if ~isnan(C(tempi,tempj))
        %                     continue
        %                 end
        %                 plot([tempj-0.4,tempj+0.4],[tempi-0.4,tempi+0.4],...
        %                     '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
        %                 hold on
        %             end
        %         end
    
        for tempi=1:size(C,1) % loc
            for tempj=1:size(C,2) % seq
                if isnan(C(tempi,tempj))
                    plot([tempj-0.4,tempj+0.4],[tempi-0.4,tempi+0.4],...
                        '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
                    hold on
                elseif temp_repeatBoolIndex(tempi,tempj) == true
                    %                     plot([tempj-0.4,tempj+0.4],[tempi-0.4,tempi+0.4],...
                    %                         '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
                    %                     hold on
                    %
                    %                     plot([tempj+0.4,tempj-0.4],[tempi-0.4,tempi+0.4],...
                    %                         '-','LineWidth',1,'Color',[1 1 1]*0.7); %#ok<*UNRCH>
                    %                     hold on
                    
                else
                    if if_seqIdentity_default0_tunedLoc1 == 0
                        temptempBoolIndex = C(tempi,tempj)>prctile(C,50,'all');
                    elseif if_seqIdentity_default0_tunedLoc1 == 1
                        %temptempBoolIndex = C(tempi,tempj)>prctile(C,50,'all');
                        temptempBoolIndex = (C(tempi,tempj)>prctile(C,50,'all')) | (ismember(tempi,TunedLoc));
                    end     
                    
                    %if temptempBoolIndex == true
                        %temp_txt = '566';
                        temp_txt = num2str(seqSet_inOne_inOne(temp_seqIndex_eachLoc(tempi,tempj)));
                        text(tempj+0.05,tempi-0.05,temp_txt,...
                            'color',[1 1 1]*0,'HorizontalAlignment','center','rotation',0,'fontsize',6);%6
                    %end
                end                
            end
        end
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 8) %12
        
        
        
    tempLabelStr_res = string(1:3);
    tempLabelStr_res(1) = "Len1";
    tempLabelStr_res(2) = "Len2";
    tempLabelStr_res(3) = "Len3";

    xtl=string(tempLabelStr_res);
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    %xtext_xp=xt;
    xtext_xp = nan(1,3);
    %for tempi=1:length(xtext_xp)
    %    xtext_xp(tempi) = sum(numSeq(1:(tempi-1))) + numSeq(tempi)/2;
    %end
    
    xtext_xp(1) = 0.6;
    xtext_xp(2) = 3.5;
    xtext_xp(3) = 10.5;
    
    %xtext_xp = [xtext_xp xtext_xp(end)+20/2];
    xtext_xp = xtext_xp+0.5;
    
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(tempLabelStr_res))+0.75+0.30+0.3;
    % 设置ttext的x坐标位置
    %ytext_xp=(xt(1))*ones(1,length(tempLabelStr_res))-1.6;%-0.75
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',8);%8
    set(gca,'xticklabel','');

    
        %set(gca,'xticklabel','');
        
        %yticklabels('');
        
        set(gca,'YTick',[1:1:numFrames]); %#ok<*NBRAK>
        
        xtickangle(0);
        
        set(gca,'TickLength',[0 0]);
        
        ylabel('Location', 'FontSize', 8, 'FontWeight', 'normal');
        
        %c = colorbar('FontSize',6.5);%9
        c = colorbar('southoutside','FontSize',6.5);%9        
        %c.Position = c.Position+[-0.036 0.065 -0.038 -0.1];
        %c.Position = c.Position+[-0.029 0.065 -0.038 -0.1];   
        
        c.Position = c.Position+[0.61 0.16 -0.62 -0.07];   
        
        %c.Ticks = [0 C_max];
        c.Ticks = [0 floor(C_max)];
        
        if if_colormap_loadEnhanced == 1
            load('parula_enhanced');
            colormap(parula_enhanced);
            %load('gray_enhanced');
            %colormap(gray_enhanced);
        elseif if_colormap_loadEnhanced == 0
            colormap parula
            %colormap gray
        end
        
        %temp_title = title(sprintf('Behavior model'),'FontSize',11);
        %temp_title.Interpreter = 'none';
        
        
        
    end    
    
    
end



if false
% cellID_F_dff_length1_sample = squeeze(F_dff_length1_sample(temp_id2,:,:));
% cellID_F_dff_length2_sample = squeeze(F_dff_length2_sample(temp_id2,:,:));
% cellID_F_dff_length3_sample = squeeze(F_dff_length3_sample(temp_id2,:,:));


F_dff_decisionBin_demo = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);

temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
F_dff_baselineBin_demo = mean(F_dff_baselinePeriod(:,:,temptemp_range),3);


decodingDataSimplified.F_dff_decisionBin_demo = F_dff_decisionBin_demo;
decodingDataSimplified.F_dff_baselineBin_demo = F_dff_baselineBin_demo;


temp_F_dff_exampleNeuron_memoryX = struct;

temp_F_dff_exampleNeuron_memoryX.F_dff_length1_sample_exampleNeuron_memory = cellID_F_dff_length1_sample;
temp_F_dff_exampleNeuron_memoryX.F_dff_length2_sample_exampleNeuron_memory = cellID_F_dff_length2_sample;
temp_F_dff_exampleNeuron_memoryX.F_dff_length3_sample_exampleNeuron_memory = cellID_F_dff_length3_sample;
temp_F_dff_exampleNeuron_memoryX.F_dff_decisionPeriodA_exampleNeuron_memory = cellID_F_dff_decisionPeriodA;


% F_dff_exampleNeuron_memory1 = temp_F_dff_exampleNeuron_memoryX;
% decodingDataSimplified.F_dff_exampleNeuron_memory1 = F_dff_exampleNeuron_memory1;

F_dff_exampleNeuron_memory2 = temp_F_dff_exampleNeuron_memoryX;
decodingDataSimplified.F_dff_exampleNeuron_memory2 = F_dff_exampleNeuron_memory2;

end

%% End