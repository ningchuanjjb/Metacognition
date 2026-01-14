function rProb_glm_output = fun_rProb_glm_frame_v6(rProb_glm_options)
%% Initialization
% close all

min_continueFrame = 6;%3-->4-->5-->6
significantThreshold = 0.05;
significantThreshold_lowBound = -1;


offloadingProb_inOne = rProb_glm_options.offloadingProb_inOne;
precision_inOne = rProb_glm_options.precision_inOne;
numSeq = rProb_glm_options.numSeq;
target_seqSet = rProb_glm_options.target_seqSet;

temp_F_dff = rProb_glm_options.temp_F_dff;
trial_para = rProb_glm_options.trial_para;

T = size(temp_F_dff,3);

seqIndex = zeros(1,trial_para.trial_count);
for tempi=1:trial_para.trial_count
    currentSequence = trial_para.currentSequence{tempi};
    temp_seq_length = length(currentSequence);
    for tempj=1:numSeq(temp_seq_length)
        if sum(ismember(currentSequence,target_seqSet{temp_seq_length}{tempj})) == temp_seq_length
            break
        end
    end
    seqIndex(tempi) = sum(numSeq(1:temp_seq_length-1)) + tempj;
end
offloadingProb_inSeq = offloadingProb_inOne;
precision_inSeq = precision_inOne;

roiNum = size(temp_F_dff,1);


temp_F_dff_seqMerged = zeros(roiNum,sum(numSeq),T);
for tempi=1:sum(numSeq)
    temp_dff = temp_F_dff(:,seqIndex==tempi,:);
    temp_F_dff_seqMerged(:,tempi,:) = mean(temp_dff,2);
end


temp_F_dff_frame = temp_F_dff_seqMerged;
temp_F_dff_mean = mean(temp_F_dff_seqMerged,3);
temp_offloadingProb = offloadingProb_inSeq;
temp_precision = precision_inSeq;

%% figglm for mean (rProb)
x = temp_F_dff_mean';
y = temp_offloadingProb';

r2_rProb_glm = zeros(roiNum,1);
p_rProb_glm = zeros(roiNum,1);
parfor tempi=1:roiNum
    % for tempi=1:roiNum
    warning('off');
    temp_mdl = fitglm(x(:,tempi),y);
    r2_rProb_glm(tempi) = temp_mdl.Rsquared.Adjusted;
    p_rProb_glm(tempi) = temp_mdl.Coefficients.pValue(2);
    warning('on');
end
r2_rProb_glm; %#ok<*VUNUS>
p_rProb_glm;


%% figglm for mean (precision)
x = temp_F_dff_mean';
y = temp_precision';

r2_precision_glm = zeros(roiNum,1);
p_precision_glm = zeros(roiNum,1);
parfor tempi=1:roiNum
    % for tempi=1:roiNum
    warning('off');
    temp_mdl = fitglm(x(:,tempi),y);
    r2_precision_glm(tempi) = temp_mdl.Rsquared.Adjusted;
    p_precision_glm(tempi) = temp_mdl.Coefficients.pValue(2);
    warning('on');
end
r2_precision_glm; %#ok<*VUNUS>
p_precision_glm;


%% rProb
if true
    %% figglm for frame (rProb)
    temp_F_dff_frame;
    y = temp_offloadingProb';
    
    p_frame = zeros(roiNum,T);
    
    % t0 = tic;
    parfor tempi=1:roiNum
        % for tempi=1:roiNum
        for tempj=1:T
            %for tempj=1:5
            x = temp_F_dff_frame(tempi,:,tempj)';
            [~,p_frame(tempi,tempj)] = corr(x(~isnan(x)),y(~isnan(x)));
        end
    end
    p_frame;
    % toc(t0);
    
    %min_continueFrame = 6;%3-->4-->5-->6
    %significantThreshold = 0.05;
    %significantThreshold_lowBound = -1;
    
    
    selectiveCellBoolIndex_raw = ...
        (p_frame<significantThreshold) & (p_frame>significantThreshold_lowBound);
    selectiveCellBoolIndex_raw2 = selectiveCellBoolIndex_raw;
    
    %% filter selectiveCellBoolIndex_raw (rProb)
    for tempi=1:size(selectiveCellBoolIndex_raw,1)
        temp_boolIndex = selectiveCellBoolIndex_raw(tempi,:);
        if sum(temp_boolIndex) == 0
            continue
        end
        
        temp_index = find(temp_boolIndex==true);
        
        temp_eventIndex = [];
        temp_eventDuration = [];
        
        for tempj=1:length(temp_index)
            temp_index2 = temp_index(tempj);
            if temp_index2 == 1 || temp_boolIndex(temp_index2-1) == false
                temp_eventIndex = [temp_eventIndex temp_index2]; %#ok<*AGROW>
                tempk = 1;
                while true
                    if (temp_index2+tempk) > length(temp_boolIndex)
                        break
                    end
                    if temp_boolIndex(temp_index2+tempk) == true
                        tempk = tempk + 1;
                    else
                        break
                    end
                end
                temp_eventDuration = [temp_eventDuration tempk];
            end
        end
        tempBoolIndex_longFrame = temp_eventDuration>=min_continueFrame;
        temp_eventIndex_long = temp_eventIndex(tempBoolIndex_longFrame);
        temp_eventDuration_long = temp_eventDuration(tempBoolIndex_longFrame);
        
        temp_boolIndex_long = false(1,length(temp_boolIndex));
        if isempty(temp_eventIndex_long) == false
            for tempj=1:length(temp_eventIndex_long)
                temp_eventIndex_long_2 = temp_eventIndex_long(tempj);
                temp_eventDuration_long_2 = temp_eventDuration_long(tempj);
                temp_range = temp_eventIndex_long_2:(temp_eventIndex_long_2+temp_eventDuration_long_2-1);
                temp_boolIndex_long(temp_range) = true;
            end
        end
        
        selectiveCellBoolIndex_raw(tempi,:) = temp_boolIndex_long;
    end
    
    %% others (rProb)
    selectiveCellBoolIndex_percent = sum(selectiveCellBoolIndex_raw,2)./size(selectiveCellBoolIndex_raw,2);
    
    temp1 = selectiveCellBoolIndex_percent(selectiveCellBoolIndex_percent>0);
    temp1_median = median(temp1);
    
    [idx,C] = kmeans(selectiveCellBoolIndex_percent,2,'Replicates',10);
    [~,I] = max(C);
    [temp_min,~] = bounds(selectiveCellBoolIndex_percent(idx==I));
    temp_threshold_raw = temp_min;
    
    temp_threshold = mean([temp_threshold_raw,temp1_median]);
    
    % selectiveCellBoolIndex_rProb_glm = selectiveCellBoolIndex_percent > temp_threshold;
    selectiveCellBoolIndex_rProb_glm = selectiveCellBoolIndex_percent >= temp_threshold;
    selectiveCellNum_rProb_glm = sum(selectiveCellBoolIndex_rProb_glm); %#ok<*NASGU>
    
    
    if false
        h=histogram(selectiveCellBoolIndex_percent); %#ok<*UNRCH>
        hold on
        plot(temp_threshold_raw*[1 1],[min(h.Values),max(h.Values)],'Color','k');
        hold on
        plot(temp_threshold*[1 1],[min(h.Values),max(h.Values)],'Color','r');
        hold on
    end
    
end




%% precision
if true
    %% figglm for frame (precision)
    temp_F_dff_frame;
    y = temp_precision';
    
    p_frame = zeros(roiNum,T);
    
    % t0 = tic;
    parfor tempi=1:roiNum
        % for tempi=1:roiNum
        for tempj=1:T
            %for tempj=1:5
            x = temp_F_dff_frame(tempi,:,tempj)';
            [~,p_frame(tempi,tempj)] = corr(x(~isnan(x)),y(~isnan(x)));
        end
    end
    p_frame;
    % toc(t0);
    
    %min_continueFrame = 6;%3-->4-->5-->6
    %significantThreshold = 0.05;
    %significantThreshold_lowBound = -1;
    
    
    selectiveCellBoolIndex_raw = ...
        (p_frame<significantThreshold) & (p_frame>significantThreshold_lowBound);
    selectiveCellBoolIndex_raw2 = selectiveCellBoolIndex_raw;
    
    %% filter selectiveCellBoolIndex_raw (precision)
    for tempi=1:size(selectiveCellBoolIndex_raw,1)
        temp_boolIndex = selectiveCellBoolIndex_raw(tempi,:);
        if sum(temp_boolIndex) == 0
            continue
        end
        
        temp_index = find(temp_boolIndex==true);
        
        temp_eventIndex = [];
        temp_eventDuration = [];
        
        for tempj=1:length(temp_index)
            temp_index2 = temp_index(tempj);
            if temp_index2 == 1 || temp_boolIndex(temp_index2-1) == false
                temp_eventIndex = [temp_eventIndex temp_index2]; %#ok<*AGROW>
                tempk = 1;
                while true
                    if (temp_index2+tempk) > length(temp_boolIndex)
                        break
                    end
                    if temp_boolIndex(temp_index2+tempk) == true
                        tempk = tempk + 1;
                    else
                        break
                    end
                end
                temp_eventDuration = [temp_eventDuration tempk];
            end
        end
        tempBoolIndex_longFrame = temp_eventDuration>=min_continueFrame;
        temp_eventIndex_long = temp_eventIndex(tempBoolIndex_longFrame);
        temp_eventDuration_long = temp_eventDuration(tempBoolIndex_longFrame);
        
        temp_boolIndex_long = false(1,length(temp_boolIndex));
        if isempty(temp_eventIndex_long) == false
            for tempj=1:length(temp_eventIndex_long)
                temp_eventIndex_long_2 = temp_eventIndex_long(tempj);
                temp_eventDuration_long_2 = temp_eventDuration_long(tempj);
                temp_range = temp_eventIndex_long_2:(temp_eventIndex_long_2+temp_eventDuration_long_2-1);
                temp_boolIndex_long(temp_range) = true;
            end
        end
        
        selectiveCellBoolIndex_raw(tempi,:) = temp_boolIndex_long;
    end
    
    %% others (precision)
    selectiveCellBoolIndex_percent = sum(selectiveCellBoolIndex_raw,2)./size(selectiveCellBoolIndex_raw,2);
    
    temp1 = selectiveCellBoolIndex_percent(selectiveCellBoolIndex_percent>0);
    temp1_median = median(temp1);
    
    [idx,C] = kmeans(selectiveCellBoolIndex_percent,2,'Replicates',10);
    [~,I] = max(C);
    [temp_min,~] = bounds(selectiveCellBoolIndex_percent(idx==I));
    temp_threshold_raw = temp_min;
    
    temp_threshold = mean([temp_threshold_raw,temp1_median]);
    
    % selectiveCellBoolIndex_precision_glm = selectiveCellBoolIndex_percent > temp_threshold;
    selectiveCellBoolIndex_precision_glm = selectiveCellBoolIndex_percent >= temp_threshold;
    selectiveCellNum_precision_glm = sum(selectiveCellBoolIndex_precision_glm); %#ok<*NASGU>
    
    
    if false
        h=histogram(selectiveCellBoolIndex_percent); %#ok<*UNRCH>
        hold on
        plot(temp_threshold_raw*[1 1],[min(h.Values),max(h.Values)],'Color','k');
        hold on
        plot(temp_threshold*[1 1],[min(h.Values),max(h.Values)],'Color','r');
        hold on
    end
    
end

%%
rProb_glm_output = struct;
rProb_glm_output.r2_rProb_glm = r2_rProb_glm;
rProb_glm_output.p_rProb_glm = p_rProb_glm;
rProb_glm_output.selectiveCellBoolIndex_rProb_glm = selectiveCellBoolIndex_rProb_glm;
rProb_glm_output.r2_precision_glm = r2_precision_glm;
rProb_glm_output.p_precision_glm = p_precision_glm;
rProb_glm_output.selectiveCellBoolIndex_precision_glm = selectiveCellBoolIndex_precision_glm;
rProb_glm_output.temp_F_dff_mean = temp_F_dff_mean;
rProb_glm_output.temp_F_dff_frame = temp_F_dff_frame;

%% End