function [selectiveCellBoolIndex,selectiveCellNum] = ...
    frameAnova1_v5(dff,trialLabel,significantThreshold,significantThreshold_lowBound,temp_target_seqSet,if_selective_seq0_loc1)

if exist('temp_target_seqSet','var') == 1
    
    seq_length = length(temp_target_seqSet{1});
    if if_selective_seq0_loc1 == 1 && seq_length > 1
        dff_raw = dff; %#ok<*NASGU>
        trialLabel_raw = trialLabel;
        
        [roi_num,trial_num,frame_num] = size(dff);
        
        
        dff_loc = zeros(roi_num, trial_num*seq_length, frame_num, 'single');
        trialLabel_loc = zeros(1, trial_num*seq_length);
        for tempi=1:trial_num
            temp_seq = temp_target_seqSet{trialLabel_raw(tempi)};
            temp_dff = dff_raw(:,tempi,:);
            for tempj=1:seq_length
                dff_loc(:,(tempi-1)*seq_length+tempj,:) = temp_dff;
                trialLabel_loc((tempi-1)*seq_length+tempj) = temp_seq(tempj);
            end
        end
        
        dff = dff_loc;
        trialLabel = trialLabel_loc;
        
        %dff = dff_raw;
        %trialLabel = trialLabel_raw;
    end
    
end

p_anova = zeros(size(dff,1),size(dff,3));
parfor tempi=1:size(dff,1)
    temp_p = zeros(1,size(dff,3));
    temp_dff = squeeze(dff(tempi,:,:));
    for tempj=1:size(dff,3)
        [temp_p(tempj),~,~] = anova1(temp_dff(:,tempj), trialLabel,'off');
    end
    p_anova(tempi,:) = temp_p;
end

min_continueFrame = 6;%3-->4-->5-->6

selectiveCellBoolIndex_raw = ...
    (p_anova<significantThreshold) & (p_anova>significantThreshold_lowBound);

selectiveCellBoolIndex_raw2 = selectiveCellBoolIndex_raw;

%% filter selectiveCellBoolIndex_raw
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
            %if temp_index2 < length(temp_boolIndex)
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
            %end            
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
%find(selectiveCellBoolIndex==true)'
% sum_raw2 = sum(selectiveCellBoolIndex_raw2,2);
% sum_raw = sum(selectiveCellBoolIndex_raw,2);
% find(sum_raw2>0)'
% find(sum_raw>0)'

%% others
selectiveCellBoolIndex_percent = sum(selectiveCellBoolIndex_raw,2)./size(selectiveCellBoolIndex_raw,2);

temp1 = selectiveCellBoolIndex_percent(selectiveCellBoolIndex_percent>0);
temp1_median = median(temp1);

[idx,C] = kmeans(selectiveCellBoolIndex_percent,2,'Replicates',10);
[~,I] = max(C);
[temp_min,~] = bounds(selectiveCellBoolIndex_percent(idx==I));
temp_threshold_raw = temp_min;

temp_threshold = mean([temp_threshold_raw,temp1_median]);

% selectiveCellBoolIndex = selectiveCellBoolIndex_percent > temp_threshold;
selectiveCellBoolIndex = selectiveCellBoolIndex_percent >= temp_threshold;
selectiveCellNum = sum(selectiveCellBoolIndex);


if false
    h=histogram(selectiveCellBoolIndex_percent); %#ok<*UNRCH>
    hold on
    plot(temp_threshold_raw*[1 1],[min(h.Values),max(h.Values)],'Color','k');
    hold on
    plot(temp_threshold*[1 1],[min(h.Values),max(h.Values)],'Color','r');
    hold on
end

%% End