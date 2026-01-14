% Chuan's 15th script (20251214)
% This script: Time courses of memory strength & meta-memory.
close all;
%% Initialization

if_4typesFilter_all0_cMC1_cF2_cME3 = 0;%0

if if_4typesFilter_all0_cMC1_cF2_cME3 == 0
    tempBoolIndex_4typesFilter = true(length(choiceMemoryCorrectBoolIndex),1);
    
elseif if_4typesFilter_all0_cMC1_cF2_cME3 == 1
    tempBoolIndex_4typesFilter = choiceMemoryCorrectBoolIndex';
    
elseif if_4typesFilter_all0_cMC1_cF2_cME3 == 2
    tempBoolIndex_4typesFilter = choiceOffloadBoolIndex';
    
elseif if_4typesFilter_all0_cMC1_cF2_cME3 == 3
    tempBoolIndex_4typesFilter = choiceMemoryErrorBoolIndex';
end



%% Preparation
trialBoolIndex_memoryPrecisionLow_metaLow_choice;
trialBoolIndex_memoryPrecisionHigh_metaHigh_choice;
trialBoolIndex_memoryPrecisionLow_metaHigh_choice;
trialBoolIndex_memoryPrecisionHigh_metaLow_choice;

% if if_memeoryPrecision_stimuli0_response1 == 0
%     temp_seqIndex = seqIndex;
% elseif if_memeoryPrecision_stimuli0_response1 == 1
%     temp_seqIndex = seqIndex_response;
% end

for tempi=1:3
    temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    
    %temp1 = ismember(seqIndex,temp_range);
    %temp2 = ismember(seqIndex_response,temp_range);
    
    %temptempBoolIndex = ismember(temp_seqIndex,temp_range)';
    temptempBoolIndex = ismember(seqIndex,temp_range)';
    
    if tempi==1
        trialBoolIndex_length1 = temptempBoolIndex;
    elseif tempi==2
        trialBoolIndex_length2 = temptempBoolIndex;
    elseif tempi==3
        trialBoolIndex_length3 = temptempBoolIndex;
    end
end
trialBoolIndex_length1;
trialBoolIndex_length2;
trialBoolIndex_length3;


%% meta_trialLevel_crossTime_multiPeriod_4types, each length
for tempj=1:4
    
    if tempj==1
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
    elseif tempj==2
        temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
    elseif tempj==3
        temp_trialBoolIndex_lengthx = trialBoolIndex_length3;
    elseif tempj==4
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
    end
    temp_trialBoolIndex_lengthx = temp_trialBoolIndex_lengthx & tempBoolIndex_4typesFilter;
    
    meta_trialLevel_crossTime_multiPeriod_4types_lengthx = cell(3,4);
    meta_trialLevel_crossTimeMean_multiPeriod_4types_lengthx = cell(3,4);
    for tempi=1:3
        if tempi==1
            temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_baseline;
        elseif tempi==2
            if tempj==1
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_length1;
            elseif tempj==2
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_length2;
            elseif tempj==3
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_length3;
            elseif tempj==4
                temp_meta_trialLevel_crossTime = nan(size(meta_trialLevel_crossTime_length3));
                
                for temptempi=1:size(temp_meta_trialLevel_crossTime,1)
                    if trialBoolIndex_length1(temptempi) == true
                        temp1 = meta_trialLevel_crossTime_length1(temptempi,:);
                        temp_interp_factor = 3;
                        
                        if isnan(temp1(1)) == true
                            continue
                        end
                        
                        x = 1:length(temp1);
                        x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);
                        
                        temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');
                        
                        for temptempj=1:ceil((temp_interp_factor-1))
                            temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
                        end
                        
                        temp_meta_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length2(temptempi) == true
                        temp1 = meta_trialLevel_crossTime_length2(temptempi,:);
                        temp_interp_factor = 1.5;
                        
                        if isnan(temp1(1)) == true
                            continue
                        end
                        
                        x = 1:length(temp1);
                        x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);
                        
                        temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');
                        
                        for temptempj=1:ceil((temp_interp_factor-1))
                            temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
                        end
                        
                        temp_meta_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length3(temptempi) == true
                        temp_meta_trialLevel_crossTime(temptempi,:) = meta_trialLevel_crossTime_length3(temptempi,:);
                    end
                    
                end
                temp_meta_trialLevel_crossTime;
                
                
            end
        elseif tempi==3
            temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_delay1;
        end
        
        temp_meta_trialLevel_crossTime_highMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice & temp_trialBoolIndex_lengthx,:);
        temp_meta_trialLevel_crossTime_lowMatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaLow_choice & temp_trialBoolIndex_lengthx,:);
        temp_meta_trialLevel_crossTime_overMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaHigh_choice & temp_trialBoolIndex_lengthx,:);
        temp_meta_trialLevel_crossTime_underMismatch = temp_meta_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaLow_choice & temp_trialBoolIndex_lengthx,:);
        
        meta_trialLevel_crossTime_multiPeriod_4types_lengthx{tempi,1} = temp_meta_trialLevel_crossTime_highMatch;
        meta_trialLevel_crossTime_multiPeriod_4types_lengthx{tempi,2} = temp_meta_trialLevel_crossTime_lowMatch;
        meta_trialLevel_crossTime_multiPeriod_4types_lengthx{tempi,3} = temp_meta_trialLevel_crossTime_overMismatch;
        meta_trialLevel_crossTime_multiPeriod_4types_lengthx{tempi,4} = temp_meta_trialLevel_crossTime_underMismatch;
        
        meta_trialLevel_crossTimeMean_multiPeriod_4types_lengthx{tempi,1} = mean(temp_meta_trialLevel_crossTime_highMatch,1,'omitnan');
        meta_trialLevel_crossTimeMean_multiPeriod_4types_lengthx{tempi,2} = mean(temp_meta_trialLevel_crossTime_lowMatch,1,'omitnan');
        meta_trialLevel_crossTimeMean_multiPeriod_4types_lengthx{tempi,3} = mean(temp_meta_trialLevel_crossTime_overMismatch,1,'omitnan');
        meta_trialLevel_crossTimeMean_multiPeriod_4types_lengthx{tempi,4} = mean(temp_meta_trialLevel_crossTime_underMismatch,1,'omitnan');
        
    end
    
    if tempj==1
        meta_trialLevel_crossTime_multiPeriod_4types_length1 = meta_trialLevel_crossTime_multiPeriod_4types_lengthx;
        meta_trialLevel_crossTimeMean_multiPeriod_4types_length1 = meta_trialLevel_crossTimeMean_multiPeriod_4types_lengthx;
    elseif tempj==2
        meta_trialLevel_crossTime_multiPeriod_4types_length2 = meta_trialLevel_crossTime_multiPeriod_4types_lengthx;
        meta_trialLevel_crossTimeMean_multiPeriod_4types_length2 = meta_trialLevel_crossTimeMean_multiPeriod_4types_lengthx;
    elseif tempj==3
        meta_trialLevel_crossTime_multiPeriod_4types_length3 = meta_trialLevel_crossTime_multiPeriod_4types_lengthx;
        meta_trialLevel_crossTimeMean_multiPeriod_4types_length3 = meta_trialLevel_crossTimeMean_multiPeriod_4types_lengthx;
    elseif tempj==4
        meta_trialLevel_crossTime_multiPeriod_4types_length123 = meta_trialLevel_crossTime_multiPeriod_4types_lengthx;
        meta_trialLevel_crossTimeMean_multiPeriod_4types_length123 = meta_trialLevel_crossTimeMean_multiPeriod_4types_lengthx;
    end
end

meta_trialLevel_crossTime_multiPeriod_4types_length1;
meta_trialLevel_crossTime_multiPeriod_4types_length2;
meta_trialLevel_crossTime_multiPeriod_4types_length3;
meta_trialLevel_crossTime_multiPeriod_4types_length123;

meta_trialLevel_crossTimeMean_multiPeriod_4types_length1;
meta_trialLevel_crossTimeMean_multiPeriod_4types_length2;
meta_trialLevel_crossTimeMean_multiPeriod_4types_length3;
meta_trialLevel_crossTimeMean_multiPeriod_4types_length123;


%% precision_trialLevel_crossTime_multiPeriod_4types, each length
for tempj=1:4
    
    if tempj==1
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
    elseif tempj==2
        temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
    elseif tempj==3
        temp_trialBoolIndex_lengthx = trialBoolIndex_length3;
    elseif tempj==4
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
    end
    temp_trialBoolIndex_lengthx = temp_trialBoolIndex_lengthx & tempBoolIndex_4typesFilter;
    
    precision_trialLevel_crossTime_multiPeriod_4types_lengthx = cell(3,4);
    precision_trialLevel_crossTimeMean_multiPeriod_4types_lengthx = cell(3,4);
    for tempi=1:3
        if tempi==1
            temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_baseline;
        elseif tempi==2
            if tempj==1
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_length1;
            elseif tempj==2
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_length2;
            elseif tempj==3
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_length3;
            elseif tempj==4
                temp_memoryPrecision_trialLevel_crossTime = nan(size(memoryPrecision_trialLevel_crossTime_length3));
                
                for temptempi=1:size(temp_memoryPrecision_trialLevel_crossTime,1)
                    if trialBoolIndex_length1(temptempi) == true
                        temp1 = memoryPrecision_trialLevel_crossTime_length1(temptempi,:);
                        temp_interp_factor = 3;
                        
                        if isnan(temp1(1)) == true
                            continue
                        end
                        
                        x = 1:length(temp1);
                        x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);
                        
                        temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');
                        
                        for temptempj=1:ceil((temp_interp_factor-1))
                            temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
                        end
                        
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length2(temptempi) == true
                        temp1 = memoryPrecision_trialLevel_crossTime_length2(temptempi,:);
                        temp_interp_factor = 1.5;
                        
                        if isnan(temp1(1)) == true
                            continue
                        end
                        
                        x = 1:length(temp1);
                        x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);
                        
                        temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');
                        
                        for temptempj=1:ceil((temp_interp_factor-1))
                            temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
                        end
                        
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length3(temptempi) == true
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = memoryPrecision_trialLevel_crossTime_length3(temptempi,:);
                    end
                    
                end
                temp_memoryPrecision_trialLevel_crossTime;
                
            end
        elseif tempi==3
            temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_delay1;
        end
        
        
        temp_memoryPrecision_trialLevel_crossTime_highMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaHigh_choice & temp_trialBoolIndex_lengthx,:);
        temp_memoryPrecision_trialLevel_crossTime_lowMatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaLow_choice & temp_trialBoolIndex_lengthx,:);
        temp_memoryPrecision_trialLevel_crossTime_overMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionLow_metaHigh_choice & temp_trialBoolIndex_lengthx,:);
        temp_memoryPrecision_trialLevel_crossTime_underMismatch = temp_memoryPrecision_trialLevel_crossTime(trialBoolIndex_memoryPrecisionHigh_metaLow_choice & temp_trialBoolIndex_lengthx,:);
        
        precision_trialLevel_crossTime_multiPeriod_4types_lengthx{tempi,1} = temp_memoryPrecision_trialLevel_crossTime_highMatch;
        precision_trialLevel_crossTime_multiPeriod_4types_lengthx{tempi,2} = temp_memoryPrecision_trialLevel_crossTime_lowMatch;
        precision_trialLevel_crossTime_multiPeriod_4types_lengthx{tempi,3} = temp_memoryPrecision_trialLevel_crossTime_overMismatch;
        precision_trialLevel_crossTime_multiPeriod_4types_lengthx{tempi,4} = temp_memoryPrecision_trialLevel_crossTime_underMismatch;
        
        precision_trialLevel_crossTimeMean_multiPeriod_4types_lengthx{tempi,1} = mean(temp_memoryPrecision_trialLevel_crossTime_highMatch,1,'omitnan');
        precision_trialLevel_crossTimeMean_multiPeriod_4types_lengthx{tempi,2} = mean(temp_memoryPrecision_trialLevel_crossTime_lowMatch,1,'omitnan');
        precision_trialLevel_crossTimeMean_multiPeriod_4types_lengthx{tempi,3} = mean(temp_memoryPrecision_trialLevel_crossTime_overMismatch,1,'omitnan');
        precision_trialLevel_crossTimeMean_multiPeriod_4types_lengthx{tempi,4} = mean(temp_memoryPrecision_trialLevel_crossTime_underMismatch,1,'omitnan');
        
    end
    
    if tempj==1
        precision_trialLevel_crossTime_multiPeriod_4types_length1 = precision_trialLevel_crossTime_multiPeriod_4types_lengthx;
        precision_trialLevel_crossTimeMean_multiPeriod_4types_length1 = precision_trialLevel_crossTimeMean_multiPeriod_4types_lengthx;
    elseif tempj==2
        precision_trialLevel_crossTime_multiPeriod_4types_length2 = precision_trialLevel_crossTime_multiPeriod_4types_lengthx;
        precision_trialLevel_crossTimeMean_multiPeriod_4types_length2 = precision_trialLevel_crossTimeMean_multiPeriod_4types_lengthx;
    elseif tempj==3
        precision_trialLevel_crossTime_multiPeriod_4types_length3 = precision_trialLevel_crossTime_multiPeriod_4types_lengthx;
        precision_trialLevel_crossTimeMean_multiPeriod_4types_length3 = precision_trialLevel_crossTimeMean_multiPeriod_4types_lengthx;
    elseif tempj==4
        precision_trialLevel_crossTime_multiPeriod_4types_length123 = precision_trialLevel_crossTime_multiPeriod_4types_lengthx;
        precision_trialLevel_crossTimeMean_multiPeriod_4types_length123 = precision_trialLevel_crossTimeMean_multiPeriod_4types_lengthx;
    end
end

precision_trialLevel_crossTime_multiPeriod_4types_length1;
precision_trialLevel_crossTime_multiPeriod_4types_length2;
precision_trialLevel_crossTime_multiPeriod_4types_length3;
precision_trialLevel_crossTime_multiPeriod_4types_length123;

precision_trialLevel_crossTimeMean_multiPeriod_4types_length1;
precision_trialLevel_crossTimeMean_multiPeriod_4types_length2;
precision_trialLevel_crossTimeMean_multiPeriod_4types_length3;
precision_trialLevel_crossTimeMean_multiPeriod_4types_length123;

%% x_trialLevel_multiPeriod_4types_lengthx_collapsed
length1_crossTime_interval = [];
length2_crossTime_interval = [];
length3_crossTime_interval = [];

for tempi=1:4
    
    if tempi == 1
        temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_4types_length1;
        temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_4types_length1;
    elseif tempi == 2
        temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_4types_length2;
        temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_4types_length2;
    elseif tempi == 3
        temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_4types_length3;
        temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_4types_length3;
    elseif tempi == 4
        temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_4types_length123;
        temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_4types_length123;
    end
    
    temp_meta_trialLevel_size = [length(temp_meta_trialLevel{1,1}) length(temp_meta_trialLevel{2,1}) length(temp_meta_trialLevel{3,1})];
    temp_crossTime_interval = [1 temp_meta_trialLevel_size(1)+1 sum(temp_meta_trialLevel_size(1:2))+1 sum(temp_meta_trialLevel_size(1:3))];
    
    temp_meta_trialLevel_collapsed = nan(temp_crossTime_interval(end),4);
    temp_memoryPrecision_trialLevel_collapsed = nan(temp_crossTime_interval(end),4);
    for tempj=1:4
        temp1 = [];
        for tempk=1:3
            temp1 = [temp1 temp_meta_trialLevel{tempk,tempj}]; %#ok<*AGROW>
        end
        temp_meta_trialLevel_collapsed(:,tempj) = temp1;
        
        temp1 = [];
        for tempk=1:3
            temp1 = [temp1 temp_memoryPrecision_trialLevel{tempk,tempj}]; %#ok<*AGROW>
        end
        temp_memoryPrecision_trialLevel_collapsed(:,tempj) = temp1;
    end
    
    temp_memoryPrecision_trialLevel_collapsed_mean = mean(temp_memoryPrecision_trialLevel_collapsed,2);
    %temp_memoryPrecision_trialLevel_collapsed_mean = temp_memoryPrecision_trialLevel_collapsed(:,1);
    [M,I] = max(temp_memoryPrecision_trialLevel_collapsed_mean);
    
    
    if tempi == 1
        meta_trialLevel_multiPeriod_4types_length1_collapsed = temp_meta_trialLevel_collapsed;
        precision_trialLevel_multiPeriod_4types_length1_collapsed = temp_memoryPrecision_trialLevel_collapsed;
        length1_crossTime_interval = temp_crossTime_interval;
        length1_crossTime_max = I;
    elseif tempi == 2
        meta_trialLevel_multiPeriod_4types_length2_collapsed = temp_meta_trialLevel_collapsed;
        precision_trialLevel_multiPeriod_4types_length2_collapsed = temp_memoryPrecision_trialLevel_collapsed;
        length2_crossTime_interval = temp_crossTime_interval;
        length2_crossTime_max = I;
    elseif tempi == 3
        meta_trialLevel_multiPeriod_4types_length3_collapsed = temp_meta_trialLevel_collapsed;
        precision_trialLevel_multiPeriod_4types_length3_collapsed = temp_memoryPrecision_trialLevel_collapsed;
        length3_crossTime_interval = temp_crossTime_interval;
        length3_crossTime_max = I;
    elseif tempi == 4
        meta_trialLevel_multiPeriod_4types_length123_collapsed = temp_meta_trialLevel_collapsed;
        precision_trialLevel_multiPeriod_4types_length123_collapsed = temp_memoryPrecision_trialLevel_collapsed;
        length3_crossTime_interval = temp_crossTime_interval;
        length123_crossTime_max = I;
    end
end

%% meta_trialLevel_crossTime_multiPeriod_3types, each length
for tempj=1:4
    
    if tempj==1
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
    elseif tempj==2
        temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
    elseif tempj==3
        temp_trialBoolIndex_lengthx = trialBoolIndex_length3;
    elseif tempj==4
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
    end
    
    
    meta_trialLevel_crossTime_multiPeriod_3types_lengthx = cell(3,3);
    meta_trialLevel_crossTimeMean_multiPeriod_3types_lengthx = cell(3,3);
    for tempi=1:3
        if tempi==1
            temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_baseline;
        elseif tempi==2
            if tempj==1
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_length1;
            elseif tempj==2
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_length2;
            elseif tempj==3
                temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_length3;
            elseif tempj==4
                temp_meta_trialLevel_crossTime = nan(size(meta_trialLevel_crossTime_length3));
                
                for temptempi=1:size(temp_meta_trialLevel_crossTime,1)
                    if trialBoolIndex_length1(temptempi) == true
                        temp1 = meta_trialLevel_crossTime_length1(temptempi,:);
                        temp_interp_factor = 3;
                        
                        if isnan(temp1(1)) == true
                            continue
                        end
                        
                        x = 1:length(temp1);
                        x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);
                        
                        temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');
                        
                        for temptempj=1:ceil((temp_interp_factor-1))
                            temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
                        end
                        
                        temp_meta_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length2(temptempi) == true
                        temp1 = meta_trialLevel_crossTime_length2(temptempi,:);
                        temp_interp_factor = 1.5;
                        
                        if isnan(temp1(1)) == true
                            continue
                        end
                        
                        x = 1:length(temp1);
                        x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);
                        
                        temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');
                        
                        for temptempj=1:ceil((temp_interp_factor-1))
                            temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
                        end
                        
                        temp_meta_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length3(temptempi) == true
                        temp_meta_trialLevel_crossTime(temptempi,:) = meta_trialLevel_crossTime_length3(temptempi,:);
                    end
                    
                end
                temp_meta_trialLevel_crossTime;
                
            end
        elseif tempi==3
            temp_meta_trialLevel_crossTime = meta_trialLevel_crossTime_delay1;
        end
        
        temp_meta_trialLevel_crossTime_choiceMemoryCorrect = temp_meta_trialLevel_crossTime(choiceMemoryCorrectBoolIndex' & temp_trialBoolIndex_lengthx,:);
        temp_meta_trialLevel_crossTime_choiceOffload = temp_meta_trialLevel_crossTime(choiceOffloadBoolIndex' & temp_trialBoolIndex_lengthx,:);
        temp_meta_trialLevel_crossTime_choiceMemoryError = temp_meta_trialLevel_crossTime(choiceMemoryErrorBoolIndex' & temp_trialBoolIndex_lengthx,:);
        
        meta_trialLevel_crossTime_multiPeriod_3types_lengthx{tempi,1} = temp_meta_trialLevel_crossTime_choiceMemoryCorrect;
        meta_trialLevel_crossTime_multiPeriod_3types_lengthx{tempi,2} = temp_meta_trialLevel_crossTime_choiceOffload;
        meta_trialLevel_crossTime_multiPeriod_3types_lengthx{tempi,3} = temp_meta_trialLevel_crossTime_choiceMemoryError;
        
        meta_trialLevel_crossTimeMean_multiPeriod_3types_lengthx{tempi,1} = mean(temp_meta_trialLevel_crossTime_choiceMemoryCorrect,1,'omitnan');
        meta_trialLevel_crossTimeMean_multiPeriod_3types_lengthx{tempi,2} = mean(temp_meta_trialLevel_crossTime_choiceOffload,1,'omitnan');
        meta_trialLevel_crossTimeMean_multiPeriod_3types_lengthx{tempi,3} = mean(temp_meta_trialLevel_crossTime_choiceMemoryError,1,'omitnan');
        
    end
    
    if tempj==1
        meta_trialLevel_crossTime_multiPeriod_3types_length1 = meta_trialLevel_crossTime_multiPeriod_3types_lengthx;
        meta_trialLevel_crossTimeMean_multiPeriod_3types_length1 = meta_trialLevel_crossTimeMean_multiPeriod_3types_lengthx;
    elseif tempj==2
        meta_trialLevel_crossTime_multiPeriod_3types_length2 = meta_trialLevel_crossTime_multiPeriod_3types_lengthx;
        meta_trialLevel_crossTimeMean_multiPeriod_3types_length2 = meta_trialLevel_crossTimeMean_multiPeriod_3types_lengthx;
    elseif tempj==3
        meta_trialLevel_crossTime_multiPeriod_3types_length3 = meta_trialLevel_crossTime_multiPeriod_3types_lengthx;
        meta_trialLevel_crossTimeMean_multiPeriod_3types_length3 = meta_trialLevel_crossTimeMean_multiPeriod_3types_lengthx;
    elseif tempj==4
        meta_trialLevel_crossTime_multiPeriod_3types_length123 = meta_trialLevel_crossTime_multiPeriod_3types_lengthx;
        meta_trialLevel_crossTimeMean_multiPeriod_3types_length123 = meta_trialLevel_crossTimeMean_multiPeriod_3types_lengthx;
    end
end

meta_trialLevel_crossTime_multiPeriod_3types_length1;
meta_trialLevel_crossTime_multiPeriod_3types_length2;
meta_trialLevel_crossTime_multiPeriod_3types_length3;
meta_trialLevel_crossTime_multiPeriod_3types_length123;

meta_trialLevel_crossTimeMean_multiPeriod_3types_length1;
meta_trialLevel_crossTimeMean_multiPeriod_3types_length2;
meta_trialLevel_crossTimeMean_multiPeriod_3types_length3;
meta_trialLevel_crossTimeMean_multiPeriod_3types_length123;


%% precision_trialLevel_crossTime_multiPeriod_3types, each length
for tempj=1:4
    
    if tempj==1
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1;
    elseif tempj==2
        temp_trialBoolIndex_lengthx = trialBoolIndex_length2;
    elseif tempj==3
        temp_trialBoolIndex_lengthx = trialBoolIndex_length3;
    elseif tempj==4
        temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
    end
    
    precision_trialLevel_crossTime_multiPeriod_3types_lengthx = cell(3,3);
    precision_trialLevel_crossTimeMean_multiPeriod_3types_lengthx = cell(3,3);
    for tempi=1:3
        if tempi==1
            temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_baseline;
        elseif tempi==2
            if tempj==1
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_length1;
            elseif tempj==2
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_length2;
            elseif tempj==3
                temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_length3;
            elseif tempj==4
                temp_memoryPrecision_trialLevel_crossTime = nan(size(memoryPrecision_trialLevel_crossTime_length3));
                
                for temptempi=1:size(temp_memoryPrecision_trialLevel_crossTime,1)
                    if trialBoolIndex_length1(temptempi) == true
                        temp1 = memoryPrecision_trialLevel_crossTime_length1(temptempi,:);
                        temp_interp_factor = 3;
                        
                        if isnan(temp1(1)) == true
                            continue
                        end
                        
                        x = 1:length(temp1);
                        x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);
                        
                        temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');
                        
                        for temptempj=1:ceil((temp_interp_factor-1))
                            temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
                        end
                        
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length2(temptempi) == true
                        temp1 = memoryPrecision_trialLevel_crossTime_length2(temptempi,:);
                        temp_interp_factor = 1.5;
                        
                        if isnan(temp1(1)) == true
                            continue
                        end
                        
                        x = 1:length(temp1);
                        x_interp = (1:temp_interp_factor:length(temp1)*temp_interp_factor)+(temp_interp_factor-1);
                        
                        temp1_interp = interp1(x_interp,temp1,1:length(temp1)*temp_interp_factor,'makima');
                        
                        for temptempj=1:ceil((temp_interp_factor-1))
                            temp1_interp(temptempj) = temp1_interp(ceil(temp_interp_factor));
                        end
                        
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1_interp;
                    end
                    if trialBoolIndex_length3(temptempi) == true
                        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = memoryPrecision_trialLevel_crossTime_length3(temptempi,:);
                    end
                    
                end
                temp_memoryPrecision_trialLevel_crossTime;
                
            end
        elseif tempi==3
            temp_memoryPrecision_trialLevel_crossTime = memoryPrecision_trialLevel_crossTime_delay1;
        end
        
        
        temp_memoryPrecision_trialLevel_crossTime_choiceMemoryCorrect = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryCorrectBoolIndex' & temp_trialBoolIndex_lengthx,:);
        temp_memoryPrecision_trialLevel_crossTime_choiceOffload = temp_memoryPrecision_trialLevel_crossTime(choiceOffloadBoolIndex' & temp_trialBoolIndex_lengthx,:);
        temp_memoryPrecision_trialLevel_crossTime_choiceMemoryError = temp_memoryPrecision_trialLevel_crossTime(choiceMemoryErrorBoolIndex' & temp_trialBoolIndex_lengthx,:);
        
        precision_trialLevel_crossTime_multiPeriod_3types_lengthx{tempi,1} = temp_memoryPrecision_trialLevel_crossTime_choiceMemoryCorrect;
        precision_trialLevel_crossTime_multiPeriod_3types_lengthx{tempi,2} = temp_memoryPrecision_trialLevel_crossTime_choiceOffload;
        precision_trialLevel_crossTime_multiPeriod_3types_lengthx{tempi,3} = temp_memoryPrecision_trialLevel_crossTime_choiceMemoryError;
        
        precision_trialLevel_crossTimeMean_multiPeriod_3types_lengthx{tempi,1} = mean(temp_memoryPrecision_trialLevel_crossTime_choiceMemoryCorrect,1,'omitnan');
        precision_trialLevel_crossTimeMean_multiPeriod_3types_lengthx{tempi,2} = mean(temp_memoryPrecision_trialLevel_crossTime_choiceOffload,1,'omitnan');
        precision_trialLevel_crossTimeMean_multiPeriod_3types_lengthx{tempi,3} = mean(temp_memoryPrecision_trialLevel_crossTime_choiceMemoryError,1,'omitnan');
        
    end
    
    if tempj==1
        precision_trialLevel_crossTime_multiPeriod_3types_length1 = precision_trialLevel_crossTime_multiPeriod_3types_lengthx;
        precision_trialLevel_crossTimeMean_multiPeriod_3types_length1 = precision_trialLevel_crossTimeMean_multiPeriod_3types_lengthx;
    elseif tempj==2
        precision_trialLevel_crossTime_multiPeriod_3types_length2 = precision_trialLevel_crossTime_multiPeriod_3types_lengthx;
        precision_trialLevel_crossTimeMean_multiPeriod_3types_length2 = precision_trialLevel_crossTimeMean_multiPeriod_3types_lengthx;
    elseif tempj==3
        precision_trialLevel_crossTime_multiPeriod_3types_length3 = precision_trialLevel_crossTime_multiPeriod_3types_lengthx;
        precision_trialLevel_crossTimeMean_multiPeriod_3types_length3 = precision_trialLevel_crossTimeMean_multiPeriod_3types_lengthx;
    elseif tempj==4
        precision_trialLevel_crossTime_multiPeriod_3types_length123 = precision_trialLevel_crossTime_multiPeriod_3types_lengthx;
        precision_trialLevel_crossTimeMean_multiPeriod_3types_length123 = precision_trialLevel_crossTimeMean_multiPeriod_3types_lengthx;
    end
end

precision_trialLevel_crossTime_multiPeriod_3types_length1;
precision_trialLevel_crossTime_multiPeriod_3types_length2;
precision_trialLevel_crossTime_multiPeriod_3types_length3;
precision_trialLevel_crossTime_multiPeriod_3types_length123;

precision_trialLevel_crossTimeMean_multiPeriod_3types_length1;
precision_trialLevel_crossTimeMean_multiPeriod_3types_length2;
precision_trialLevel_crossTimeMean_multiPeriod_3types_length3;
precision_trialLevel_crossTimeMean_multiPeriod_3types_length123;


%% x_trialLevel_multiPeriod_3types_lengthx_collapsed
length1_crossTime_interval_3types = [];
length2_crossTime_interval_3types = [];
length3_crossTime_interval_3types = [];

for tempi=1:4
    
    if tempi == 1
        temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_3types_length1;
        temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_3types_length1;
        temp_meta_trialLevel_eachTrial = meta_trialLevel_crossTime_multiPeriod_3types_length1;
        temp_memoryPrecision_trialLevel_eachTrial = precision_trialLevel_crossTime_multiPeriod_3types_length1;
    elseif tempi == 2
        temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_3types_length2;
        temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_3types_length2;
        temp_meta_trialLevel_eachTrial = meta_trialLevel_crossTime_multiPeriod_3types_length2;
        temp_memoryPrecision_trialLevel_eachTrial = precision_trialLevel_crossTime_multiPeriod_3types_length2;
    elseif tempi == 3
        temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_3types_length3;
        temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_3types_length3;
        temp_meta_trialLevel_eachTrial = meta_trialLevel_crossTime_multiPeriod_3types_length3;
        temp_memoryPrecision_trialLevel_eachTrial = precision_trialLevel_crossTime_multiPeriod_3types_length3;
    elseif tempi == 4
        temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_3types_length123;
        temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_3types_length123;
        temp_meta_trialLevel_eachTrial = meta_trialLevel_crossTime_multiPeriod_3types_length123;
        temp_memoryPrecision_trialLevel_eachTrial = precision_trialLevel_crossTime_multiPeriod_3types_length123;
    end
    
    temp_meta_trialLevel_size = [length(temp_meta_trialLevel{1,1}) length(temp_meta_trialLevel{2,1}) length(temp_meta_trialLevel{3,1})];
    temp_crossTime_interval = [1 temp_meta_trialLevel_size(1)+1 sum(temp_meta_trialLevel_size(1:2))+1 sum(temp_meta_trialLevel_size(1:3))];
    
    % temp_x_trialLevel_collapsed
    temp_meta_trialLevel_collapsed = nan(temp_crossTime_interval(end),3);
    temp_memoryPrecision_trialLevel_collapsed = nan(temp_crossTime_interval(end),3);
    for tempj=1:3
        temp1 = [];
        for tempk=1:3
            temp1 = [temp1 temp_meta_trialLevel{tempk,tempj}]; %#ok<*AGROW>
        end
        temp_meta_trialLevel_collapsed(:,tempj) = temp1;
        
        temp1 = [];
        for tempk=1:3
            temp1 = [temp1 temp_memoryPrecision_trialLevel{tempk,tempj}]; %#ok<*AGROW>
        end
        temp_memoryPrecision_trialLevel_collapsed(:,tempj) = temp1;
    end
    
    temp_memoryPrecision_trialLevel_collapsed_mean = mean(temp_memoryPrecision_trialLevel_collapsed,2);
    [M,I] = max(temp_memoryPrecision_trialLevel_collapsed_mean);
    
    % temp_meta_trialLevel_eachTrial_collapsed
    temp_meta_trialLevel_eachTrial_collapsed = cell(1,3);
    temp_memoryPrecision_trialLevel_eachTrial_collapsed = cell(1,3);
    for tempj=1:3
        temp_meta_trialLevel_eachTrial_collapsed{tempj} = nan(temp_crossTime_interval(end),size(temp_meta_trialLevel_eachTrial{1,tempj},1));
        
        temp1 = [];
        for tempk=1:3
            temp1 = [temp1 temp_meta_trialLevel_eachTrial{tempk,tempj}]; %#ok<*AGROW>
        end
        temp_meta_trialLevel_eachTrial_collapsed{tempj} = temp1';
        
        temp1 = [];
        for tempk=1:3
            temp1 = [temp1 temp_memoryPrecision_trialLevel_eachTrial{tempk,tempj}]; %#ok<*AGROW>
        end
        temp_memoryPrecision_trialLevel_eachTrial_collapsed{tempj} = temp1';
    end
    
    
    % end
    if tempi == 1
        meta_trialLevel_multiPeriod_3types_length1_collapsed = temp_meta_trialLevel_collapsed;
        precision_trialLevel_multiPeriod_3types_length1_collapsed = temp_memoryPrecision_trialLevel_collapsed;
        meta_eachTrial_multiPeriod_3types_length1_collapsed = temp_meta_trialLevel_eachTrial_collapsed;
        precision_eachTrial_multiPeriod_3types_length1_collapsed = temp_memoryPrecision_trialLevel_eachTrial_collapsed;
        length1_crossTime_interval_3types = temp_crossTime_interval;
        length1_crossTime_max_3types = I;
    elseif tempi == 2
        meta_trialLevel_multiPeriod_3types_length2_collapsed = temp_meta_trialLevel_collapsed;
        precision_trialLevel_multiPeriod_3types_length2_collapsed = temp_memoryPrecision_trialLevel_collapsed;
        meta_eachTrial_multiPeriod_3types_length2_collapsed = temp_meta_trialLevel_eachTrial_collapsed;
        precision_eachTrial_multiPeriod_3types_length2_collapsed = temp_memoryPrecision_trialLevel_eachTrial_collapsed;
        length2_crossTime_interval_3types = temp_crossTime_interval;
        length2_crossTime_max_3types = I;
    elseif tempi == 3
        meta_trialLevel_multiPeriod_3types_length3_collapsed = temp_meta_trialLevel_collapsed;
        precision_trialLevel_multiPeriod_3types_length3_collapsed = temp_memoryPrecision_trialLevel_collapsed;
        meta_eachTrial_multiPeriod_3types_length3_collapsed = temp_meta_trialLevel_eachTrial_collapsed;
        precision_eachTrial_multiPeriod_3types_length3_collapsed = temp_memoryPrecision_trialLevel_eachTrial_collapsed;
        length3_crossTime_interval_3types = temp_crossTime_interval;
        length3_crossTime_max_3types = I;
    elseif tempi == 4
        meta_trialLevel_multiPeriod_3types_length123_collapsed = temp_meta_trialLevel_collapsed;
        precision_trialLevel_multiPeriod_3types_length123_collapsed = temp_memoryPrecision_trialLevel_collapsed;
        meta_eachTrial_multiPeriod_3types_length123_collapsed = temp_meta_trialLevel_eachTrial_collapsed;
        precision_eachTrial_multiPeriod_3types_length123_collapsed = temp_memoryPrecision_trialLevel_eachTrial_collapsed;
        length3_crossTime_interval_3types = temp_crossTime_interval;
        length123_crossTime_max_3types = I;
    end
end


%% Get new boundary
meta_trialLevel_crossTime_multiPeriod_4types_length123;
% precision_trialLevel_crossTime_multiPeriod_4types_length123;

length3_crossTime_interval;
length123_crossTime_max;


temp_time_max = length123_crossTime_max-length3_crossTime_interval(3)+1;

temp1 = meta_trialLevel_crossTime_multiPeriod_4types_length123{3,1}(:,temp_time_max);
temp2 = meta_trialLevel_crossTime_multiPeriod_4types_length123{3,2}(:,temp_time_max);
temp3 = meta_trialLevel_crossTime_multiPeriod_4types_length123{3,3}(:,temp_time_max);
temp4 = meta_trialLevel_crossTime_multiPeriod_4types_length123{3,4}(:,temp_time_max);
temp1234 = [temp1;temp2;temp3;temp4];
temp_meta_median_length123 = median(temp1234,'omitnan');
temp_meta_mean_length123 = mean(temp1234,'omitnan');


temp1 = precision_trialLevel_crossTime_multiPeriod_4types_length123{3,1}(:,temp_time_max);
temp2 = precision_trialLevel_crossTime_multiPeriod_4types_length123{3,2}(:,temp_time_max);
temp3 = precision_trialLevel_crossTime_multiPeriod_4types_length123{3,3}(:,temp_time_max);
temp4 = precision_trialLevel_crossTime_multiPeriod_4types_length123{3,4}(:,temp_time_max);
temp1234 = [temp1;temp2;temp3;temp4];
temp_precision_median_length123 = median(temp1234,'omitnan');
temp_precision_mean_length123 = mean(temp1234,'omitnan');

%%
temp_markerList = ['o','^','s','x'];

%% Plot for 4types
if false
    fig = figure('Name','asd','NumberTitle','off'); %#ok<*NASGU>
    %set(gcf,'Position',[510 90 720 480*1.2]);
    %t = tiledlayout(3,4,'TileSpacing','Compact','Padding','Compact');
    set(gcf,'Position',[510 90 720 480*1.2*1.2*0.93]);
    t = tiledlayout(4,4,'TileSpacing','Compact','Padding','Compact');
    
    x = [precision_trialLevel_multiPeriod_4types_length1_collapsed(1:length1_crossTime_max,:);...
        precision_trialLevel_multiPeriod_4types_length2_collapsed(1:length2_crossTime_max,:);...
        precision_trialLevel_multiPeriod_4types_length3_collapsed(1:length3_crossTime_max,:);...
        precision_trialLevel_multiPeriod_4types_length123_collapsed(1:length123_crossTime_max,:)];
    y = [meta_trialLevel_multiPeriod_4types_length1_collapsed(1:length1_crossTime_max,:);...
        meta_trialLevel_multiPeriod_4types_length2_collapsed(1:length2_crossTime_max,:);...
        meta_trialLevel_multiPeriod_4types_length3_collapsed(1:length3_crossTime_max,:);...
        meta_trialLevel_multiPeriod_4types_length123_collapsed(1:length123_crossTime_max,:)];
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:));
    
    
    % temp_markerList = ['o','x','s','|'];
    % temp_markerList = ['o','x','s','+'];
    %temp_markerList = ['o','^','s','x'];
    
    for tempi=1:4
        if tempi == 1
            temp_strA = 'Len1';
            
            temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_4types_length1_collapsed;
            temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_4types_length1_collapsed;
            temp_crossTime_interval = length1_crossTime_interval;
            temp_crossTime_max = length1_crossTime_max;
            
        elseif tempi == 2
            temp_strA = 'Len2';
            
            temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_4types_length2_collapsed;
            temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_4types_length2_collapsed;
            temp_crossTime_interval = length2_crossTime_interval;
            temp_crossTime_max = length2_crossTime_max;
            
        elseif tempi == 3
            temp_strA = 'Len3';
            
            temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_4types_length3_collapsed;
            temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_4types_length3_collapsed;
            temp_crossTime_interval = length3_crossTime_interval;
            temp_crossTime_max = length3_crossTime_max;
        elseif tempi == 4
            temp_strA = 'Len123';
            
            temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_4types_length123_collapsed;
            temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_4types_length123_collapsed;
            temp_crossTime_interval = length3_crossTime_interval;
            temp_crossTime_max = length123_crossTime_max;
        end
        
        for tempPeriodIndex=1:4
            nexttile
            
            if tempPeriodIndex == 1
                temp_strB = 'Baseline';
                x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):(temp_crossTime_interval(2)),:);
                y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):(temp_crossTime_interval(2)),:);
                temp_crossTime_marker = [temp_crossTime_interval(1) (temp_crossTime_interval(2))];
                temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
                %temptemp_markerList = temp_markerList([tempPeriodIndex end]);
                temptemp_markerList = temp_markerList([tempPeriodIndex tempPeriodIndex+1]);
            elseif tempPeriodIndex == 2
                temp_strB = 'Sample';
                x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(2):(temp_crossTime_interval(3)),:);
                y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(2):(temp_crossTime_interval(3)),:);
                temp_crossTime_marker = [temp_crossTime_interval(2) (temp_crossTime_interval(3))];
                temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
                %temptemp_markerList = temp_markerList([tempPeriodIndex end]);
                temptemp_markerList = temp_markerList([tempPeriodIndex tempPeriodIndex+1]);
            elseif tempPeriodIndex == 3
                temp_strB = 'Delay1';
                %x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(3):(temp_crossTime_interval(4)),:);
                %y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(3):(temp_crossTime_interval(4)),:);
                %temp_crossTime_marker = [temp_crossTime_interval(3) temp_crossTime_interval(4)];
                
                x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(3):temp_crossTime_max,:);
                y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(3):temp_crossTime_max,:);
                temp_crossTime_marker = [temp_crossTime_interval(3) temp_crossTime_max];
                
                temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
                %temptemp_markerList = temp_markerList([tempPeriodIndex end]);
                temptemp_markerList = temp_markerList([tempPeriodIndex tempPeriodIndex+1]);
            elseif tempPeriodIndex == 4
                temp_strB = 'All';
                %x = temp_memoryPrecision_trialLevel_collapsed;
                %y = temp_meta_trialLevel_collapsed;
                %temp_crossTime_marker = temp_crossTime_interval;
                
                x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
                y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
                temp_crossTime_marker = temp_crossTime_interval;
                temp_crossTime_marker(end) = temp_crossTime_max;
                
                temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
                temptemp_markerList = temp_markerList;
            end
            
            
            %         plot(x(:,1), y(:,1), '-', 'LineWidth', 0.5, 'Color', color_choiceMemoryHigh);
            %         hold on
            %
            %         plot(x(:,2), y(:,2), '-', 'LineWidth', 0.5, 'Color', color_choiceOffloadLow);
            %         hold on
            %
            %         plot(x(:,3), y(:,3), '-', 'LineWidth', 0.5, 'Color', color_choiceMemoryLow);
            %         hold on
            %
            %         plot(x(:,4), y(:,4), '-', 'LineWidth', 0.5, 'Color', color_choiceOffloadHigh);
            %         hold on
            
            
            plot(x(:,1), y(:,1), ':', 'LineWidth', 0.75, 'Color', color_choiceMemoryHigh);
            hold on
            
            plot(x(:,2), y(:,2), ':', 'LineWidth', 0.75, 'Color', color_choiceOffloadLow);
            hold on
            
            plot(x(:,3), y(:,3), ':', 'LineWidth', 0.75, 'Color', color_choiceMemoryLow);
            hold on
            
            plot(x(:,4), y(:,4), ':', 'LineWidth', 0.75, 'Color', color_choiceOffloadHigh);
            hold on
            
            
            %         plot(x(temp_crossTime_marker,1), y(temp_crossTime_marker,1), '-', 'LineWidth', 0.5, 'Color', color_choiceMemoryHigh);
            %         hold on
            %
            %         plot(x(temp_crossTime_marker,2), y(temp_crossTime_marker,2), '-', 'LineWidth', 0.5, 'Color', color_choiceOffloadLow);
            %         hold on
            %
            %         plot(x(temp_crossTime_marker,3), y(temp_crossTime_marker,3), '-', 'LineWidth', 0.5, 'Color', color_choiceMemoryLow);
            %         hold on
            %
            %         plot(x(temp_crossTime_marker,4), y(temp_crossTime_marker,4), '-', 'LineWidth', 0.5, 'Color', color_choiceOffloadHigh);
            %         hold on
            
            
            for temptempi=1:length(temp_crossTime_marker)
                temp_marker = temptemp_markerList(temptempi);
                
                temptemp_markerTime = temp_crossTime_marker(temptempi);
                
                temp_size = 13;%20
                if strcmp(temp_marker,'x')
                    temp_size = 50;%20
                end
                if strcmp(temp_marker,'+')
                    temp_size = 20;
                end
                temp_alpha = 1;
                temp_LineWidth = 1.5;
                
                scatter(x(temptemp_markerTime,1),y(temptemp_markerTime,1),...
                    temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryHigh,'LineWidth',temp_LineWidth,...
                    'MarkerEdgeAlpha',temp_alpha);
                hold on
                scatter(x(temptemp_markerTime,2),y(temptemp_markerTime,2),...
                    temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffloadLow,'LineWidth',temp_LineWidth,...
                    'MarkerEdgeAlpha',temp_alpha);
                hold on
                scatter(x(temptemp_markerTime,3),y(temptemp_markerTime,3),...
                    temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryLow,'LineWidth',temp_LineWidth,...
                    'MarkerEdgeAlpha',temp_alpha);
                hold on
                scatter(x(temptemp_markerTime,4),y(temptemp_markerTime,4),...
                    temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffloadHigh,'LineWidth',temp_LineWidth,...
                    'MarkerEdgeAlpha',temp_alpha);
                hold on
                
            end
            
            
            xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
            ylim([temp_ymin-(temp_ymax-temp_ymin)*0.08 temp_ymax+(temp_ymax-temp_ymin)*0.08]);
            %ylim([0 1]);
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 10)
            set(gca,'box','off');% 取消右、上边框
            xlabel(sprintf('Memory precision'), 'FontSize', 10, 'FontWeight', 'normal');
            ylabel(sprintf('Meta-memory'), 'FontSize', 10, 'FontWeight', 'normal');
            
            temp_title = title(sprintf('%s, %s',temp_strA,temp_strB), 'FontSize', 10, 'FontWeight', 'normal');
            temp_title.Interpreter = 'none';
            
        end
    end
    
end


% color_choiceMemoryHigh_B = [128,205,193]/255; %[146,197,222]/255
% color_choiceMemoryLow_B = [1,133,113]/255; %[5,113,176]/255
% color_choiceOffloadHigh_B = [166,97,26]/255; %[202,0,32]/255
% color_choiceOffloadLow_B = [223,194,125]/255; %[244,165,130]/255

color_choiceMemoryHigh_B = [217,240,236]/255;
color_choiceMemoryLow_B = [178,218,212]/255;
color_choiceOffloadHigh_B = [228,207,186]/255;
color_choiceOffloadLow_B = [245,237,216]/255;

%% Plot for 4types (3d), x is time
if true
    fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
    %set(gcf,'Position',[10 150 450*0.8*0.8*0.8*0.95 450*0.8*0.8*0.8*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    set(gcf,'Position',[10 150 450*0.8*0.8*0.8*0.95*1.63*1.08*1.04 450*0.8*0.8*0.8*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    t = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');
    
    temp_strA = 'Len123';
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_4types_length123_collapsed;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_4types_length123_collapsed;
    temp_crossTime_interval = length3_crossTime_interval;
    temp_crossTime_max = length123_crossTime_max;
    
    
    temp_strB = 'All';
    y = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    z = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    temptemp_markerList = temp_markerList;
    
    x = 1:temp_crossTime_marker(end);
    
    windowSize = 8;%8
    y = smoothdata(y,1,'gaussian',windowSize);
    z = smoothdata(z,1,'gaussian',windowSize);
    
    temp_LineWidth = 0.75;
    temp_Linestyle = '-';
    
    plot3(x, y(:,1),z(:,1), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryHigh);
    hold on
    
    plot3(x, y(:,2),z(:,2), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceOffloadLow);
    hold on
    
    plot3(x, y(:,3),z(:,3), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryLow);
    hold on
    
    plot3(x, y(:,4),z(:,4), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceOffloadHigh);
    hold on
    
    
    
    XYBinLimits = [0 1];
    
    temp_offset = 0.003;%0.02
    temp_offset2 = 0.015;%0.015
    temp_lineWidth = 1;%4,1.5
    
    temp_threshold_memoryPrecision = temp_precision_mean_length123;
    temp_threshold_meta = temp_meta_mean_length123;
    
    %     if if_4typesFilter_all0_cMC1_cF2_cME3 == 3
    %         temp_threshold_memoryPrecision = lowThreshold_memoryPrecision;
    %         temp_threshold_meta = lowThreshold_meta;
    %     end
    
    temp_facealpha = 0.9;%0.3
    
    %         % high-match
    %         temp2 = [temp_threshold_memoryPrecision temp_threshold_memoryPrecision 1 1];
    %         temp3 = [temp_threshold_meta 1 1 temp_threshold_meta];
    %         temp1 = [1 1 1 1].*temp_crossTime_marker(2);
    %         h = fill3(temp1,temp2,temp3,color_choiceMemoryHigh_B);
    %         set(h,'edgealpha',0,'facealpha',temp_facealpha)
    %
    %         % under-mismatch
    %         temp2 = [temp_threshold_memoryPrecision temp_threshold_memoryPrecision 1 1];
    %         temp3 = [0 temp_threshold_meta temp_threshold_meta 0];
    %         temp1 = [1 1 1 1].*temp_crossTime_marker(2);
    %         h = fill3(temp1,temp2,temp3,color_choiceOffloadHigh_B);
    %         set(h,'edgealpha',0,'facealpha',temp_facealpha)
    %
    %         % low-match
    %         temp2 = [0 0 temp_threshold_memoryPrecision temp_threshold_memoryPrecision];
    %         temp3 = [0 temp_threshold_meta temp_threshold_meta 0];
    %         temp1 = [1 1 1 1].*temp_crossTime_marker(2);
    %         h = fill3(temp1,temp2,temp3,color_choiceOffloadLow_B);
    %         set(h,'edgealpha',0,'facealpha',temp_facealpha)
    %
    %         % over-mismatch
    %         temp2 = [0 0 temp_threshold_memoryPrecision temp_threshold_memoryPrecision];
    %         temp3 = [temp_threshold_meta 1 1 temp_threshold_meta];
    %         temp1 = [1 1 1 1].*temp_crossTime_marker(2);
    %         h = fill3(temp1,temp2,temp3,color_choiceMemoryLow_B);
    %         set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    
    % high-match
    temp2 = [temp_threshold_memoryPrecision temp_threshold_memoryPrecision 1 1];
    temp3 = [temp_threshold_meta 1 1 temp_threshold_meta];
    temp1 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceMemoryHigh_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % under-mismatch
    temp2 = [temp_threshold_memoryPrecision temp_threshold_memoryPrecision 1 1];
    temp3 = [0 temp_threshold_meta temp_threshold_meta 0];
    temp1 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceOffloadHigh_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % low-match
    temp2 = [0 0 temp_threshold_memoryPrecision temp_threshold_memoryPrecision];
    temp3 = [0 temp_threshold_meta temp_threshold_meta 0];
    temp1 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceOffloadLow_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % over-mismatch
    temp2 = [0 0 temp_threshold_memoryPrecision temp_threshold_memoryPrecision];
    temp3 = [temp_threshold_meta 1 1 temp_threshold_meta];
    temp1 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceMemoryLow_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    
    for temptempi=1:length(temp_crossTime_marker)
        
        %if temptempi==2 || temptempi==3
        if temptempi==3
            continue
        end
        
        temp_marker = temptemp_markerList(temptempi);
        
        temptemp_markerTime = temp_crossTime_marker(temptempi);
        
        temp_size = 13*1.5*1.3;%20
        if strcmp(temp_marker,'x')
            temp_size = 50*1.5;%20
        end
        temp_alpha = 1;
        temp_LineWidth = 1.75;%1.5
        
        scatter3(temptemp_markerTime,y(temptemp_markerTime,1),z(temptemp_markerTime,1),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryHigh,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temptemp_markerTime,y(temptemp_markerTime,2),z(temptemp_markerTime,2),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffloadLow,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temptemp_markerTime,y(temptemp_markerTime,3),z(temptemp_markerTime,3),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryLow,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temptemp_markerTime,y(temptemp_markerTime,4),z(temptemp_markerTime,4),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffloadHigh,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        
    end
    
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:));
    [temp_zmin,temp_zmax] = bounds(z(:));
    
    %tempSize = 3;%10
    %scatter3(x,y,z,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    %hold on
    
    %axis equal
    
    %xticks([]);
    %yticks([]);
    
    %temp_xmin = 0;
    
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
    xticklabels({'Fixation','Sample','Delay1'});
    
    %grid on
    
    set(gca,'ydir','reverse');
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.02 temp_xmax+(temp_xmax-temp_xmin)*0.02]);%0.08
    %ylim([temp_ymin-(temp_ymax-temp_ymin)*0.12 temp_ymax+(temp_ymax-temp_ymin)*0.12]);
    ylim([temp_ymin temp_ymax+(temp_ymax-temp_ymin)*0.12]);
    zlim([temp_zmin-(temp_zmax-temp_zmin)*0.12 temp_zmax+(temp_zmax-temp_zmin)*0.12]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 7.5)%9
    set(gca,'box','off');% 取消右、上边框
    %xlabel('Memory precision', 'FontSize', 9, 'FontWeight', 'normal');
    %ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
    xlabel('Time', 'FontSize', 9, 'FontWeight', 'normal');
    
end

%% Plot for 4types (3d), z is time
if false
% if true
    fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
    set(gcf,'Position',[10 150 450*0.8*0.8*0.8 450*0.8*0.8*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    t = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');
    
    temp_strA = 'Len123';
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_4types_length123_collapsed;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_4types_length123_collapsed;
    temp_crossTime_interval = length3_crossTime_interval;
    %temp_crossTime_max = length123_crossTime_max;
    temp_crossTime_max = 121;%89,121
    
    temp_strB = 'All';
    x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    temptemp_markerList = temp_markerList;
    
    z = 1:temp_crossTime_marker(end);
    
    windowSize = 8;%8
    %x_raw = x;
    x = smoothdata(x,1,'gaussian',windowSize);
    y = smoothdata(y,1,'gaussian',windowSize);
    
    %temp_LineWidth = 1;%1
    %temp_Linestyle = ':';%':'
    
    temp_LineWidth = 0.75;
    temp_Linestyle = '-';
    
    plot3(x(:,1), y(:,1),z, temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryHigh);
    hold on
    
    plot3(x(:,2), y(:,2),z, temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceOffloadLow);
    hold on
    
    plot3(x(:,3), y(:,3),z, temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryLow);
    hold on
    
    plot3(x(:,4), y(:,4),z, temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceOffloadHigh);
    hold on
    
    
    
    XYBinLimits = [0 1];
    
    temp_offset = 0.003;%0.02
    temp_offset2 = 0.015;%0.015
    temp_lineWidth = 1;%4,1.5
    
    %temp_threshold_memoryPrecision = lowThreshold_memoryPrecision;
    %temp_threshold_meta = lowThreshold_meta;
    %temp_threshold_memoryPrecision = temp_precision_median_length123;
    %temp_threshold_meta = temp_meta_median_length123;
    temp_threshold_memoryPrecision = temp_precision_mean_length123;
    temp_threshold_meta = temp_meta_mean_length123;
    
    temp_facealpha = 0.9;%0.3
    
    % high-match
    temp1 = [temp_threshold_memoryPrecision temp_threshold_memoryPrecision 1 1];
    temp2 = [temp_threshold_meta 1 1 temp_threshold_meta];
    temp3 = [1 1 1 1].*temp_crossTime_marker(2);
    h = fill3(temp1,temp2,temp3,color_choiceMemoryHigh_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % under-mismatch
    temp1 = [temp_threshold_memoryPrecision temp_threshold_memoryPrecision 1 1];
    temp2 = [0 temp_threshold_meta temp_threshold_meta 0];
    temp3 = [1 1 1 1].*temp_crossTime_marker(2);
    h = fill3(temp1,temp2,temp3,color_choiceOffloadHigh_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % low-match
    temp1 = [0 0 temp_threshold_memoryPrecision temp_threshold_memoryPrecision];
    temp2 = [0 temp_threshold_meta temp_threshold_meta 0];
    temp3 = [1 1 1 1].*temp_crossTime_marker(2);
    h = fill3(temp1,temp2,temp3,color_choiceOffloadLow_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % over-mismatch
    temp1 = [0 0 temp_threshold_memoryPrecision temp_threshold_memoryPrecision];
    temp2 = [temp_threshold_meta 1 1 temp_threshold_meta];
    temp3 = [1 1 1 1].*temp_crossTime_marker(2);
    h = fill3(temp1,temp2,temp3,color_choiceMemoryLow_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    
    % high-match
    temp1 = [temp_threshold_memoryPrecision temp_threshold_memoryPrecision 1 1];
    temp2 = [temp_threshold_meta 1 1 temp_threshold_meta];
    temp3 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceMemoryHigh_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % under-mismatch
    temp1 = [temp_threshold_memoryPrecision temp_threshold_memoryPrecision 1 1];
    temp2 = [0 temp_threshold_meta temp_threshold_meta 0];
    temp3 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceOffloadHigh_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % low-match
    temp1 = [0 0 temp_threshold_memoryPrecision temp_threshold_memoryPrecision];
    temp2 = [0 temp_threshold_meta temp_threshold_meta 0];
    temp3 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceOffloadLow_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % over-mismatch
    temp1 = [0 0 temp_threshold_memoryPrecision temp_threshold_memoryPrecision];
    temp2 = [temp_threshold_meta 1 1 temp_threshold_meta];
    temp3 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceMemoryLow_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    
    %             plot3([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
    %                 [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
    %             hold on
    %             plot3([1 1].*(XYBinLimits(1)+temp_offset2),...
    %                 [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
    %             hold on
    %             plot3([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
    %                 [1 1].*(temp_threshold_meta+temp_offset),...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
    %             hold on
    %             plot3([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
    %                 [1 1].*(XYBinLimits(2)-temp_offset2),...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceMemoryLow);
    %             hold on
    %
    %
    %             plot3([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
    %                 [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
    %             hold on
    %             plot3([1 1].*(XYBinLimits(2)-temp_offset2),...
    %                 [temp_threshold_meta+temp_offset XYBinLimits(2)-temp_offset2],...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
    %             hold on
    %             plot3([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
    %                 [1 1].*(temp_threshold_meta+temp_offset),...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
    %             hold on
    %             plot3([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
    %                 [1 1].*(XYBinLimits(2)-temp_offset2),...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceMemoryHigh);
    %             hold on
    %
    %
    %
    %
    %             plot3([1 1].*(temp_threshold_memoryPrecision-temp_offset),...
    %                 [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
    %             hold on
    %             plot3([1 1].*(XYBinLimits(1)+temp_offset2),...
    %                 [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
    %             hold on
    %             plot3([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
    %                 [1 1].*(temp_threshold_meta-temp_offset),...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
    %             hold on
    %             plot3([temp_offset2+XYBinLimits(1) temp_threshold_memoryPrecision-temp_offset],...
    %                 [1 1].*(XYBinLimits(1)+temp_offset2),...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceOffloadLow);
    %             hold on
    %
    %
    %             plot3([1 1].*(temp_threshold_memoryPrecision+temp_offset),...
    %                 [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
    %             hold on
    %             plot3([1 1].*(XYBinLimits(2)-temp_offset2),...
    %                 [XYBinLimits(1)+temp_offset2 temp_threshold_meta-temp_offset],...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
    %             hold on
    %             plot3([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
    %                 [1 1].*(temp_threshold_meta-temp_offset),...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
    %             hold on
    %             plot3([temp_threshold_memoryPrecision+temp_offset XYBinLimits(2)-temp_offset2],...
    %                 [1 1].*(XYBinLimits(1)+temp_offset2),...
    %                 [1 1].*temp_crossTime_marker(end),...
    %                 'LineWidth',temp_lineWidth,'color',color_choiceOffloadHigh);
    %             hold on
    
    
    
    for temptempi=1:length(temp_crossTime_marker)
        temp_marker = temptemp_markerList(temptempi);
        
        temptemp_markerTime = temp_crossTime_marker(temptempi);
        
        temp_size = 13*1.5*1.3;%20
        if strcmp(temp_marker,'x')
            temp_size = 50*1.5;%20
        end
        temp_alpha = 1;
        temp_LineWidth = 1.75;%1.5
        
        scatter3(x(temptemp_markerTime,1),y(temptemp_markerTime,1),temptemp_markerTime,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryHigh,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(x(temptemp_markerTime,2),y(temptemp_markerTime,2),temptemp_markerTime,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffloadLow,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(x(temptemp_markerTime,3),y(temptemp_markerTime,3),temptemp_markerTime,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryLow,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(x(temptemp_markerTime,4),y(temptemp_markerTime,4),temptemp_markerTime,...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffloadHigh,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        
    end
    
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:));
    [temp_zmin,temp_zmax] = bounds(z(:));
    
    %tempSize = 3;%10
    %scatter3(x,y,z,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    %hold on
    
    %axis equal
    
    %xticks([]);
    %yticks([]);
    
    %temp_xmin = 0;
    
    temp1 = temp_xmin:(temp_xmax-temp_xmin)/4:temp_xmax;
    xticks([temp1(2) temp1(4)]);
    
    temp2 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
    yticks([temp2(2) temp2(4)]);
    
    xticklabels({'Low','High'});
    yticklabels({'Low','High'});
    
    xtickangle(0);
    ytickangle(0);
    
    %zticks([]);
    zticks(temp_crossTime_marker(1:3));
    zticklabels({'Fixation','Sample','Delay1'});
    
    %grid on
    
    set(gca,'zdir','reverse');
    %xlim([temp_xmin-(temp_xmax-temp_xmin)*0.18 temp_xmax+(temp_xmax-temp_xmin)*0.18]);%0.08
    %xlim([0 temp_xmax+(temp_xmax-temp_xmin)*0.18]);%0.08
    %xlim([0 1]);
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.12 temp_xmax+(temp_xmax-temp_xmin)*0.12]);%0.08
    
    
    %ylim([temp_ymin-(temp_ymax-temp_ymin)*0.18 temp_ymax+(temp_ymax-temp_ymin)*0.18]);
    %ylim([0 temp_ymax+(temp_ymax-temp_ymin)*0.18]);
    %ylim([0 1]);
    ylim([temp_ymin-(temp_ymax-temp_ymin)*0.12 temp_ymax+(temp_ymax-temp_ymin)*0.12]);
    
    
    %zlim([temp_zmin-(temp_zmax-temp_zmin)*0.08 temp_zmax+(temp_zmax-temp_zmin)*0.08]);
    zlim([temp_zmin-(temp_zmax-temp_zmin)*0.02 temp_zmax+(temp_zmax-temp_zmin)*0.02]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 7.5)%9
    set(gca,'box','off');% 取消右、上边框
    %xlabel('Memory precision', 'FontSize', 9, 'FontWeight', 'normal');
    %ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
    zlabel('Time', 'FontSize', 9, 'FontWeight', 'normal');
    
end

%% Plot for 3types (only length 123 in one panel)
if true
    fig = figure('Name','asd','NumberTitle','off'); %#ok<*UNRCH>
    %set(gcf,'Position',[510 90 240*0.75*0.95*1.16*0.85 240*0.75*0.95*1.16]);
    set(gcf,'Position',[510 90 240*0.75*0.95*1.16*0.85*1.08 240*0.75*0.95*1.16*0.75*1.15*0.95*0.975]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    x = [precision_trialLevel_multiPeriod_3types_length1_collapsed(1:length1_crossTime_max,:);...
        precision_trialLevel_multiPeriod_3types_length2_collapsed(1:length2_crossTime_max,:);...
        precision_trialLevel_multiPeriod_3types_length3_collapsed(1:length3_crossTime_max,:);...
        precision_trialLevel_multiPeriod_3types_length123_collapsed(1:length123_crossTime_max,:)];
    y = [meta_trialLevel_multiPeriod_3types_length1_collapsed(1:length1_crossTime_max,:);...
        meta_trialLevel_multiPeriod_3types_length2_collapsed(1:length2_crossTime_max,:);...
        meta_trialLevel_multiPeriod_3types_length3_collapsed(1:length3_crossTime_max,:);...
        meta_trialLevel_multiPeriod_3types_length123_collapsed(1:length123_crossTime_max,:)];
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:));
    
    temp_markerList = ['o','^','s','x'];
    
    %for tempi=1:4
    if true
        tempi = 4;
        if tempi == 4
            temp_strA = 'Len123';
            
            temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_length123_collapsed;
            temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_length123_collapsed;
            temp_crossTime_interval = length3_crossTime_interval;
            %temp_crossTime_max = length123_crossTime_max;
            temp_crossTime_max = 115;%89,121
            
        end
        
        %for tempPeriodIndex=1:4
        if true
            tempPeriodIndex = 4;
            nexttile
            
            if tempPeriodIndex == 4
                temp_strB = 'All';
                x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
                y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
                temp_crossTime_marker = temp_crossTime_interval;
                temp_crossTime_marker(end) = temp_crossTime_max;
                
                temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
                temptemp_markerList = temp_markerList;
            end
            
            
            %             plot(x(:,1), y(:,1), '-', 'LineWidth', 0.5, 'Color', color_choiceMemory);
            %             hold on
            %
            %             plot(x(:,2), y(:,2), '-', 'LineWidth', 0.5, 'Color', color_choiceOffload);
            %             hold on
            %
            %             plot(x(:,3), y(:,3), '-', 'LineWidth', 0.5, 'Color', color_choiceMemoryError);
            %             hold on
            
            temp_LineWidth = 0.75;%0.75
            plot(x(:,1), y(:,1), ':', 'LineWidth', temp_LineWidth, 'Color', color_choiceMemory);
            hold on
            
            plot(x(:,2), y(:,2), ':', 'LineWidth', temp_LineWidth, 'Color', color_choiceOffload);
            hold on
            
            plot(x(:,3), y(:,3), ':', 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryError);
            hold on
            
            
            %         plot(x(temp_crossTime_marker,1), y(temp_crossTime_marker,1), '-', 'LineWidth', 0.5, 'Color', color_choiceMemory);
            %         hold on
            %
            %         plot(x(temp_crossTime_marker,2), y(temp_crossTime_marker,2), '-', 'LineWidth', 0.5, 'Color', color_choiceOffload);
            %         hold on
            %
            %         plot(x(temp_crossTime_marker,3), y(temp_crossTime_marker,3), '-', 'LineWidth', 0.5, 'Color', color_choiceMemoryError);
            %         hold on
            
            
            for temptempi=1:length(temp_crossTime_marker)
                
                if temptempi == 2 || temptempi == 3
                    continue
                end
                
                temp_marker = temptemp_markerList(temptempi);
                
                temptemp_markerTime = temp_crossTime_marker(temptempi);
                
                temp_size = 13*2;%20
                if strcmp(temp_marker,'x')
                    temp_size = 50*1.3;%20
                end
                temp_alpha = 1;
                temp_LineWidth = 1.5;
                
                scatter(x(temptemp_markerTime,1),y(temptemp_markerTime,1),...
                    temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'LineWidth',temp_LineWidth,...
                    'MarkerEdgeAlpha',temp_alpha);
                hold on
                scatter(x(temptemp_markerTime,2),y(temptemp_markerTime,2),...
                    temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'LineWidth',temp_LineWidth,...
                    'MarkerEdgeAlpha',temp_alpha);
                hold on
                scatter(x(temptemp_markerTime,3),y(temptemp_markerTime,3),...
                    temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth,...
                    'MarkerEdgeAlpha',temp_alpha);
                hold on
                
            end
            
            
            xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
            ylim([temp_ymin-(temp_ymax-temp_ymin)*0.08 temp_ymax+(temp_ymax-temp_ymin)*0.08]);
            
            %xlim([temp_xmin-(temp_xmax-temp_xmin)*0.02 temp_xmax-(temp_xmax-temp_xmin)*0.25]);
            %ylim([temp_ymin+(temp_ymax-temp_ymin)*0.08 temp_ymax-(temp_ymax-temp_ymin)*0.04]);
            
            %xlim([0 1]);
            %ylim([0 1]);
            
            set(gca,'linewidth',1.5);
            set(gca, 'FontSize', 8);%10
            set(gca,'box','off');% 取消右、上边框
            xlabel(sprintf('Memory precision'), 'FontSize', 9, 'FontWeight', 'normal');%10
            ylabel(sprintf('Meta-memory'), 'FontSize', 9, 'FontWeight', 'normal');%10
            
            temp_title = title(sprintf('%s, %s',temp_strA,temp_strB), 'FontSize', 10, 'FontWeight', 'normal');
            temp_title.Interpreter = 'none';
            
        end
    end
    
end


%% Plot for 3types
if false
    fig = figure('Name','asd','NumberTitle','off'); %#ok<*UNRCH>
    set(gcf,'Position',[510 90 720 480*1.2*1.2*0.93]);
    t = tiledlayout(4,4,'TileSpacing','Compact','Padding','Compact');
    
    x = [precision_trialLevel_multiPeriod_3types_length1_collapsed(1:length1_crossTime_max,:);...
        precision_trialLevel_multiPeriod_3types_length2_collapsed(1:length2_crossTime_max,:);...
        precision_trialLevel_multiPeriod_3types_length3_collapsed(1:length3_crossTime_max,:);...
        precision_trialLevel_multiPeriod_3types_length123_collapsed(1:length123_crossTime_max,:)];
    y = [meta_trialLevel_multiPeriod_3types_length1_collapsed(1:length1_crossTime_max,:);...
        meta_trialLevel_multiPeriod_3types_length2_collapsed(1:length2_crossTime_max,:);...
        meta_trialLevel_multiPeriod_3types_length3_collapsed(1:length3_crossTime_max,:);...
        meta_trialLevel_multiPeriod_3types_length123_collapsed(1:length123_crossTime_max,:)];
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:));
    
    temp_markerList = ['o','^','s','x'];
    
    for tempi=1:4
        if tempi == 1
            temp_strA = 'Len1';
            temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_3types_length1;
            temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_3types_length1;
            
            temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_length1_collapsed;
            temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_length1_collapsed;
            temp_crossTime_interval = length1_crossTime_interval;
            %temp_crossTime_max = length1_crossTime_max_3types;
            temp_crossTime_max = length1_crossTime_max;
        elseif tempi == 2
            temp_strA = 'Len2';
            temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_3types_length2;
            temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_3types_length2;
            
            temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_length2_collapsed;
            temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_length2_collapsed;
            temp_crossTime_interval = length2_crossTime_interval;
            %temp_crossTime_max = length2_crossTime_max_3types;
            temp_crossTime_max = length2_crossTime_max;
        elseif tempi == 3
            temp_strA = 'Len3';
            temp_meta_trialLevel = meta_trialLevel_crossTimeMean_multiPeriod_3types_length3;
            temp_memoryPrecision_trialLevel = precision_trialLevel_crossTimeMean_multiPeriod_3types_length3;
            
            temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_length3_collapsed;
            temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_length3_collapsed;
            temp_crossTime_interval = length3_crossTime_interval;
            %temp_crossTime_max = length3_crossTime_max_3types;
            temp_crossTime_max = length3_crossTime_max;
        elseif tempi == 4
            temp_strA = 'Len123';
            
            temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_length123_collapsed;
            temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_length123_collapsed;
            temp_crossTime_interval = length3_crossTime_interval;
            temp_crossTime_max = length123_crossTime_max;
        end
        
        for tempPeriodIndex=1:4
            nexttile
            
            if tempPeriodIndex == 1
                temp_strB = 'Baseline';
                x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):(temp_crossTime_interval(2)),:);
                y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):(temp_crossTime_interval(2)),:);
                temp_crossTime_marker = [temp_crossTime_interval(1) (temp_crossTime_interval(2))];
                temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
                %temptemp_markerList = temp_markerList([tempPeriodIndex end]);
                temptemp_markerList = temp_markerList([tempPeriodIndex tempPeriodIndex+1]);
            elseif tempPeriodIndex == 2
                temp_strB = 'Sample';
                x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(2):(temp_crossTime_interval(3)),:);
                y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(2):(temp_crossTime_interval(3)),:);
                temp_crossTime_marker = [temp_crossTime_interval(2) (temp_crossTime_interval(3))];
                temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
                %temptemp_markerList = temp_markerList([tempPeriodIndex end]);
                temptemp_markerList = temp_markerList([tempPeriodIndex tempPeriodIndex+1]);
            elseif tempPeriodIndex == 3
                temp_strB = 'Delay1';
                %x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(3):(temp_crossTime_interval(4)),:);
                %y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(3):(temp_crossTime_interval(4)),:);
                %temp_crossTime_marker = [temp_crossTime_interval(3) temp_crossTime_interval(4)];
                
                x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(3):temp_crossTime_max,:);
                y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(3):temp_crossTime_max,:);
                temp_crossTime_marker = [temp_crossTime_interval(3) temp_crossTime_max];
                
                temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
                %temptemp_markerList = temp_markerList([tempPeriodIndex end]);
                temptemp_markerList = temp_markerList([tempPeriodIndex tempPeriodIndex+1]);
            elseif tempPeriodIndex == 4
                temp_strB = 'All';
                %x = temp_memoryPrecision_trialLevel_collapsed;
                %y = temp_meta_trialLevel_collapsed;
                %temp_crossTime_marker = temp_crossTime_interval;
                
                x = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
                y = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
                temp_crossTime_marker = temp_crossTime_interval;
                temp_crossTime_marker(end) = temp_crossTime_max;
                
                temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
                temptemp_markerList = temp_markerList;
            end
            
            
            %         plot(x(:,1), y(:,1), '-', 'LineWidth', 0.5, 'Color', color_choiceMemory);
            %         hold on
            %
            %         plot(x(:,2), y(:,2), '-', 'LineWidth', 0.5, 'Color', color_choiceOffload);
            %         hold on
            %
            %         plot(x(:,3), y(:,3), '-', 'LineWidth', 0.5, 'Color', color_choiceMemoryError);
            %         hold on
            
            
            plot(x(:,1), y(:,1), ':', 'LineWidth', 0.75, 'Color', color_choiceMemory);
            hold on
            
            plot(x(:,2), y(:,2), ':', 'LineWidth', 0.75, 'Color', color_choiceOffload);
            hold on
            
            plot(x(:,3), y(:,3), ':', 'LineWidth', 0.75, 'Color', color_choiceMemoryError);
            hold on
            
            
            %         plot(x(temp_crossTime_marker,1), y(temp_crossTime_marker,1), '-', 'LineWidth', 0.5, 'Color', color_choiceMemory);
            %         hold on
            %
            %         plot(x(temp_crossTime_marker,2), y(temp_crossTime_marker,2), '-', 'LineWidth', 0.5, 'Color', color_choiceOffload);
            %         hold on
            %
            %         plot(x(temp_crossTime_marker,3), y(temp_crossTime_marker,3), '-', 'LineWidth', 0.5, 'Color', color_choiceMemoryError);
            %         hold on
            
            
            for temptempi=1:length(temp_crossTime_marker)
                temp_marker = temptemp_markerList(temptempi);
                
                temptemp_markerTime = temp_crossTime_marker(temptempi);
                
                temp_size = 13;%20
                if strcmp(temp_marker,'x')
                    temp_size = 50;%20
                end
                if strcmp(temp_marker,'+')
                    temp_size = 20;
                end
                temp_alpha = 1;
                temp_LineWidth = 1.5;
                
                scatter(x(temptemp_markerTime,1),y(temptemp_markerTime,1),...
                    temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'LineWidth',temp_LineWidth,...
                    'MarkerEdgeAlpha',temp_alpha);
                hold on
                scatter(x(temptemp_markerTime,2),y(temptemp_markerTime,2),...
                    temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'LineWidth',temp_LineWidth,...
                    'MarkerEdgeAlpha',temp_alpha);
                hold on
                scatter(x(temptemp_markerTime,3),y(temptemp_markerTime,3),...
                    temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth,...
                    'MarkerEdgeAlpha',temp_alpha);
                hold on
                
            end
            
            
            xlim([temp_xmin-(temp_xmax-temp_xmin)*0.08 temp_xmax+(temp_xmax-temp_xmin)*0.08]);
            ylim([temp_ymin-(temp_ymax-temp_ymin)*0.08 temp_ymax+(temp_ymax-temp_ymin)*0.08]);
            %ylim([0 1]);
            
            set(gca,'linewidth',1.5)
            set(gca, 'FontSize', 10)
            set(gca,'box','off');% 取消右、上边框
            xlabel(sprintf('Memory precision'), 'FontSize', 10, 'FontWeight', 'normal');
            ylabel(sprintf('Meta-memory'), 'FontSize', 10, 'FontWeight', 'normal');
            
            temp_title = title(sprintf('%s, %s',temp_strA,temp_strB), 'FontSize', 10, 'FontWeight', 'normal');
            temp_title.Interpreter = 'none';
            
        end
    end
    
end

%% Plot for single trial of choice-memory error (3d), x is time
temp1 = meta_eachTrial_multiPeriod_3types_length123_collapsed{3};
temp2 = precision_eachTrial_multiPeriod_3types_length123_collapsed{3};

temptempBoolIndex = (~isnan(temp1(1,:))) & (~isnan(temp2(1,:)));

temp1 = temp1(:,temptempBoolIndex);
temp2 = temp2(:,temptempBoolIndex);

if false
    fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
    set(gcf,'Position',[10 150 450*0.8*0.8*0.8 450*0.8*0.8*0.8*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    t = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');
    
    temp_strA = 'Len123';
    temp_meta_trialLevel_collapsed = temp1;
    temp_memoryPrecision_trialLevel_collapsed = temp2;
    temp_crossTime_interval = length3_crossTime_interval;
    temp_crossTime_max = length123_crossTime_max;
    
    
    temp_strB = 'All';
    y = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    z = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    temptemp_markerList = temp_markerList;
    
    x = 1:temp_crossTime_marker(end);
    
    windowSize = 8;%8
    y = smoothdata(y,1,'gaussian',windowSize);
    z = smoothdata(z,1,'gaussian',windowSize);
    
    temp_LineWidth = 0.5;%0.75
    temp_Linestyle = '-';
    
    for tempi=1:size(y,2)
        temp_color = color_choiceMemoryError;
        
        plot3(x, y(:,tempi),z(:,tempi), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', temp_color);
        hold on
    end
    
    
    
    XYBinLimits = [0 1];
    
    temp_offset = 0.003;%0.02
    temp_offset2 = 0.015;%0.015
    temp_lineWidth = 1;%4,1.5
    
    temp_threshold_memoryPrecision = temp_precision_mean_length123;
    temp_threshold_meta = temp_meta_mean_length123;
    
    temp_facealpha = 0.9;%0.3
    
    % high-match
    temp2 = [temp_threshold_memoryPrecision temp_threshold_memoryPrecision 1 1];
    temp3 = [temp_threshold_meta 1 1 temp_threshold_meta];
    temp1 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceMemoryHigh_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % under-mismatch
    temp2 = [temp_threshold_memoryPrecision temp_threshold_memoryPrecision 1 1];
    temp3 = [0 temp_threshold_meta temp_threshold_meta 0];
    temp1 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceOffloadHigh_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % low-match
    temp2 = [0 0 temp_threshold_memoryPrecision temp_threshold_memoryPrecision];
    temp3 = [0 temp_threshold_meta temp_threshold_meta 0];
    temp1 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceOffloadLow_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    % over-mismatch
    temp2 = [0 0 temp_threshold_memoryPrecision temp_threshold_memoryPrecision];
    temp3 = [temp_threshold_meta 1 1 temp_threshold_meta];
    temp1 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,color_choiceMemoryLow_B);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:));
    [temp_zmin,temp_zmax] = bounds(z(:));
    
    temp1 = temp_ymin:(temp_ymax-temp_ymin)/4:temp_ymax;
    yticks([temp1(2) temp1(4)]);
    
    temp2 = temp_zmin:(temp_zmax-temp_zmin)/4:temp_zmax;
    zticks([temp2(2) temp2(4)]);
    
    yticklabels({'Low','High'});
    zticklabels({'Low','High'});
    
    xtickangle(0);
    ytickangle(0);
    ztickangle(0);
    
    xticks(temp_crossTime_marker(1:3));
    xticklabels({'Fixation','Sample','Delay1'});
    
    set(gca,'ydir','reverse');
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.02 temp_xmax+(temp_xmax-temp_xmin)*0.02]);%0.08
    ylim([temp_ymin-(temp_ymax-temp_ymin)*0.12 temp_ymax+(temp_ymax-temp_ymin)*0.12]);
    zlim([temp_zmin-(temp_zmax-temp_zmin)*0.12 temp_zmax+(temp_zmax-temp_zmin)*0.12]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 7.5)%9
    set(gca,'box','off');% 取消右、上边框
    xlabel('Time', 'FontSize', 9, 'FontWeight', 'normal');
    
end



%% Plot for 3types (3d), x is time
if true
    fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
    %set(gcf,'Position',[10 150 450*0.8*0.8*0.8*0.95 450*0.8*0.8*0.8*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    set(gcf,'Position',[10 550 450*0.8*0.8*0.8*0.95*1.63*1.08*1.04*0.9*0.95 450*0.8*0.8*0.8*0.8*0.9*0.9*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    t = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');
    
    temp_strA = 'Len123';
    temp_meta_trialLevel_collapsed = meta_trialLevel_multiPeriod_3types_length123_collapsed;
    temp_memoryPrecision_trialLevel_collapsed = precision_trialLevel_multiPeriod_3types_length123_collapsed;
    temp_crossTime_interval = length3_crossTime_interval;
    %temp_crossTime_max = length123_crossTime_max;    
    temp_crossTime_max = 121;%89,121,94,115
    
    
    temp_strB = 'All';
    y = temp_memoryPrecision_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    z = temp_meta_trialLevel_collapsed(temp_crossTime_interval(1):temp_crossTime_max,:);
    temp_crossTime_marker = temp_crossTime_interval;
    temp_crossTime_marker(end) = temp_crossTime_max;
    
    temp_crossTime_marker = temp_crossTime_marker - temp_crossTime_marker(1) + 1;
    temptemp_markerList = temp_markerList;
    
    x = 1:temp_crossTime_marker(end);
    
    windowSize = 8;%8
    y = smoothdata(y,1,'gaussian',windowSize);
    z = smoothdata(z,1,'gaussian',windowSize);
    
    temp_LineWidth = 0.75;
    temp_Linestyle = '-';
    
    plot3(x, y(:,1),z(:,1), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemory);
    hold on
    
    plot3(x, y(:,2),z(:,2), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceOffload);
    hold on
    
    plot3(x, y(:,3),z(:,3), temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', color_choiceMemoryError);
    hold on
    
    plot3([1 1]*temp_crossTime_marker(end), [y(temp_crossTime_marker(end),1) y(temp_crossTime_marker(end),2)]-0.01,...
        [z(temp_crossTime_marker(end),1) z(temp_crossTime_marker(end),2)], temp_Linestyle, 'LineWidth', temp_LineWidth, 'Color', [1 1 1]*0.5);
    hold on
    
    
    XYBinLimits = [0 1];
    
    temp_offset = 0.003;%0.02
    temp_offset2 = 0.015;%0.015
    temp_lineWidth = 1;%4,1.5
    
    temp_threshold_memoryPrecision = temp_precision_mean_length123;
    temp_threshold_meta = temp_meta_mean_length123;
    
    temp_facealpha = 0.9;%0.3
    
    temp2 = [0 0 1 1];
    temp3 = [0 1 1 0];
    temp1 = [1 1 1 1].*temp_crossTime_marker(end);
    h = fill3(temp1,temp2,temp3,[1 1 1]*0.9);
    set(h,'edgealpha',0,'facealpha',temp_facealpha)
    
    
    for temptempi=1:length(temp_crossTime_marker)
        
        if temptempi==3
            continue
        end
        
        temp_marker = temptemp_markerList(temptempi);
        
        temptemp_markerTime = temp_crossTime_marker(temptempi);
        
        temp_size = 13*1.5*1.3;%20
        if strcmp(temp_marker,'x')
            temp_size = 50*1.5;%20
        end
        temp_alpha = 1;
        temp_LineWidth = 1.75;%1.5
        
        scatter3(temptemp_markerTime,y(temptemp_markerTime,1),z(temptemp_markerTime,1),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemory,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temptemp_markerTime,y(temptemp_markerTime,2),z(temptemp_markerTime,2),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceOffload,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        scatter3(temptemp_markerTime,y(temptemp_markerTime,3),z(temptemp_markerTime,3),...
            temp_size,temp_marker,'MarkerEdgeColor',color_choiceMemoryError,'LineWidth',temp_LineWidth,...
            'MarkerEdgeAlpha',temp_alpha);
        hold on
        
    end
    
    
    [temp_xmin,temp_xmax] = bounds(x(:));
    [temp_ymin,temp_ymax] = bounds(y(:));
    [temp_zmin,temp_zmax] = bounds(z(:));
    
    
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
    xticklabels({'Fixation','Sample','Delay1'});
    
    %grid on
    
    set(gca,'ydir','reverse');
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.02 temp_xmax+(temp_xmax-temp_xmin)*0.02]);%0.08
    ylim([temp_ymin temp_ymax+(temp_ymax-temp_ymin)*0.12]);
    zlim([temp_zmin-(temp_zmax-temp_zmin)*0.12 temp_zmax+(temp_zmax-temp_zmin)*0.12]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 7.5)%9
    set(gca,'box','off');% 取消右、上边框
    %xlabel('Memory precision', 'FontSize', 9, 'FontWeight', 'normal');
    %ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
    xlabel('Time', 'FontSize', 9, 'FontWeight', 'normal');
    
end





%% meta_trialLevel_crossTime_lastT
lengthx_sample_interval = length1_sample_interval;
temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;


temp_meta_trialLevel_crossTime = nan(size(meta_trialLevel_crossTime_length1));

for temptempi=1:size(temp_meta_trialLevel_crossTime,1)
    if trialBoolIndex_length1(temptempi) == true
        temp1 = meta_trialLevel_crossTime_length1(temptempi,:);
        if isnan(temp1(1)) == true
            continue
        end
        temp_meta_trialLevel_crossTime(temptempi,:) = temp1;
    end
    if trialBoolIndex_length2(temptempi) == true
        temp1 = meta_trialLevel_crossTime_length2(temptempi,:);
        if isnan(temp1(1)) == true
            continue
        end
        temp1 = temp1(end-lengthx_sample_interval(end)+1:end);
        temp_meta_trialLevel_crossTime(temptempi,:) = temp1;
    end
    if trialBoolIndex_length3(temptempi) == true
        temp1 = meta_trialLevel_crossTime_length3(temptempi,:);
        if isnan(temp1(1)) == true
            continue
        end
        temp1 = temp1(end-lengthx_sample_interval(end)+1:end);
        temp_meta_trialLevel_crossTime(temptempi,:) = temp1;
    end
    
end
temp_meta_trialLevel_crossTime;

meta_trialLevel_crossTime_lastT = temp_meta_trialLevel_crossTime;

meta_trialLevel_crossTime_lastTDelay1 = [meta_trialLevel_crossTime_lastT,meta_trialLevel_crossTime_delay1];

%% memoryPrecision_trialLevel_crossTime_lastT
lengthx_sample_interval = length1_sample_interval;
temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;


temp_memoryPrecision_trialLevel_crossTime = nan(size(memoryPrecision_trialLevel_crossTime_length1));

for temptempi=1:size(temp_memoryPrecision_trialLevel_crossTime,1)
    if trialBoolIndex_length1(temptempi) == true
        temp1 = memoryPrecision_trialLevel_crossTime_length1(temptempi,:);
        if isnan(temp1(1)) == true
            continue
        end
        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1;
    end
    if trialBoolIndex_length2(temptempi) == true
        temp1 = memoryPrecision_trialLevel_crossTime_length2(temptempi,:);
        if isnan(temp1(1)) == true
            continue
        end
        temp1 = temp1(end-lengthx_sample_interval(end)+1:end);
        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1;
    end
    if trialBoolIndex_length3(temptempi) == true
        temp1 = memoryPrecision_trialLevel_crossTime_length3(temptempi,:);
        if isnan(temp1(1)) == true
            continue
        end
        temp1 = temp1(end-lengthx_sample_interval(end)+1:end);
        temp_memoryPrecision_trialLevel_crossTime(temptempi,:) = temp1;
    end
    
end
temp_memoryPrecision_trialLevel_crossTime;

memoryPrecision_trialLevel_crossTime_lastT = temp_memoryPrecision_trialLevel_crossTime;

memoryPrecision_trialLevel_crossTime_lastTDelay1 = [memoryPrecision_trialLevel_crossTime_lastT,memoryPrecision_trialLevel_crossTime_delay1];

%% Cross correlation
memoryPrecision_trialLevel_crossTime_lastTDelay1;
meta_trialLevel_crossTime_lastTDelay1;

if_cross_trial0_seq1 = 0;

temptempBoolIndex = choiceBoolIndex;
% temptempBoolIndex = choiceMemoryBoolIndex;
% temptempBoolIndex = choiceMemoryCorrectBoolIndex;
% temptempBoolIndex = choiceMemoryErrorBoolIndex;
% temptempBoolIndex = choiceOffloadBoolIndex;


% temp1 = memoryPrecision_trialLevel_crossTime_lastTDelay1;
% temp2 = meta_trialLevel_crossTime_lastTDelay1;
% temp2 = memoryPrecision_trialLevel_crossTime_lastTDelay1;
% temp2(:,1:50) = temp1(:,11:60);%temp1 slower than temp2
% temp2(:,11:60) = temp1(:,1:50);%temp1 faster than temp2

% temp1 = memoryPrecision_trialLevel_crossTime_lastT;
% temp2 = meta_trialLevel_crossTime_lastT;
% temp2 = memoryPrecision_trialLevel_crossTime_lastT;
% temp2(:,1:14) = temp1(:,5:18);%temp1 slower than temp2
% temp2(:,5:18) = temp1(:,1:14);%temp1 faster than temp2

% temp1 = memoryPrecision_trialLevel_crossTime_delay1;
% temp2 = meta_trialLevel_crossTime_delay1;
% temp2 = memoryPrecision_trialLevel_crossTime_delay1;
% temp2(:,1:36) = temp1(:,11:46);%temp1 slower than temp2
% temp2(:,11:46) = temp1(:,1:36);%temp1 faster than temp2


% temp1 = memoryPrecision_trialLevel_crossTime_length1;
% temp2 = meta_trialLevel_crossTime_length1;
% temp2 = memoryPrecision_trialLevel_crossTime_length1;
% temp2(:,1:14) = temp1(:,5:18);%temp1 slower than temp2
% temp2(:,5:18) = temp1(:,1:14);%temp1 faster than temp2

% temp1 = memoryPrecision_trialLevel_crossTime_length2;
% temp2 = meta_trialLevel_crossTime_length2;
% temp2 = memoryPrecision_trialLevel_crossTime_length2;
% temp2(:,1:14) = temp1(:,5:18);%temp1 slower than temp2
% temp2(:,5:18) = temp1(:,1:14);%temp1 faster than temp2

temp1 = memoryPrecision_trialLevel_crossTime_length3;
temp2 = meta_trialLevel_crossTime_length3;
% temp2 = memoryPrecision_trialLevel_crossTime_length3;
% temp2(:,1:14) = temp1(:,5:18);%temp1 slower than temp2
% temp2(:,5:18) = temp1(:,1:14);%temp1 faster than temp2

% temp1 = memoryPrecision_trialLevel_crossTime_baseline;
% temp2 = meta_trialLevel_crossTime_baseline;
% temp2 = memoryPrecision_trialLevel_crossTime_baseline;
% temp2(:,1:14) = temp1(:,5:18);%temp1 slower than temp2
% temp2(:,5:18) = temp1(:,1:14);%temp1 faster than temp2


temp1_seqLevel = zeros(sum(numSeq(1:3)),size(temp1,2));
for temptempi=1:sum(numSeq(1:3))
    temptempBoolIndex_B = temp_seqIndex==temptempi;
    %temptempBoolIndex_C = temptempBoolIndex_B;
    temptempBoolIndex_C = temptempBoolIndex_B & temptempBoolIndex;
    temp1_seqLevel(temptempi,:) = mean(temp1(temptempBoolIndex_C,:),'omitnan');
end

temp2_seqLevel = zeros(sum(numSeq(1:3)),size(temp2,2));
for temptempi=1:sum(numSeq(1:3))
    temptempBoolIndex_B = temp_seqIndex==temptempi;
    %temptempBoolIndex_C = temptempBoolIndex_B;
    temptempBoolIndex_C = temptempBoolIndex_B & temptempBoolIndex;
    temp2_seqLevel(temptempi,:) = mean(temp2(temptempBoolIndex_C,:),'omitnan');
end


temp1 = temp1(temptempBoolIndex,:);
temp2 = temp2(temptempBoolIndex,:);

if if_cross_trial0_seq1 == 0
    
elseif if_cross_trial0_seq1 == 1
    temp1 = temp1_seqLevel;
    temp2 = temp2_seqLevel;
end

% crossCorr > 0 means temp1 slower than temp2
% crossCorr < 0 means temp1 faster than temp2
r_crossCorr = nan(size(temp1,1),2*size(temp1,2)-1);
for tempi=1:size(temp1,1)
    temp_r_crossCorr = xcorr(temp1(tempi,:),temp2(tempi,:),'normalized');
    %temp_r_crossCorr = xcorr(temp1(tempi,:),temp1(tempi,:),'normalized');
    r_crossCorr(tempi,:) = temp_r_crossCorr;
end
r_crossCorr_valid = r_crossCorr(~isnan(r_crossCorr(:,1)),:);


r_crossCorr_valid_mean = mean(r_crossCorr_valid,1);
r_crossCorr_valid_sem = std(r_crossCorr_valid,1,1)./sqrt(size(r_crossCorr_valid,1));
% r_crossCorr_valid_sem = std(r_crossCorr_valid,1,1);

interval_crossCorr = [1,size(temp1,2),2*size(temp1,2)-1];

interval_crossCorr_time = 1000*(interval_crossCorr-interval_crossCorr(2))/30;
interval_crossCorr_time_str = [{num2str(interval_crossCorr_time(1))} {num2str(interval_crossCorr_time(2))} {num2str(interval_crossCorr_time(3))}];

%% Plot
fig = figure('Name','asd','NumberTitle','off');

temptemp1 = 360*0.975*0.74*1.37*0.97*1.02;
temptemp2 = 200*0.8*0.8;

set(gcf,'Position',[10+400*(7-4) 50+600 temptemp1 temptemp2]);

t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

nexttile

temp1_mean = r_crossCorr_valid_mean;
temp1_sem = r_crossCorr_valid_sem;

y_min = min([temp1_mean-temp1_sem]);
y_max = max([temp1_mean+temp1_sem]);

temp_interval = interval_crossCorr;

x = 1:length(r_crossCorr_valid_mean);

[x_min,x_max] = bounds(x);

y = temp1_mean;
plot(x,y,'color',[0.3 0.3 0.3],'linewidth',1);
hold on
y_sem = temp1_sem;
patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.3 0.3 0.3],'FaceAlpha',0.1,'EdgeColor',[0.62 0.62 0.62]);
hold on

[M,I] = max(y);
x_optimal = x(I);
y_optimal = M;

time_optimal = 1000*(x_optimal-interval_crossCorr(2))/30;

plot([1 1].*x_optimal,[y_min-(y_max-y_min)*0.1 y_optimal],':','color',[0.3 0.3 0.3],'linewidth',1);
hold on

text(x_optimal+(x_max-x_min)*0.02,y_optimal-(y_max-y_min)*0.4,sprintf('%.0f ms',time_optimal), 'fontsize',8,'FontWeight','normal');


set(gca,'linewidth',1.5)
set(gca, 'FontSize', 8)%12
set(gca,'color',backgounrdColor);
set(gca,'box','off');% 取消右、上边框
xtickangle(0);

% temp_interval_all = [baselinePeriod_interval(1) ...
%     (baselinePeriod_interval(end) + lengthx_sample_interval(1:end-1)) ...
%     (baselinePeriod_interval(end) + lengthx_sample_interval(end) + decisionPeriodA_interval)];

% xticks(temp_interval_all(1:end-2));
xticks(temp_interval);
% xticklabels({'-2.2','0','2.2'});
xticklabels(interval_crossCorr_time_str);

xlim([temp_interval(1) temp_interval(end)]);
ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.1]);

%set(gca,'ydir','reverse');
xlabel('Time (s)', 'FontSize', 9, 'FontWeight', 'normal');%10
ylabel('Cross-correlation', 'FontSize', 9, 'FontWeight', 'normal');%10
% title('Memory precision VS. accuracy','FontSize',9);




%% End