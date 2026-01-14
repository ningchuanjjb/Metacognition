%This script: To conduct behavioral analysis
%% Initialization
% clear all %#ok<CLALL>
close all

%% Some flags to control the following script

% if_monkey_D0_Z1 = 0;% To decide whether dealing with Ding's data or Zelku's data

ifSaveOffloadingProb = 0;
ifSave_gAcc = 0;
ifSave_endingHold = 0;
ifSaveMapping = 0;

ifSave_paperFig2Data = 0;

ifSave_fig = zeros(1, 27);
% ifSave_fig = ones(1, 17);
% ifSave_fig(17) = 1;
% ifSave_fig(1)  = 1;

% ifSave_fig(1)  = 1;
% ifSave_fig(2)  = 1;
% ifSave_fig(3)  = 1;
% ifSave_fig(4)  = 1;
% ifSave_fig(5)  = 1;
% ifSave_fig(6)  = 0;
% ifSave_fig(7)  = 0;
% ifSave_fig(8)  = 1;
% ifSave_fig(9)  = 1;
% ifSave_fig(10) = 1;
% ifSave_fig(11) = 1;
% ifSave_fig(12) = 1;
% ifSave_fig(13) = 1;
% ifSave_fig(14) = 1;
% ifSave_fig(15) = 1;
% ifSave_fig(17) = 1;
% ifSave_fig(18) = 1;
% ifSave_fig(19) = 1;
% ifSave_fig(20) = 1;
% ifSave_fig(21) = 1;
% ifSave_fig(22) = 1;
% ifSave_fig(26) = 1;
% ifSave_fig(27) = 1;

path_preprocessed = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\preprocessed';


if if_monkey_D0_Z1 == 0
    monkey_name = 'D';
    
    % searchName = 'preprocessed_DEYECopy1_of_';
    
    %searchName = 'preprocessed_DnumFrames6_A';
%     searchName = 'preprocessed_DnumFrames8_A';    
%     fileMergeNum = 4;
%     searchName = 'preprocessed_DnumFrames7_A';
%     fileMergeNum = 1;    
%     searchName = 'preprocessed_DnumFrames6_Train1';
    %searchName = 'preprocessed_DnumFrames6_Train2';
    
    searchName = 'preprocessed_DnumFrames6_2p';
    fileMergeNum = 3;   
    
elseif if_monkey_D0_Z1 == 1
    monkey_name = 'Z';
    
%     searchName = 'preprocessed_ZnumFrames6_A';  
%     fileMergeNum = 3;    

%     searchName = 'preprocessed_ZnumFrames7_A';  
%     fileMergeNum = 1;    
    
%     searchName = 'preprocessed_ZnumFrames6_B2';
%     searchName = 'preprocessed_ZnumFrames6_Train1';
%     searchName = 'preprocessed_ZnumFrames6_Train2';
%     searchName = 'preprocessed_ZnumFrames6_Train3';
%     searchName = 'preprocessed_ZnumFrames6_C';  
%     searchName = 'preprocessed_ZnumFrames6_D';  
%     searchName = 'preprocessed_ZnumFrames6_early2p';
    searchName = 'preprocessed_ZnumFrames6_all2p';
%     searchName = 'preprocessed_ZnumFrames6_almost2p';
    fileMergeNum = 3;    
end

%% To load data from files
% Load all sessions' data files
[MAT_file_load, MAT_file_name_raw] = loadMat_multi(searchName, path_preprocessed);
fileSize = length(MAT_file_name_raw);

MAT_file_name = MAT_file_name_raw;
for tempi=1:fileSize
    %MAT_file_name{tempi} = MAT_file_name_raw{tempi};
    temp_slash = strfind(MAT_file_name_raw{tempi},'\');
    MAT_file_name{tempi} = MAT_file_name_raw{tempi}(temp_slash(end)+1:end);
end


fileMerged_size = floor(fileSize/fileMergeNum);

scatter_hslColor = zeros(fileMerged_size, 3);
scatter_hslColor(:,2) = 0.5;
scatter_hslColor(:,3) = 0.5;
for tempi=1:fileMerged_size
    scatter_hslColor(tempi,1) = tempi/fileMerged_size;    
end
scatter_rgbColor = hsl2rgb(scatter_hslColor);


%% Get data from each file
initialization = 0;
initialization_2 = 0;
subject_name = cell(1, fileSize);
trial_num_file = zeros(1, fileSize);

if if_monkey_D0_Z1 == 1
    trial_num_valid = 63*3*2;
end
for file_index=1:fileSize
    temp_load = eval(['MAT_file_load.file',sprintf('%d',file_index)]); 
    
    
    file_name = cell2mat(MAT_file_name(file_index));
    temp_tail = strfind(file_name,'train');
    subject_name(file_index) = {file_name(length(searchName)+1+6: temp_tail-1)};
    
    monkey_name = file_name(temp_tail+4+1);
    
    
    
    trial_num = temp_load.trial_para.trial_count;
    
    if if_monkey_D0_Z1 == 1
        %trial_num = min([trial_num,trial_num_valid]);
    end
    
    trial_num_file(file_index) = temp_load.trial_para.trial_count;
    
    currentSequence = temp_load.trial_para.currentSequence;

    
    seq_length_rangeHead = temp_load.basic_para.seq_length_rangeHead;
    seq_length_rangeTail = temp_load.basic_para.seq_length_rangeTail;
    pointKindsNum = seq_length_rangeTail-seq_length_rangeHead+1;
    numFrames = temp_load.basic_para.numFrames_rangeHead;


    
    if initialization == 0
        target_trialNum = zeros(fileSize, pointKindsNum);
        target_seq = cell(fileSize, pointKindsNum);
        target_trialIndex = cell(fileSize, pointKindsNum);
        target_numSeq = 0;
        target_seqSet = cell(1, pointKindsNum);
        target_seqSetIndex = cell(fileSize, pointKindsNum);
        target_seqSetNum = cell(fileSize, pointKindsNum);
        
        offloading_trial_count_newSeq = cell(fileSize, pointKindsNum);
        offloadingCorrect_trial_count_newSeq = cell(fileSize, pointKindsNum);
        internal_trial_count_newSeq = cell(fileSize, pointKindsNum);
        offloadingProb = cell(fileSize, pointKindsNum);
        
        gCorrect_trial_count = cell(fileSize, pointKindsNum);
        g_trial_count = cell(fileSize, pointKindsNum);
        gAcc = cell(fileSize, pointKindsNum);
        
        gCorrect_choice_trial_count = cell(fileSize, pointKindsNum);
        g_choice_trial_count = cell(fileSize, pointKindsNum);
        gAcc_choice = cell(fileSize, pointKindsNum);
        
        gCorrect_noChoice_trial_count = cell(fileSize, pointKindsNum);
        g_noChoice_trial_count = cell(fileSize, pointKindsNum);
        gAcc_noChoice = cell(fileSize, pointKindsNum);
        
        gCorrect_noChoice_trial_count_inOrder = cell(fileSize, pointKindsNum);
        g_noChoice_trial_count_inOrder = cell(fileSize, pointKindsNum);        
        gAcc_noChoice_inOrder = cell(fileSize, pointKindsNum);
        
        gDistri_noChoice_trial_count = cell(fileSize, pointKindsNum);
        gDistri_choice_trial_count = cell(fileSize, pointKindsNum);
        gDistri_length_noChoice_trial_count = cell(fileSize, pointKindsNum);
        gDistri_length_choice_trial_count = cell(fileSize, pointKindsNum);
        gDistri_lengthDetail_noChoice_trial_count = cell(fileSize, pointKindsNum);
        gDistri_lengthDetail_choice_g_trial_count = cell(fileSize, pointKindsNum);        
        
        endingHold_correct = cell(fileSize, pointKindsNum);
        endingHold_error = cell(fileSize, pointKindsNum);

        
        selecting_RT = cell(fileSize, 1);
        
        delay1_file = cell(fileSize, 1);
        delay2_file = cell(fileSize, 1);
        delay12_file = cell(fileSize, 1);
        isCorrect_file = cell(fileSize, 1);
        ifSelectOffloading_file = cell(fileSize, 1);
        choiceCondition_flag_file = cell(fileSize, 1);
        
        trial_num_length_file = zeros(fileSize, pointKindsNum);

        
        initialization = 1;
    end
    
    trial_num_length_file(file_index, :); %#ok<NOEFF>
    
    for tempi=1:pointKindsNum
        trial_num_length_file(file_index, tempi) = sum(temp_load.trial_para.seq_length==tempi);                
    end
    
    selecting_RT{file_index} = temp_load.trial_para.selecting_RT;
    
    delay1_file{file_index} = temp_load.basic_para.fixationTime;
    delay2_file{file_index} = temp_load.basic_para.twoFixationTime;
    delay12_file{file_index} = temp_load.basic_para.fixationTime + temp_load.basic_para.twoFixationTime;
    isCorrect_file{file_index} = temp_load.trial_para.isCorrect;
    ifSelectOffloading_file{file_index} = temp_load.trial_para.ifSelectOffloading;
    choiceCondition_flag_file{file_index} = temp_load.trial_para.choiceCondition_flag;
    
    for target_seqLength=1:pointKindsNum
        

        for tempi=1:trial_num
            if length(currentSequence{tempi}) == target_seqLength
                target_trialNum(file_index, target_seqLength) = target_trialNum(file_index, target_seqLength) + 1;
            end
        end
        target_seq{file_index, target_seqLength} = cell(1, target_trialNum(file_index, target_seqLength));
        target_trialIndex{file_index, target_seqLength} = zeros(1, target_trialNum(file_index, target_seqLength));
        tempj = 0;
        for tempi=1:trial_num
            if length(currentSequence{tempi}) == target_seqLength
                tempj = tempj + 1;
                target_seq{file_index, target_seqLength, tempj} = currentSequence{tempi};
                target_trialIndex{file_index, target_seqLength, tempj} = tempi;
            end
        end
        
        if target_seqLength == 1
            target_numSeq = nchoosek(temp_load.basic_para.numFrames_rangeTail, 1);
        elseif target_seqLength == 2
            target_numSeq = nchoosek(temp_load.basic_para.numFrames_rangeTail, 2);
        elseif target_seqLength == 3
            target_numSeq = nchoosek(temp_load.basic_para.numFrames_rangeTail, 3);
        elseif target_seqLength == 4
            target_numSeq = nchoosek(temp_load.basic_para.numFrames_rangeTail, 4);
        end
        
        target_seqSet{target_seqLength}= cell(target_numSeq, 1);
        temp_index = 0;
        if target_seqLength == 1
            for tempi=1:temp_load.basic_para.numFrames_rangeTail+1-target_seqLength
                temp_index = temp_index + 1;
                target_seqSet{target_seqLength}{temp_index} = [tempi]; %#ok<*NBRAK>
            end
        elseif target_seqLength == 2
            for tempi=1:temp_load.basic_para.numFrames_rangeTail+1-target_seqLength
                for tempj=tempi+1:temp_load.basic_para.numFrames_rangeTail+1-target_seqLength+1
                    temp_index = temp_index + 1;
                    target_seqSet{target_seqLength}{temp_index} = [tempi tempj];
                end
            end
        elseif target_seqLength == 3
            for tempi=1:temp_load.basic_para.numFrames_rangeTail+1-target_seqLength
                for tempj=tempi+1:temp_load.basic_para.numFrames_rangeTail+1-target_seqLength+1
                    for tempk=tempj+1:temp_load.basic_para.numFrames_rangeTail+1-target_seqLength+2
                        temp_index = temp_index + 1;
                        target_seqSet{target_seqLength}{temp_index} = [tempi tempj tempk];
                    end
                end
            end
            
        elseif target_seqLength == 4
            for tempi=1:temp_load.basic_para.numFrames_rangeTail+1-target_seqLength
                for tempj=tempi+1:temp_load.basic_para.numFrames_rangeTail+1-target_seqLength+1
                    for tempk=tempj+1:temp_load.basic_para.numFrames_rangeTail+1-target_seqLength+2
                        for tempz=tempk+1:temp_load.basic_para.numFrames_rangeTail+1-target_seqLength+3
                            temp_index = temp_index + 1;
                            target_seqSet{target_seqLength}{temp_index} = [tempi tempj tempk tempz];
                        end
                    end
                end
            end
        end
        
        
        
        target_seqSetIndex{file_index, target_seqLength} = cell(target_numSeq, 1);
        for tempi=1:target_trialNum(file_index, target_seqLength)
            for tempj=1:target_numSeq
                if target_seq{file_index, target_seqLength, tempi} == target_seqSet{target_seqLength}{tempj}
                    target_seqSetIndex{file_index, target_seqLength}{tempj} =...
                        [target_seqSetIndex{file_index, target_seqLength}{tempj} ...
                        target_trialIndex{file_index, target_seqLength, tempi}];
                end
            end
        end
        %2020-11-05 mark
        %2020-11-06 mark
        
        
    end
    
    if initialization_2 == 0
        numSeq = zeros(1, pointKindsNum);
        for tempi=1:pointKindsNum
            numSeq(tempi) = size(target_seqSet{1, tempi}, 1);
        end
        
        endingHold = struct;
        endingHold.offload = struct;
        endingHold.offload.all = cell(1, pointKindsNum);
        endingHold.offload.correct = cell(1, pointKindsNum);
        endingHold.offload.error = cell(1, pointKindsNum);
        endingHold.choiceMemory = struct;
        endingHold.choiceMemory.all = cell(1, pointKindsNum);
        endingHold.choiceMemory.correct = cell(1, pointKindsNum);
        endingHold.choiceMemory.error = cell(1, pointKindsNum);        
        endingHold.noChoiceMemory = struct;
        endingHold.noChoiceMemory.all = cell(1, pointKindsNum);
        endingHold.noChoiceMemory.correct = cell(1, pointKindsNum);
        endingHold.noChoiceMemory.error = cell(1, pointKindsNum);    
        
        for tempi=1:pointKindsNum
            endingHold.offload.all{tempi} = cell(1, numSeq(tempi));
            endingHold.offload.correct{tempi} = cell(1, numSeq(tempi));
            endingHold.offload.error{tempi} = cell(1, numSeq(tempi));
            endingHold.choiceMemory.all{tempi} = cell(1, numSeq(tempi));
            endingHold.choiceMemory.correct{tempi} = cell(1, numSeq(tempi));
            endingHold.choiceMemory.error{tempi} = cell(1, numSeq(tempi));
            endingHold.noChoiceMemory.all{tempi} = cell(1, numSeq(tempi));
            endingHold.noChoiceMemory.correct{tempi} = cell(1, numSeq(tempi));
            endingHold.noChoiceMemory.error{tempi} = cell(1, numSeq(tempi));            
        end        

        initialization_2 = 1;
    end
       
    internal_trial_count_newSeq(file_index, :) = target_seqSet;
    offloading_trial_count_newSeq(file_index, :) = target_seqSet;
    offloadingCorrect_trial_count_newSeq(file_index, :) = target_seqSet;
    offloadingProb(file_index, :) = target_seqSet;
    for target_seqLength=1:pointKindsNum
        current_numSeq = length(target_seqSet{target_seqLength});
        
        internal_trial_count_newSeq{file_index, target_seqLength} = zeros(current_numSeq, 1);
        offloading_trial_count_newSeq{file_index, target_seqLength} = zeros(current_numSeq, 1);
        offloadingCorrect_trial_count_newSeq{file_index, target_seqLength} = zeros(current_numSeq, 1);        
        offloadingProb{file_index, target_seqLength} = zeros(current_numSeq, 1);

        for tempi=1:current_numSeq
            for tempj=1:length(target_seqSetIndex{file_index, target_seqLength}{tempi})
                
                trial_count = target_seqSetIndex{file_index, target_seqLength}{tempi}(tempj);
                SL = temp_load.basic_para.cumulativeCorrectLimit_seqLengthSwitch;
                if temp_load.trial_para.selectFlag_newSeq(trial_count) == 1 && ...
                        temp_load.trial_para.choiceCondition_flag(trial_count) == 2 &&...
                        floor((trial_count-1)/SL)*SL == trial_count-1
                    
                    if temp_load.trial_para.ifSelectOffloading(trial_count) == 1
                        offloading_trial_count_newSeq{file_index, target_seqLength}(tempi) = ...
                            offloading_trial_count_newSeq{file_index, target_seqLength}(tempi) + 1;
                        
                        if temp_load.trial_para.isCorrect(trial_count) == 1
                            offloadingCorrect_trial_count_newSeq{file_index, target_seqLength}(tempi) = ...
                                offloadingCorrect_trial_count_newSeq{file_index, target_seqLength}(tempi) + 1;
                        end
                    elseif temp_load.trial_para.ifSelectOffloading(trial_count) == -1
                        internal_trial_count_newSeq{file_index, target_seqLength}(tempi) = ...
                            internal_trial_count_newSeq{file_index, target_seqLength}(tempi) + 1;                        
                    end                                        
                end
            end                                             
        end
        offloadingProb{file_index, target_seqLength} = offloading_trial_count_newSeq{file_index, target_seqLength}./ ...
            (offloading_trial_count_newSeq{file_index, target_seqLength}+internal_trial_count_newSeq{file_index, target_seqLength});
                
               
    end
    

end

numSeq = zeros(1, pointKindsNum);
for target_seqLength=1:pointKindsNum
    numSeq(target_seqLength) = size(target_seqSet{target_seqLength}, 1);
end

for file_index=1:fileSize
    for target_seqLength=1:pointKindsNum
        current_numSeq = length(target_seqSet{target_seqLength});
        for tempi=1:current_numSeq
        target_seqSetNum{file_index, target_seqLength}(tempi) = length(target_seqSetIndex{file_index, target_seqLength}{tempi});
        end
    end
end


offloading_trial_count_newSeq_merged = offloading_trial_count_newSeq(1:fileMerged_size, :);
internal_trial_count_newSeq_merged = offloading_trial_count_newSeq_merged;
offloadingProb_merged = offloading_trial_count_newSeq_merged;
for fileMerged_index=1:fileMerged_size
    for target_seqLength=1:pointKindsNum
        offloading_trial_count_newSeq_merged{fileMerged_index, target_seqLength} = zeros(length(offloading_trial_count_newSeq_merged{fileMerged_index, target_seqLength}),1);
        internal_trial_count_newSeq_merged{fileMerged_index, target_seqLength} = zeros(length(offloading_trial_count_newSeq_merged{fileMerged_index, target_seqLength}),1);
        for temp_index=1:fileMergeNum
            
            offloading_trial_count_newSeq_merged{fileMerged_index, target_seqLength} = offloading_trial_count_newSeq_merged{fileMerged_index, target_seqLength} ...
                + offloading_trial_count_newSeq{(fileMerged_index-1)*fileMergeNum+temp_index, target_seqLength};
            internal_trial_count_newSeq_merged{fileMerged_index, target_seqLength} = internal_trial_count_newSeq_merged{fileMerged_index, target_seqLength} ...
                + internal_trial_count_newSeq{(fileMerged_index-1)*fileMergeNum+temp_index, target_seqLength};

        end
        
        offloadingProb_merged{fileMerged_index, target_seqLength} = offloading_trial_count_newSeq_merged{fileMerged_index, target_seqLength}./ ...
            (offloading_trial_count_newSeq_merged{fileMerged_index, target_seqLength}+internal_trial_count_newSeq_merged{fileMerged_index, target_seqLength});
    end
end

% delay variability
delay1_file;
delay2_file;
delay12_file;
isCorrect_file;
ifSelectOffloading_file;
choiceCondition_flag_file;

delay1_collapsed = [];
delay2_collapsed = [];
delay12_collapsed = [];
isCorrect_collapsed = [];
ifSelectOffloading_collapsed = [];
choiceCondition_flag_collapsed = [];

% file_index = 3;

for file_index=1:fileSize
    delay1_collapsed = [delay1_collapsed delay1_file{file_index}];
    delay2_collapsed = [delay2_collapsed delay2_file{file_index}];
    delay12_collapsed = [delay12_collapsed delay12_file{file_index}];
    isCorrect_collapsed = [isCorrect_collapsed isCorrect_file{file_index}];
    ifSelectOffloading_collapsed = [ifSelectOffloading_collapsed ifSelectOffloading_file{file_index}];
    choiceCondition_flag_collapsed = [choiceCondition_flag_collapsed choiceCondition_flag_file{file_index}];    
end
a = 1;

temp_delay1_binSize = 35;%15-->25
temp_delay12_binSize = 70;%30-->50

temp_trialNum = length(delay1_collapsed);
[temp_delay1Min,temp_delay1Max] = bounds(delay1_collapsed);

temp_delay1_binBound = linspace(temp_delay1Min,temp_delay1Max,ceil((temp_delay1Max-temp_delay1Min)/temp_delay1_binSize));
temp_delay1_binNum = length(temp_delay1_binBound)-1;
% temp_delay1_binSize = temp_delay1_binBound(2)-temp_delay1_binBound(1);

temp_delay1_binMean = zeros(1, temp_delay1_binNum);
for tempi=1:temp_delay1_binNum
    temp_delay1_binMean(tempi) = mean([temp_delay1_binBound(tempi), temp_delay1_binBound(tempi+1)]);
end


delay1_binCategory = zeros(1, temp_trialNum);

for tempi=1:temp_trialNum
    delay1_binCategory(tempi) = find(delay1_collapsed(tempi)>=temp_delay1_binBound,1,'last');
end
delay1_binCategory(delay1_binCategory>temp_delay1_binNum) = temp_delay1_binNum;

isCorrect_delay1_binCategory = cell(1, temp_delay1_binNum);
ifSelectOffloading_delay1_binCategory = cell(1, temp_delay1_binNum);
choiceCondition_flag_delay1_binCategory = cell(1, temp_delay1_binNum);
for tempi=1:temp_delay1_binNum
    isCorrect_delay1_binCategory{tempi} = isCorrect_collapsed(delay1_binCategory == tempi);    
    ifSelectOffloading_delay1_binCategory{tempi} = ifSelectOffloading_collapsed(delay1_binCategory == tempi);    
    choiceCondition_flag_delay1_binCategory{tempi} = choiceCondition_flag_collapsed(delay1_binCategory == tempi);        
end

gAcc_noChoice_delay1_binCategory = zeros(1, temp_delay1_binNum);
rProb_delay1_binCategory = zeros(1, temp_delay1_binNum);
for tempi=1:temp_delay1_binNum
    tempIndex_noChoice = choiceCondition_flag_delay1_binCategory{tempi} == 0;
    temp_isCorrect_noChoice = isCorrect_delay1_binCategory{tempi}(tempIndex_noChoice);
    gAcc_noChoice_delay1_binCategory(tempi) = sum(temp_isCorrect_noChoice==1)/length(temp_isCorrect_noChoice);
    
    tempIndex_choice = choiceCondition_flag_delay1_binCategory{tempi} == 2;
    rProb_delay1_binCategory(tempi) = sum(ifSelectOffloading_delay1_binCategory{tempi}(tempIndex_choice)==1)/sum(tempIndex_choice);
end




[temp_delay12Min,temp_delay12Max] = bounds(delay12_collapsed);

temp_delay12_binBound = linspace(temp_delay12Min,temp_delay12Max,ceil((temp_delay12Max-temp_delay12Min)/temp_delay12_binSize));
temp_delay12_binNum = length(temp_delay12_binBound)-1;
% temp_delay12_binSize = temp_delay12_binBound(2)-temp_delay12_binBound(1);

temp_delay12_binMean = zeros(1, temp_delay12_binNum);
for tempi=1:temp_delay12_binNum
    temp_delay12_binMean(tempi) = mean([temp_delay12_binBound(tempi), temp_delay12_binBound(tempi+1)]);
end


delay12_binCategory = zeros(1, temp_trialNum);

for tempi=1:temp_trialNum
    delay12_binCategory(tempi) = find(delay12_collapsed(tempi)>=temp_delay12_binBound,1,'last');
end
delay12_binCategory(delay12_binCategory>temp_delay12_binNum) = temp_delay12_binNum;

isCorrect_delay12_binCategory = cell(1, temp_delay12_binNum);
ifSelectOffloading_delay12_binCategory = cell(1, temp_delay12_binNum);
choiceCondition_flag_delay12_binCategory = cell(1, temp_delay12_binNum);
for tempi=1:temp_delay12_binNum
    isCorrect_delay12_binCategory{tempi} = isCorrect_collapsed(delay12_binCategory == tempi);    
    ifSelectOffloading_delay12_binCategory{tempi} = ifSelectOffloading_collapsed(delay12_binCategory == tempi);    
    choiceCondition_flag_delay12_binCategory{tempi} = choiceCondition_flag_collapsed(delay12_binCategory == tempi);        
end

gAcc_noChoice_delay12_binCategory = zeros(1, temp_delay12_binNum);
for tempi=1:temp_delay12_binNum
    tempIndex_noChoice = choiceCondition_flag_delay12_binCategory{tempi} == 0;
    temp_isCorrect_noChoice = isCorrect_delay12_binCategory{tempi}(tempIndex_noChoice);
    gAcc_noChoice_delay12_binCategory(tempi) = sum(temp_isCorrect_noChoice==1)/length(temp_isCorrect_noChoice);    
end



fig111 = figure(111);
set(gcf,'Position',[0 50 1100 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

nexttile
plot(1-gAcc_noChoice_delay1_binCategory,'-o', 'Color', [0.25 0.25 0.25 0.75],'MarkerSize',5);
hold on
plot(rProb_delay1_binCategory,'-o', 'Color', [0.6350 0.0780 0.1840],'MarkerSize',5);

hl = legend('error rate',...
    'offloading rate',...
    'Location','southeast','fontsize',11);

set(hl,'box','off');

ylim([0 1]);
set(gca, 'FontSize', 12)
set(gca,'XLim',[0 temp_delay1_binNum+1]);%X轴的数据显示范围
set(gca,'XTick',1:3:temp_delay1_binNum);%设置要显示坐标刻度的范围
set(gca,'XTickLabel',round(temp_delay1_binMean(1:3:temp_delay1_binNum)));%给坐标加标签
set(gca,'box','off');% 取消右、上边框
xlabel('Delay1 duration (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 14, 'FontWeight', 'bold');


nexttile

plot(1-gAcc_noChoice_delay12_binCategory,'-o', 'Color', [0.25 0.25 0.25 0.75],'MarkerSize',5);

hl = legend('error rate',...
    'Location','southeast','fontsize',11);

set(hl,'box','off');

ylim([0 1]);
set(gca, 'FontSize', 12)
set(gca,'XLim',[0 temp_delay12_binNum+1]);%X轴的数据显示范围
set(gca,'XTick',1:3:temp_delay12_binNum);%设置要显示坐标刻度的范围
set(gca,'XTickLabel',round(temp_delay12_binMean(1:3:temp_delay12_binNum)));%给坐标加标签
set(gca,'box','off');% 取消右、上边框
xlabel('Delay12 duration (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 14, 'FontWeight', 'bold');



%% %%%%%%%%%%%%%%%%%%%% ProbOffload %%%%%%%%%%%%%%%%%%%
fig1 = figure (1);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
% movegui(gcf, 'northwest');

% tiledlayout(4,1,'Padding','Compact');
% tiledlayout(4,1);
%tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact');

offloading_trial_count_newSeq_all = offloadingProb(1, :);
offloadingCorrect_trial_count_newSeq_all = offloadingProb(1, :);
internal_trial_count_newSeq_all = offloadingProb(1, :);
offloadingProb_all = offloadingProb(1, :);
offloadingAcc_all = offloadingProb(1, :);

for target_seqLength=1:pointKindsNum    
    current_numSeq = length(target_seqSet{target_seqLength});
    for tempi=1:current_numSeq
        offloading_trial_count_newSeq_all{target_seqLength} = zeros(current_numSeq, 1);
        offloadingCorrect_trial_count_newSeq_all{target_seqLength} = zeros(current_numSeq, 1);        
        internal_trial_count_newSeq_all{target_seqLength} = zeros(current_numSeq, 1);
    end
end

seqSet_inOne = offloadingProb(1, :);
kClusters = offloadingProb(1, :);
KCentroids = offloadingProb(1, :);
for target_seqLength=1:pointKindsNum
    current_numSeq = length(target_seqSet{target_seqLength});
    
    for file_index=1:fileSize
        offloading_trial_count_newSeq_all{target_seqLength} = ...
            offloading_trial_count_newSeq_all{target_seqLength} + ...
            offloading_trial_count_newSeq{file_index, target_seqLength};
        
        internal_trial_count_newSeq_all{target_seqLength} = ...
            internal_trial_count_newSeq_all{target_seqLength} + ...
            internal_trial_count_newSeq{file_index, target_seqLength};
        
        offloadingCorrect_trial_count_newSeq_all{target_seqLength} = ...
            offloadingCorrect_trial_count_newSeq_all{target_seqLength} + ...
            offloadingCorrect_trial_count_newSeq{file_index, target_seqLength};
        
    end
    
    offloadingProb_all{target_seqLength} = offloading_trial_count_newSeq_all{target_seqLength}./ ...
        (offloading_trial_count_newSeq_all{target_seqLength}+internal_trial_count_newSeq_all{target_seqLength});
    
    offloadingProb_all{target_seqLength}(isnan(offloadingProb_all{target_seqLength})) = 0;
    
    offloadingAcc_all{target_seqLength} = offloadingCorrect_trial_count_newSeq_all{target_seqLength}./ ...
        offloading_trial_count_newSeq_all{target_seqLength};
    
    offloadingAcc_all{target_seqLength}(isnan(offloadingAcc_all{target_seqLength})) = 0;
    
    subplot(pointKindsNum, 1, target_seqLength);
    %nexttile
    %plot(1:current_numSeq, offloadingProb_all{target_seqLength},'--o', 'LineWidth', 2);
%     plot(1:current_numSeq, offloadingProb_all{target_seqLength},'--s', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');    
    plot(1:current_numSeq, offloadingProb_all{target_seqLength},'-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.75 0.25 0.25]);    

    hold on
    
    for fileMerged_index=1:fileMerged_size
%         scatter(1:current_numSeq, offloadingProb_merged{fileMerged_index, target_seqLength}, 25, 'filled', ...
%             'MarkerFaceColor', scatter_rgbColor(fileMerged_index, :), ...
%             'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);  

%         scatter(1:current_numSeq, offloadingProb_merged{fileMerged_index, target_seqLength}, 17, 'filled', ...
%             'MarkerFaceColor', scatter_rgbColor(fileMerged_index, :), ...
%             'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);   
        scatter(1:current_numSeq, offloadingProb_merged{fileMerged_index, target_seqLength}, 10, 'filled', ...
            'MarkerFaceColor', scatter_rgbColor(fileMerged_index, :), ...
            'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);   
        hold on
    end     
    
    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
    
    %XTickLabel
    seqSet_inOne{target_seqLength} = zeros(1, current_numSeq);
    for tempi=1:current_numSeq
        for tempj=1:target_seqLength
            seqSet_inOne{target_seqLength}(tempi) = 10^(target_seqLength-tempj) * ...
                target_seqSet{target_seqLength}{tempi}(tempj) + seqSet_inOne{target_seqLength}(tempi);
        end
    end
    
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength});%给坐标加标签
    ylabel('Offloading rate', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('ProbOffload-SeqLength=%d, All-%ddays', target_seqLength, fileSize), 'FontSize',15);    
    %title(sprintf('ProbOffload-SeqLength=%d, All-%ddays', target_seqLength, fileSize), 'FontSize',12);
    % 获取xticklabel的忿
    xtl=get(gca,'XTickLabel');
    % 获取xtick的忿
    xt=get(gca,'XTick');
    % 获取ytick的忿
    yt=get(gca,'YTick');
    % 设置text的x坐标位置仿
    xtextp=xt;
    % 设置text的y坐标位置仿
    ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.08;
    %ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))*0.4-0.05;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',13,'FontWeight','bold');
    %text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',11,'FontWeight','bold');
    % 取消原始ticklabel
    set(gca,'xticklabel','');    
    % 取消右、上边框
    set(gca,'box','off');
    
    kNum = 4;
    [kClusters{target_seqLength}, KCentroids{target_seqLength}] = kmeans(offloadingProb_all{target_seqLength},kNum,'Replicates',10);
    tempVar = sort(KCentroids{target_seqLength});
    tempSort = tempVar;
    tempUsed = [];
    tempValid = 1;
    for tempk=1:kNum
        for tempi=1:kNum
            if tempVar(tempk) == KCentroids{target_seqLength}(tempi) && tempValid == 1 && ...
                ismember(tempi, tempUsed) == 0
                
                tempSort(tempk) = tempi;
                tempUsed = [tempUsed tempi];
                tempValid = 0;
            end 
        end
        tempValid = 1;
    end
    KCentroids{target_seqLength} = tempVar;
        
    for tempi=1:current_numSeq
        kClusters{target_seqLength}(tempi) = find(tempSort==kClusters{target_seqLength}(tempi));        
    end
end
if ifSaveOffloadingProb == 1
    currentName = ['offloadingProb', '_'];    
    fileName_rProb = getFileNameMonkey_MAT(subject_name, currentName, monkey_name);
    save(fileName_rProb, 'offloadingProb'); 
    save(fileName_rProb, 'offloadingProb_all', '-append');
    
    if false
        if if_monkey_D0_Z1 == 0
            fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\D_offloadingMatrix_2.mat';            
        elseif if_monkey_D0_Z1 == 1
            fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\Z_offloadingMatrix_2.mat';            
        end
        %offloadingMatrix = offloading_trial_count_newSeq_all;
        offloadingMatrix = [];
        for tempi=1:pointKindsNum
            offloadingMatrix = [offloadingMatrix; offloading_trial_count_newSeq_all{tempi}];            
        end
        save(fileName, 'offloadingMatrix');
    end    
end

currentFigName = ['rProb', '_'];
% to generate a unique file name for saving figure
fileName_fig1 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(1) == 1
        exportgraphics(fig1, fileName_fig1, 'BackgroundColor', 'none');
end
%% %%%%%%%%%%%%%%%%%%%% ProbOffload, clustered %%%%%%%%%%%%%%%%%%%
clustered_offloadingProb_all = offloadingProb_all;
clustered_seqSet_inOne = seqSet_inOne;
for target_seqLength=1:pointKindsNum
    current_numSeq = length(target_seqSet{target_seqLength});
    tempi = 1;
    for tempk=1:kNum
        for tempj=1:current_numSeq
            if kClusters{target_seqLength}(tempj) == tempk
                clustered_offloadingProb_all{target_seqLength}(tempi) = ...
                    offloadingProb_all{target_seqLength}(tempj);
                
                clustered_seqSet_inOne{target_seqLength}(tempi) = ...
                    seqSet_inOne{target_seqLength}(tempj);
                
                tempi = tempi + 1;
            end                        
        end
    end    
                
end

%% %%%%%%%%%%%%%%%%%%%% ProbOffload, clustered, scattered %%%%%%%%%%%%%%%%%%%
fig2 = figure(2);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
%tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact');


for target_seqLength=1:pointKindsNum
    current_numSeq = length(target_seqSet{target_seqLength});
    subplot(pointKindsNum, 1, target_seqLength);
    %nexttile
    
    for tempk=1:kNum
        tempVar = sum(sum(kClusters{target_seqLength}==[1:tempk-1]));
        plot(tempVar+1:tempVar+sum(kClusters{target_seqLength}==tempk), offloadingProb_all{target_seqLength}(kClusters{target_seqLength}==tempk) ...
            ,'--s', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
        hold on
    end
    
    
    %offloadingProb_merged       
    for tempk=1:kNum
        tempVar = sum(sum(kClusters{target_seqLength}==[1:tempk-1]));
        for fileMerged_index=1:fileMerged_size
            scatter(tempVar+1:tempVar+sum(kClusters{target_seqLength}==tempk),...
                offloadingProb_merged{fileMerged_index, target_seqLength}(kClusters{target_seqLength}==tempk), 10, 'filled', ...
                'MarkerFaceColor', scatter_rgbColor(fileMerged_index, :), ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
    end


    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
    
    
    set(gca,'XTickLabel',clustered_seqSet_inOne{target_seqLength});%给坐标加标签
    ylabel('Offloading rate', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('ProbOffload-SeqLength=%d, All-%ddays-scattered, Clustered, ', target_seqLength, fileSize), 'FontSize',15);    

    % 获取xticklabel的忿
    xtl=get(gca,'XTickLabel');
    % 获取xtick的忿
    xt=get(gca,'XTick');
    % 获取ytick的忿
    yt=get(gca,'YTick');
    % 设置text的x坐标位置仿
    xtextp=xt;
    % 设置text的y坐标位置仿
    ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.08;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',13,'FontWeight','bold');
    % 取消原始ticklabel
    set(gca,'xticklabel','');    
    % 取消右、上边框
    set(gca,'box','off');    
        
end

currentFigName = ['kmeans', '_'];
% to generate a unique file name for saving figure
fileName_fig2 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(2) == 1
        exportgraphics(fig2, fileName_fig2, 'BackgroundColor', 'none');
end

%% %%%%%%%%%%%%%%%%%%%% ProbOffload, Collapsed in one plot %%%%%%%%%%%%%%%%%%%
fig3 = figure(3);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

current_numSeq = 0;
for target_seqLength=1:pointKindsNum
    current_numSeq = length(target_seqSet{target_seqLength}) + current_numSeq;
end

collapsed_offloadingProb = [];
for target_seqLength=1:pointKindsNum
    collapsed_offloadingProb = [collapsed_offloadingProb; offloadingProb_all{target_seqLength}]; 
end

% plot(1:current_numSeq, collapsed_offloadingProb,'--o', 'LineWidth', 2);


[sorted_collapsed_offloadingProb, sorted_index_collapsed_offloadingProb] = sort(collapsed_offloadingProb);

collapsed_seqSet_inOne = []; 
for target_seqLength=1:pointKindsNum
    collapsed_seqSet_inOne = [collapsed_seqSet_inOne seqSet_inOne{target_seqLength}]; 
end
sorted_collapsed_seqSet_inOne = collapsed_seqSet_inOne(sorted_index_collapsed_offloadingProb);

% plot(1:current_numSeq, sorted_collapsed_offloadingProb,'--o', 'LineWidth', 2);
tempFlag = zeros(1, pointKindsNum);
for tempi=1:current_numSeq
    if sorted_collapsed_seqSet_inOne(tempi) <= 9
        if tempFlag(1) == 0
            tempFlag(1) = 1;
            plot(tempi, sorted_collapsed_offloadingProb(tempi),'--s', 'LineWidth', 2, 'Color', [0.25 0.75 0.75], 'MarkerSize', 15, 'MarkerFaceColor', [0.25 0.75 0.75]);
        end
    elseif sorted_collapsed_seqSet_inOne(tempi) <= 99
        if tempFlag(1) == 1
            tempFlag(1) = 2;
            tempFlag(2) = 1;
            plot(tempi, sorted_collapsed_offloadingProb(tempi),'--s', 'LineWidth', 2, 'Color', [0.75 0.25 0.75], 'MarkerSize', 15, 'MarkerFaceColor', [0.75 0.25 0.75]);
        end
    elseif sorted_collapsed_seqSet_inOne(tempi)  <= 999
        if tempFlag(2) == 1
            tempFlag(1) = 2;
            tempFlag(2) = 2;
            tempFlag(3) = 1;
            plot(tempi, sorted_collapsed_offloadingProb(tempi),'--s', 'LineWidth', 2, 'Color', [0.75 0.75 0.25], 'MarkerSize', 15, 'MarkerFaceColor', [0.75 0.75 0.25]);
        end
    elseif sorted_collapsed_seqSet_inOne(tempi)  <= 9999
        if tempFlag(3) == 1
            tempFlag(1) = 2;
            tempFlag(2) = 2;
            tempFlag(3) = 2;
            tempFlag(4) = 1;
            plot(tempi, sorted_collapsed_offloadingProb(tempi),'--s', 'LineWidth', 2, 'Color', [0.25 0.25 0.25], 'MarkerSize', 15, 'MarkerFaceColor', [0.25 0.25 0.25]);
            break
        end
    end
    hold on
end
for tempi=1:current_numSeq
    if sorted_collapsed_seqSet_inOne(tempi) <= 9
        plot(tempi, sorted_collapsed_offloadingProb(tempi),'--s', 'LineWidth', 2, 'Color', [0.25 0.75 0.75], 'MarkerSize', 15, 'MarkerFaceColor', [0.25 0.75 0.75]);
    elseif sorted_collapsed_seqSet_inOne(tempi) <= 99
        plot(tempi, sorted_collapsed_offloadingProb(tempi),'--s', 'LineWidth', 2, 'Color', [0.75 0.25 0.75], 'MarkerSize', 15, 'MarkerFaceColor', [0.75 0.25 0.75]);
    elseif sorted_collapsed_seqSet_inOne(tempi)  <= 999
        plot(tempi, sorted_collapsed_offloadingProb(tempi),'--s', 'LineWidth', 2, 'Color', [0.75 0.75 0.25], 'MarkerSize', 15, 'MarkerFaceColor', [0.75 0.75 0.25]);
    elseif sorted_collapsed_seqSet_inOne(tempi)  <= 9999
        plot(tempi, sorted_collapsed_offloadingProb(tempi),'--s', 'LineWidth', 2, 'Color', [0.25 0.25 0.25], 'MarkerSize', 15, 'MarkerFaceColor', [0.25 0.25 0.25]);
    end
    hold on
end


for fileMerged_index=1:fileMerged_size
    
    offloadingProb_oneFile = offloadingProb_merged(fileMerged_index, :);
    collapsed_offloadingProb_oneFile =[]; 
    for target_seqLength=1:pointKindsNum
        collapsed_offloadingProb_oneFile = [collapsed_offloadingProb_oneFile; offloadingProb_oneFile{target_seqLength}]; 
    end
    
    scatter(1:current_numSeq, collapsed_offloadingProb_oneFile(sorted_index_collapsed_offloadingProb), 10, 'filled', ...
        'MarkerFaceColor', scatter_rgbColor(fileMerged_index, :), ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);       
    hold on
end

if pointKindsNum == 4
    legend('1 item','2 item','3 item','4 item','Location','southeast','fontsize',20)
elseif pointKindsNum == 3
    legend('1 item','2 item','3 item','Location','southeast','fontsize',20)
end

set(gca, 'FontSize', 20)
ylim([0 1]);
set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
set(gca,'XTickLabel',sorted_collapsed_seqSet_inOne);%给坐标加标签
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');
title(sprintf('ProbOffload, Collapsed in one plot, %ddays', fileSize), 'FontSize',20);

% 获取xticklabel的忿
xtl=get(gca,'XTickLabel');
% 获取xtick的忿
xt=get(gca,'XTick');
% 获取ytick的忿
yt=get(gca,'YTick');
% 设置text的x坐标位置仿
xtextp=xt;
% 设置text的y坐标位置仿
ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt));
% rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
% 朿3个属性忼：left，right，center
text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',11,'FontWeight','bold');
% 取消原始ticklabel
set(gca,'xticklabel','');
% 取消右、上边框
set(gca,'box','off');

currentFigName = ['kmeans_onePlot', '_'];
% to generate a unique file name for saving figure
fileName_fig3 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(3) == 1
        exportgraphics(fig3, fileName_fig3, 'BackgroundColor', 'none');
end


%% %%%%%%%%%%% ProbOffload, Collapsed, small right offset, exclude 8th %%%%%%%%%%%%%%
% figure (fileSize+1+1+1+1+1)
% set(gcf,'Position',get(0,'ScreenSize'));%full screen show
% 
% current_numSeq = 0;
% for target_seqLength=1:pointKindsNum
%     current_numSeq = length(target_seqSet{target_seqLength}) + current_numSeq;
% end
% 
% 
% 
% % temp_indexToExclude = [];
% % temp_seqSetToExclude = [];
% valid_index = [];
% valid_seqSet= [];
% for tempi=1:current_numSeq       
%     temp_unit = collapsed_seqSet_inOne(tempi)-floor(collapsed_seqSet_inOne(tempi)/10)*10;    
%     
%     if tempi <= length(target_seqSet{1})
%         tempj = 1;
%         tempk = tempi;
%     elseif tempi <= length(target_seqSet{1})+length(target_seqSet{2})
%         tempj = 2;
%         tempk = tempi-length(target_seqSet{1});
%     elseif tempi <= length(target_seqSet{1})+length(target_seqSet{2})+length(target_seqSet{3})
%         tempj = 3;
%         tempk = tempi-length(target_seqSet{1})-length(target_seqSet{2});
%     else
%         tempj = 4;
%         tempk = tempi-length(target_seqSet{1})-length(target_seqSet{2})-length(target_seqSet{3});
%     end    
%     rightOffset = sum(target_seqSet{tempj}{tempk})/tempj - 4.5;
%     
% %     if temp_unit == 8 || rightOffset > 1
% % %     if temp_unit == 8
% % %     if rightOffset > 1
% %         temp_indexToExclude = [temp_indexToExclude tempi]; %#ok<AGROW>
% %         temp_seqSetToExclude = [temp_seqSetToExclude collapsed_seqSet_inOne(tempi)]; %#ok<AGROW>       
% %     end
% %     
% % end
%     if temp_unit ~= 8 && abs(rightOffset) <= 1.5
% %     if abs(rightOffset) <= 1.5
%         valid_index = [valid_index tempi]; %#ok<AGROW>
%         valid_seqSet = [valid_seqSet collapsed_seqSet_inOne(tempi)]; %#ok<AGROW>       
%     end
%     
% end
% valid_numSeq = length(valid_seqSet);
% 
% valid_sorted_index_collapsed_offloadingProb = zeros(1, valid_numSeq);
% tempj = 1;
% for tempi=1:current_numSeq
%     if ismember(sorted_index_collapsed_offloadingProb(tempi), valid_index)
%         valid_sorted_index_collapsed_offloadingProb(tempj) = sorted_index_collapsed_offloadingProb(tempi);
%         tempj = tempj + 1;        
%     end  
% end
% 
% valid_sorted_collapsed_seqSet_inOne = collapsed_seqSet_inOne(valid_sorted_index_collapsed_offloadingProb);
% valid_sorted_collapsed_offloadingProb = collapsed_offloadingProb(valid_sorted_index_collapsed_offloadingProb);
% 
% for tempi=1:valid_numSeq
%     if valid_sorted_collapsed_seqSet_inOne(tempi) <= 9
%         plot(tempi, valid_sorted_collapsed_offloadingProb(tempi),'--o', 'LineWidth', 2, 'Color', [0 1 1]);
%     elseif valid_sorted_collapsed_seqSet_inOne(tempi) <= 99
%         plot(tempi, valid_sorted_collapsed_offloadingProb(tempi),'--o', 'LineWidth', 2, 'Color', [1 0 1]);
%     elseif valid_sorted_collapsed_seqSet_inOne(tempi)  <= 999
%         plot(tempi, valid_sorted_collapsed_offloadingProb(tempi),'--o', 'LineWidth', 2, 'Color', [0.9 0.9 0]);
%     elseif valid_sorted_collapsed_seqSet_inOne(tempi)  <= 9999
%         plot(tempi, valid_sorted_collapsed_offloadingProb(tempi),'--o', 'LineWidth', 2, 'Color', [0 0 0]);
%     end
%     hold on 
% end
% 
% for file_index=1:fileSize   
%     offloadingProb_oneFile = offloadingProb(file_index, :);
%     collapsed_offloadingProb_oneFile =[]; 
%     for target_seqLength=1:pointKindsNum
%         collapsed_offloadingProb_oneFile = [collapsed_offloadingProb_oneFile; offloadingProb_oneFile{target_seqLength}]; %#ok<AGROW>
%     end
%     
%     scatter(1:valid_numSeq, collapsed_offloadingProb_oneFile(valid_sorted_index_collapsed_offloadingProb), 'x');
%     hold on
% end
% 
% 
% ylim([0 1]);
% set(gca,'XLim',[0 valid_numSeq+1]);%X轴的数据显示范围
% set(gca,'XTick',[1:1:valid_numSeq]);%设置要显示坐标刻度的范围
% 
% 
% set(gca,'XTickLabel',valid_sorted_collapsed_seqSet_inOne);%给坐标加标签
% title(sprintf('ProbOffload, Collapsed in one plot, right offset <= 1.5, exclude 8th, %ddays', fileSize), 'FontSize',15);
% 
% % 获取xticklabel的忿
% xtl=get(gca,'XTickLabel');
% % 获取xtick的忿
% xt=get(gca,'XTick');
% % 获取ytick的忿
% yt=get(gca,'YTick');
% % 设置text的x坐标位置仿
% xtextp=xt;
% % 设置text的y坐标位置仿
% ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt));
% % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
% % 朿3个属性忼：left，right，center
% text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',10,'FontWeight','bold');
% % 取消原始ticklabel
% set(gca,'xticklabel','');
% 

%% %%%%%%%%%%%%%%%%%%%%%%%%% Accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig4 = figure(4);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
% tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact');

gCorrect_trial_count_collapsed = offloadingProb(1, :);
g_trial_count_collapsed = offloadingProb(1, :);
gAcc_collapsed = offloadingProb(1, :);


gCorrect_choice_trial_count;
g_choice_trial_count;
gAcc_choice;
gCorrect_noChoice_trial_count;
g_noChoice_trial_count;
gAcc_noChoice;

gCorrect_choice_trial_count_collapsed = offloadingProb(1, :);
g_choice_trial_count_collapsed = offloadingProb(1, :);
gAcc_choice_collapsed = offloadingProb(1, :);

gCorrect_noChoice_trial_count_collapsed = offloadingProb(1, :);
g_noChoice_trial_count_collapsed = offloadingProb(1, :);
gAcc_noChoice_collapsed = offloadingProb(1, :);

gDistri_noChoice_trial_count_collapsed = cell(1, pointKindsNum);
gDistri_choice_trial_count_collapsed = cell(1, pointKindsNum);
gDistri_length_noChoice_trial_count_collapsed = cell(1, pointKindsNum);
gDistri_length_choice_trial_count_collapsed = cell(1, pointKindsNum);
gDistri_lengthDetail_noChoice_trial_count_collapsed = cell(1, pointKindsNum);
gDistri_lengthDetail_choice_g_trial_count_collapsed = cell(1, pointKindsNum);

for file_index=1:fileSize
    temp_load = eval(['MAT_file_load.file',sprintf('%d',file_index)]); 
    
    gCorrect_trial_count(file_index, :) = target_seqSet;
    g_trial_count(file_index, :) = target_seqSet;
    gAcc(file_index, :) = target_seqSet; 
    
    gCorrect_choice_trial_count(file_index, :) = target_seqSet;
    g_choice_trial_count(file_index, :) = target_seqSet;
    gAcc_choice(file_index, :) = target_seqSet; 
    gCorrect_noChoice_trial_count(file_index, :) = target_seqSet;
    g_noChoice_trial_count(file_index, :) = target_seqSet;
    gAcc_noChoice(file_index, :) = target_seqSet; 
    
    for target_seqLength=1:pointKindsNum
        current_numSeq = length(target_seqSet{target_seqLength});  
        gCorrect_trial_count{file_index, target_seqLength} = zeros(current_numSeq, 1);
        g_trial_count{file_index, target_seqLength} = zeros(current_numSeq, 1);        
        gAcc{file_index, target_seqLength} = zeros(current_numSeq, 1);
        
        gCorrect_choice_trial_count{file_index, target_seqLength} = zeros(current_numSeq, 1);
        g_choice_trial_count{file_index, target_seqLength} = zeros(current_numSeq, 1);
        gAcc_choice{file_index, target_seqLength} = zeros(current_numSeq, 1);
        
        gCorrect_noChoice_trial_count{file_index, target_seqLength} = zeros(current_numSeq, 1);
        g_noChoice_trial_count{file_index, target_seqLength} = zeros(current_numSeq, 1);
        gAcc_noChoice{file_index, target_seqLength} = zeros(current_numSeq, 1);
        
        gCorrect_noChoice_trial_count_inOrder{file_index, target_seqLength} = zeros(target_seqLength, 1);
        g_noChoice_trial_count_inOrder{file_index, target_seqLength} = zeros(target_seqLength, 1);        
        gAcc_noChoice_inOrder{file_index, target_seqLength} = zeros(target_seqLength, 1);
        
        gDistri_noChoice_trial_count{file_index, target_seqLength} = cell(current_numSeq, 1);
        gDistri_choice_trial_count{file_index, target_seqLength} = cell(current_numSeq, 1);
        gDistri_length_noChoice_trial_count{file_index, target_seqLength} = cell(current_numSeq, 1);
        gDistri_length_choice_trial_count{file_index, target_seqLength} = cell(current_numSeq, 1);
        gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength} = cell(current_numSeq, 1);
        gDistri_lengthDetail_choice_g_trial_count{file_index, target_seqLength} = cell(current_numSeq, 1);
        
        for tempi=1:current_numSeq
            gDistri_noChoice_trial_count{file_index, target_seqLength}{tempi} = zeros(current_numSeq, 1);
            gDistri_choice_trial_count{file_index, target_seqLength}{tempi} = zeros(current_numSeq, 1);
            gDistri_length_noChoice_trial_count{file_index, target_seqLength}{tempi} = zeros(numFrames+1, 1);
            gDistri_length_choice_trial_count{file_index, target_seqLength}{tempi} = zeros(numFrames+1, 1);
            gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength}{tempi} = cell(numFrames+1, 1);
            gDistri_lengthDetail_choice_g_trial_count{file_index, target_seqLength}{tempi} = cell(numFrames+1, 1);
            
            % To initialize gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength}{tempi}
            numSeq_extend = zeros(1, numFrames-1);
            numSeq_extend(1:pointKindsNum) = numSeq;
            for tempj=1:length(numSeq_extend)
                if numSeq_extend(tempj) == 0
                   numSeq_extend(tempj) = numSeq_extend(numFrames - tempj); 
                end
            end
            for tempj=1:(numFrames+1)
                if tempj==1 || tempj == (numFrames+1)
                    temp1 = zeros(1, 1);
                else
                    temp1 = zeros(numSeq_extend(tempj-1), 1);
                end
                %gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength}{tempi, tempj} = temp1;
                gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength}{tempi}{tempj} = temp1;
                gDistri_lengthDetail_choice_g_trial_count{file_index, target_seqLength}{tempi}{tempj} = temp1;                
            end
            target_seqSet;
            target_seqSet_extend = cell(1, numFrames-1);
            target_seqSet_extend(1:pointKindsNum) = target_seqSet;
            for tempj=1:length(numSeq_extend)
                if isempty(target_seqSet_extend{tempj})
                   %target_seqSet_extend{numFrames - tempj}; 
                   %numSeq_extend(tempj);
                   for tempk=1:numSeq_extend(tempj)
                       %temp1 = target_seqSet_extend{numFrames - tempj}{numFrames-tempk+1};
                       temp1 = target_seqSet_extend{numFrames - tempj}{numSeq(numFrames - tempj)-tempk+1};
                       temp2 = 1:numFrames;
                       target_seqSet_extend{tempj}{tempk, 1} = temp2(~ismember(temp2, temp1));                       
                   end
                end
            end            
            
            for tempj=1:length(target_seqSetIndex{file_index, target_seqLength}{tempi})                
                trial_count = target_seqSetIndex{file_index, target_seqLength}{tempi}(tempj);
                %if temp_load.trial_para.taskMode_noChoice_Green0_Red1(trial_count) == 0   
                if temp_load.trial_para.ifSelectOffloading(trial_count) == -1                        
                    g_trial_count{file_index, target_seqLength}(tempi) = ...
                        g_trial_count{file_index, target_seqLength}(tempi) + 1;
                    if temp_load.trial_para.isCorrect(trial_count) == 1
                        gCorrect_trial_count{file_index, target_seqLength}(tempi) = ...
                            gCorrect_trial_count{file_index, target_seqLength}(tempi) + 1;
                    end
                end
                
                
                SL = temp_load.basic_para.cumulativeCorrectLimit_seqLengthSwitch;
                if temp_load.trial_para.selectFlag_newSeq(trial_count) == 1 && ...
                        temp_load.trial_para.choiceCondition_flag(trial_count) == 2 &&...
                        floor((trial_count-1)/SL)*SL == trial_count-1                    
                    
                    if temp_load.trial_para.ifSelectOffloading(trial_count) == -1
                        
                        g_choice_trial_count{file_index, target_seqLength}(tempi) = ...
                            g_choice_trial_count{file_index, target_seqLength}(tempi) + 1;
                        
                        if temp_load.trial_para.isCorrect(trial_count) == 1
                            gCorrect_choice_trial_count{file_index, target_seqLength}(tempi) = ...
                                gCorrect_choice_trial_count{file_index, target_seqLength}(tempi) + 1;
                        end
                        
                        current_seq = target_seqSet{target_seqLength}{tempi};
                        current_isFilled = temp_load.trial_para.isFilled{trial_count};
                        current_trial_response = find(current_isFilled==0);
                        
                        current_seq;
                        current_trial_response;
                        gDistri_choice_trial_count{file_index, target_seqLength}{tempi}; %#ok<*VUNUS>
                        if length(current_trial_response) == target_seqLength
                            current_trial_response;
                            temp_response_index = 0;
                            for tempk=1:numSeq(target_seqLength)
                                temp_seq = target_seqSet{target_seqLength}{tempk};
                                if sum(current_trial_response == temp_seq) == target_seqLength
                                    temp_response_index = tempk;
                                end
                            end
                            gDistri_choice_trial_count{file_index, target_seqLength}{tempi}(temp_response_index) = ...
                                gDistri_choice_trial_count{file_index, target_seqLength}{tempi}(temp_response_index) + 1;
                        end
                        current_trial_response;
                        gDistri_length_choice_trial_count{file_index, target_seqLength}{tempi}(length(current_trial_response)+1) = ...
                            gDistri_length_choice_trial_count{file_index, target_seqLength}{tempi}(length(current_trial_response)+1) + 1;
                        
                        
                        gDistri_lengthDetail_choice_g_trial_count{file_index, target_seqLength}{tempi};
                        current_trial_response_index = 1;
                        length_response = length(current_trial_response);
                        if length_response >= 1 && length_response <= (numFrames-1)
                            for tempk = 1:numSeq_extend(length_response)
                                temp_seq = target_seqSet_extend{length_response}{tempk};
                                if sum(current_trial_response == temp_seq) == length_response
                                    current_trial_response_index = tempk;
                                end
                            end
                        end
                        current_trial_response_index;
                        gDistri_lengthDetail_choice_g_trial_count{file_index, target_seqLength}{tempi}{length_response+1}(current_trial_response_index) = ...
                            gDistri_lengthDetail_choice_g_trial_count{file_index, target_seqLength}{tempi}{length_response+1}(current_trial_response_index) + 1;
                        
                        
                    end
                    
                elseif temp_load.trial_para.selectFlag_newSeq(trial_count) == 1 && ...
                        temp_load.trial_para.choiceCondition_flag(trial_count) == 0 &&...
                        floor((trial_count-1)/SL)*SL == trial_count-1
                    
                    g_noChoice_trial_count{file_index, target_seqLength}(tempi) = ...
                        g_noChoice_trial_count{file_index, target_seqLength}(tempi) + 1;
                    
                    if temp_load.trial_para.isCorrect(trial_count) == 1
                        gCorrect_noChoice_trial_count{file_index, target_seqLength}(tempi) = ...
                            gCorrect_noChoice_trial_count{file_index, target_seqLength}(tempi) + 1;
                    end
                        
                    % g_noChoice_trial_count_inOrder{file_index, target_seqLength};
                    current_seq = target_seqSet{target_seqLength}{tempi};
                    current_isFilled = temp_load.trial_para.isFilled{trial_count};
                    current_trial_response = find(current_isFilled==0);
                    
                    current_seq;
                    current_trial_response;
                    for tempk=1:target_seqLength
                        g_noChoice_trial_count_inOrder{file_index, target_seqLength}(tempk) = ...
                            g_noChoice_trial_count_inOrder{file_index, target_seqLength}(tempk) + 1;
                        if ismember(current_seq(tempk), current_trial_response)
                            gCorrect_noChoice_trial_count_inOrder{file_index, target_seqLength}(tempk) = ...
                                gCorrect_noChoice_trial_count_inOrder{file_index, target_seqLength}(tempk) + 1;
                        end
                    end
                    
                    % gDistri_noChoice_trial_count
                    current_seq;
                    current_trial_response;
                    gDistri_noChoice_trial_count{file_index, target_seqLength}{tempi};
                    if length(current_trial_response) == target_seqLength                        
                        current_trial_response;
                        temp_response_index = 0;
                        for tempk=1:numSeq(target_seqLength)
                            temp_seq = target_seqSet{target_seqLength}{tempk};
                            if sum(current_trial_response == temp_seq) == target_seqLength
                                temp_response_index = tempk;
                            end
                        end
                        gDistri_noChoice_trial_count{file_index, target_seqLength}{tempi}(temp_response_index) = ...
                            gDistri_noChoice_trial_count{file_index, target_seqLength}{tempi}(temp_response_index) + 1;
                    end
                    current_trial_response;
                    gDistri_length_noChoice_trial_count{file_index, target_seqLength}{tempi}(length(current_trial_response)+1) = ...
                        gDistri_length_noChoice_trial_count{file_index, target_seqLength}{tempi}(length(current_trial_response)+1) + 1;
                    
                    
                    gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength}{tempi};
                    current_trial_response_index = 1;
                    length_response = length(current_trial_response);
                    if length_response >= 1 && length_response <= (numFrames-1)
                        for tempk = 1:numSeq_extend(length_response)
                            temp_seq = target_seqSet_extend{length_response}{tempk};
                            if sum(current_trial_response == temp_seq) == length_response
                                current_trial_response_index = tempk;
                            end
                        end
                    end
                    current_trial_response_index;
                    gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength}{tempi}{length_response+1}(current_trial_response_index) = ...
                        gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength}{tempi}{length_response+1}(current_trial_response_index) + 1;
                    
                    a = 1;
                end
            end
        end
        
        gAcc{file_index, target_seqLength} = gCorrect_trial_count{file_index, target_seqLength}./...
            g_trial_count{file_index, target_seqLength};
        
        gAcc_choice{file_index, target_seqLength} = gCorrect_choice_trial_count{file_index, target_seqLength}./...
            g_choice_trial_count{file_index, target_seqLength};
        
        gAcc_noChoice{file_index, target_seqLength} = gCorrect_noChoice_trial_count{file_index, target_seqLength}./...
            g_noChoice_trial_count{file_index, target_seqLength};
        
        %gAcc_noChoice_inOrder{file_index, target_seqLength};
        for tempk=1:target_seqLength
            gAcc_noChoice_inOrder{file_index, target_seqLength}(tempk) = ...
                gCorrect_noChoice_trial_count_inOrder{file_index, target_seqLength}(tempk) ./...
                g_noChoice_trial_count_inOrder{file_index, target_seqLength}(tempk);
        end
        
        if file_index == 1
            gCorrect_trial_count_collapsed{target_seqLength} = gCorrect_trial_count{file_index, target_seqLength};
            g_trial_count_collapsed{target_seqLength} = g_trial_count{file_index, target_seqLength};
            
            gCorrect_choice_trial_count_collapsed{target_seqLength} = gCorrect_choice_trial_count{file_index, target_seqLength};
            g_choice_trial_count_collapsed{target_seqLength} = g_choice_trial_count{file_index, target_seqLength};            
            
            gCorrect_noChoice_trial_count_collapsed{target_seqLength} = gCorrect_noChoice_trial_count{file_index, target_seqLength};
            g_noChoice_trial_count_collapsed{target_seqLength} = g_noChoice_trial_count{file_index, target_seqLength};            
            
            gDistri_noChoice_trial_count_collapsed{target_seqLength} = gDistri_noChoice_trial_count{file_index, target_seqLength};
            gDistri_choice_trial_count_collapsed{target_seqLength} = gDistri_choice_trial_count{file_index, target_seqLength};
            gDistri_length_noChoice_trial_count_collapsed{target_seqLength} = gDistri_length_noChoice_trial_count{file_index, target_seqLength};
            gDistri_length_choice_trial_count_collapsed{target_seqLength} = gDistri_length_choice_trial_count{file_index, target_seqLength};
            gDistri_lengthDetail_noChoice_trial_count_collapsed{target_seqLength} = gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength};
            gDistri_lengthDetail_choice_g_trial_count_collapsed{target_seqLength} = gDistri_lengthDetail_choice_g_trial_count{file_index, target_seqLength};
            
        else
            gCorrect_trial_count_collapsed{target_seqLength} = gCorrect_trial_count_collapsed{target_seqLength}...
                + gCorrect_trial_count{file_index, target_seqLength};
            g_trial_count_collapsed{target_seqLength} = g_trial_count_collapsed{target_seqLength}...
                + g_trial_count{file_index, target_seqLength};
            
            gCorrect_choice_trial_count_collapsed{target_seqLength} = gCorrect_choice_trial_count_collapsed{target_seqLength}...
                + gCorrect_choice_trial_count{file_index, target_seqLength};
            g_choice_trial_count_collapsed{target_seqLength} = g_choice_trial_count_collapsed{target_seqLength}...
                + g_choice_trial_count{file_index, target_seqLength};
            
            gCorrect_noChoice_trial_count_collapsed{target_seqLength} = gCorrect_noChoice_trial_count_collapsed{target_seqLength}...
                + gCorrect_noChoice_trial_count{file_index, target_seqLength};
            g_noChoice_trial_count_collapsed{target_seqLength} = g_noChoice_trial_count_collapsed{target_seqLength}...
                + g_noChoice_trial_count{file_index, target_seqLength};
            
            for tempi=1:numSeq(target_seqLength)
                gDistri_noChoice_trial_count_collapsed{target_seqLength}{tempi} = ...
                    gDistri_noChoice_trial_count_collapsed{target_seqLength}{tempi} + ...
                    gDistri_noChoice_trial_count{file_index, target_seqLength}{tempi};   
                gDistri_choice_trial_count_collapsed{target_seqLength}{tempi} = ...
                    gDistri_choice_trial_count_collapsed{target_seqLength}{tempi} + ...
                    gDistri_choice_trial_count{file_index, target_seqLength}{tempi};    
                gDistri_length_noChoice_trial_count_collapsed{target_seqLength}{tempi} = ...
                    gDistri_length_noChoice_trial_count_collapsed{target_seqLength}{tempi} + ...
                    gDistri_length_noChoice_trial_count{file_index, target_seqLength}{tempi};     
                gDistri_length_choice_trial_count_collapsed{target_seqLength}{tempi} = ...
                    gDistri_length_choice_trial_count_collapsed{target_seqLength}{tempi} + ...
                    gDistri_length_choice_trial_count{file_index, target_seqLength}{tempi};                
                
                gDistri_lengthDetail_noChoice_trial_count_collapsed{target_seqLength}{tempi};
                gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength}{tempi};
                for tempj=1:(numFrames+1)
                    gDistri_lengthDetail_noChoice_trial_count_collapsed{target_seqLength}{tempi}{tempj} = ...
                        gDistri_lengthDetail_noChoice_trial_count_collapsed{target_seqLength}{tempi}{tempj} + ...
                        gDistri_lengthDetail_noChoice_trial_count{file_index, target_seqLength}{tempi}{tempj};
                end
                gDistri_lengthDetail_choice_g_trial_count_collapsed{target_seqLength}{tempi};
                gDistri_lengthDetail_choice_g_trial_count{file_index, target_seqLength}{tempi};
                for tempj=1:(numFrames+1)
                    gDistri_lengthDetail_choice_g_trial_count_collapsed{target_seqLength}{tempi}{tempj} = ...
                        gDistri_lengthDetail_choice_g_trial_count_collapsed{target_seqLength}{tempi}{tempj} + ...
                        gDistri_lengthDetail_choice_g_trial_count{file_index, target_seqLength}{tempi}{tempj};
                end                
            end
        end
        
    end
end
g_noChoice_trial_count_collapsed;
gDistri_noChoice_trial_count_collapsed;
gDistri_noChoice_collapsed = gDistri_noChoice_trial_count_collapsed;

g_choice_trial_count_collapsed;
gDistri_choice_trial_count_collapsed;
gDistri_choice_collapsed = gDistri_choice_trial_count_collapsed;

gDistri_length_noChoice_trial_count_collapsed;
gDistri_length_noChoice_collapsed = gDistri_length_noChoice_trial_count_collapsed;
gDistri_length_choice_trial_count_collapsed;
gDistri_length_choice_collapsed = gDistri_length_choice_trial_count_collapsed;

gDistri_lengthDetail_noChoice_trial_count_collapsed;
gDistri_lengthDetail_noChoice_collapsed = gDistri_lengthDetail_noChoice_trial_count_collapsed;
gDistri_lengthDetail_choice_g_trial_count_collapsed;
gDistri_lengthDetail_choice_g_collapsed = gDistri_lengthDetail_choice_g_trial_count_collapsed;
for target_seqLength=1:pointKindsNum
    gDistri_noChoice_trial_count_collapsed{target_seqLength};
    for tempi=1:numSeq(target_seqLength)
        temp_distri1 = gDistri_noChoice_trial_count_collapsed{target_seqLength}{tempi};
        gDistri_noChoice_collapsed{target_seqLength}{tempi} = ...
            temp_distri1 / g_noChoice_trial_count_collapsed{target_seqLength}(tempi);    
        
        temp_distri2 = gDistri_choice_trial_count_collapsed{target_seqLength}{tempi};
        gDistri_choice_collapsed{target_seqLength}{tempi} = ...
            temp_distri2 / g_choice_trial_count_collapsed{target_seqLength}(tempi);  
        
        temp_distri3 = gDistri_length_noChoice_trial_count_collapsed{target_seqLength}{tempi};
        gDistri_length_noChoice_collapsed{target_seqLength}{tempi} = ...
            temp_distri3 / g_noChoice_trial_count_collapsed{target_seqLength}(tempi); 
        
        temp_distri4 = gDistri_length_choice_trial_count_collapsed{target_seqLength}{tempi};
        gDistri_length_choice_collapsed{target_seqLength}{tempi} = ...
            temp_distri4 / g_choice_trial_count_collapsed{target_seqLength}(tempi);   
        
        gDistri_lengthDetail_noChoice_trial_count_collapsed{target_seqLength}{tempi};
        for tempj=1:(numFrames+1)
            gDistri_lengthDetail_noChoice_trial_count_collapsed{target_seqLength}{tempi}{tempj};
            
            gDistri_lengthDetail_noChoice_collapsed{target_seqLength}{tempi}{tempj} = ...
                gDistri_lengthDetail_noChoice_trial_count_collapsed{target_seqLength}{tempi}{tempj} ...
                / g_noChoice_trial_count_collapsed{target_seqLength}(tempi);
        end    
        gDistri_lengthDetail_choice_g_trial_count_collapsed{target_seqLength}{tempi};
        for tempj=1:(numFrames+1)
            gDistri_lengthDetail_choice_g_trial_count_collapsed{target_seqLength}{tempi}{tempj};
            
            gDistri_lengthDetail_choice_g_collapsed{target_seqLength}{tempi}{tempj} = ...
                gDistri_lengthDetail_choice_g_trial_count_collapsed{target_seqLength}{tempi}{tempj} ...
                / g_choice_trial_count_collapsed{target_seqLength}(tempi);
        end          
    end
end

target_seqSet_extend;
target_seqSet_extend_inOne = [];
for tempi=1:size(target_seqSet_extend, 2)
    target_seqSet_extend_inOne = [target_seqSet_extend_inOne; target_seqSet_extend{tempi}];    
end

gDistri_lengthDetail_noChoice_trial_count_collapsed;
gDistri_lengthDetail_noChoice_trial_count_inOne_raw = [];
gDistri_lengthDetail_choice_g_trial_count_collapsed;
gDistri_lengthDetail_choice_g_trial_count_inOne_raw = [];
for tempi=1:pointKindsNum
    gDistri_lengthDetail_noChoice_trial_count_inOne_raw = [gDistri_lengthDetail_noChoice_trial_count_inOne_raw; ...
        gDistri_lengthDetail_noChoice_trial_count_collapsed{tempi}];
    gDistri_lengthDetail_choice_g_trial_count_inOne_raw = [gDistri_lengthDetail_choice_g_trial_count_inOne_raw; ...
        gDistri_lengthDetail_choice_g_trial_count_collapsed{tempi}];    
end
gDistri_lengthDetail_noChoice_trial_count_inOne = cell(sum(numSeq), 1);
gDistri_lengthDetail_choice_g_trial_count_inOne = cell(sum(numSeq), 1);

for tempi=1:sum(numSeq)
    %gDistri_lengthDetail_noChoice_trial_count_inOne{tempi};
    %gDistri_lengthDetail_noChoice_trial_count_inOne_raw{tempi};
    for tempj=2:size(gDistri_lengthDetail_noChoice_trial_count_inOne_raw{tempi}, 1)-1    
        gDistri_lengthDetail_noChoice_trial_count_inOne_raw{tempi}{tempj};
        gDistri_lengthDetail_noChoice_trial_count_inOne{tempi} = ...
            [gDistri_lengthDetail_noChoice_trial_count_inOne{tempi}; ...
            gDistri_lengthDetail_noChoice_trial_count_inOne_raw{tempi}{tempj};];
    end
    for tempj=2:size(gDistri_lengthDetail_choice_g_trial_count_inOne_raw{tempi}, 1)-1    
        gDistri_lengthDetail_choice_g_trial_count_inOne_raw{tempi}{tempj};
        gDistri_lengthDetail_choice_g_trial_count_inOne{tempi} = ...
            [gDistri_lengthDetail_choice_g_trial_count_inOne{tempi}; ...
            gDistri_lengthDetail_choice_g_trial_count_inOne_raw{tempi}{tempj};];
    end    
end
gDistri_lengthDetail_noChoice_trial_count_inOne_array = zeros(sum(numSeq), ...
    length(gDistri_lengthDetail_noChoice_trial_count_inOne{1}));
gDistri_lengthDetail_choice_g_trial_count_inOne_array = zeros(sum(numSeq), ...
    length(gDistri_lengthDetail_choice_g_trial_count_inOne{1}));
for tempi=1:sum(numSeq)
    gDistri_lengthDetail_noChoice_trial_count_inOne_array(tempi, :) = ...
        gDistri_lengthDetail_noChoice_trial_count_inOne{tempi};
    gDistri_lengthDetail_choice_g_trial_count_inOne_array(tempi, :) = ...
        gDistri_lengthDetail_choice_g_trial_count_inOne{tempi};    
end
gDistri_lengthDetail_noChoice_trial_count_inOne_array;
gDistri_lengthDetail_choice_g_trial_count_inOne_array;
% for tempi=1:sum(numSeq)
%     gDistri_lengthDetail_noChoice_trial_count_inOne{tempi};
%     current_seq = target_seqSet_extend_inOne{tempi};
%     current_length = length(current_seq);
%     gDistri_lengthDetail_noChoice_trial_count_collapsed{}
%     
%     
% end


% It's the response distribution of specific seqs.(Althouth rule out the
% error situation that seq_length is wrong.) 
gDistri_noChoice_collapsed; 
gDistri_choice_collapsed;

gDistri_length_noChoice_collapsed;
gDistri_length_choice_collapsed;

gDistri_lengthDetail_noChoice_collapsed;
gDistri_lengthDetail_choice_g_collapsed;

% gDistri_choice_trial_count;
% gDistri_choice_trial_count_collapsed;
% gDistri_choice_collapsed = gDistri_choice_trial_count_collapsed;
% for target_seqLength=1:pointKindsNum
%     gDistri_choice_trial_count_collapsed{target_seqLength};
%     for tempi=1:numSeq(target_seqLength)
%         temp_distri = gDistri_choice_trial_count_collapsed{target_seqLength}{tempi};
%         gDistri_choice_collapsed{target_seqLength}{tempi} = ...
%             temp_distri / g_choice_trial_count_collapsed{target_seqLength}(tempi);     
%     end
% end
% gDistri_choice_collapsed;

gDistri_length_noChoice_collapsed;
gDistri_length_noChoice_collapsed_seqAVG = cell(1, pointKindsNum);
gDistri_length_choice_collapsed;
gDistri_length_choice_collapsed_seqAVG = cell(1, pointKindsNum);

gDistri_length_noChoice_collapsed_std = cell(1, pointKindsNum);
for tempi=1:pointKindsNum
    gDistri_length_noChoice_collapsed_seqAVG{tempi} = zeros(1, numFrames+1);
    temp_append1 = [];
    temp_append2 = [];
    for tempj=1:numSeq(tempi)
        gDistri_length_noChoice_collapsed{tempi}{tempj};
        temp_append1 = [temp_append1; gDistri_length_noChoice_collapsed{tempi}{tempj}'];
        temp_append2 = [temp_append2; gDistri_length_choice_collapsed{tempi}{tempj}'];        
    end
    isvalid_1 = ~isnan(temp_append1(:, 1));
    isvalid_2 = ~isnan(temp_append2(:, 1));
    gDistri_length_noChoice_collapsed_seqAVG{tempi} = mean(temp_append1(isvalid_1, :), 1);
    gDistri_length_choice_collapsed_seqAVG{tempi} = mean(temp_append2(isvalid_2, :), 1);
    
    gDistri_length_noChoice_collapsed_std{tempi} = std(temp_append1, 1);
end
gDistri_length_noChoice_collapsed_seqAVG;
gDistri_length_noChoice_collapsed_std;

gDistri_lengthDetail_noChoice_collapsed;
gDistri_length_noChoice_collapsed;
gDistri_neighborLengthError_ismember_noChoice_collapsed = cell(1, pointKindsNum);
gDistri_neighborLengthError_ismember_raw_noChoice_collapsed = cell(1, pointKindsNum);
target_seqSet_extend;
for target_seqLength=1:pointKindsNum
    gDistri_neighborLengthError_ismember_noChoice_collapsed{target_seqLength} = ...
        zeros(numSeq(target_seqLength), 3); 
    gDistri_neighborLengthError_ismember_raw_noChoice_collapsed{target_seqLength} = ...
        zeros(numSeq(target_seqLength), 1);     
    gDistri_lengthDetail_noChoice_collapsed{target_seqLength};
    for tempi=1:numSeq(target_seqLength)
        gDistri_lengthDetail_noChoice_collapsed{target_seqLength}{tempi};
        temp_seq = target_seqSet_extend{target_seqLength}{tempi};

        %temp_neighborLength = target_seqLength-1:target_seqLength+1;
        temp_neighborLength_down = target_seqLength-1;
        temp_neighborLength_up = target_seqLength+1;
        
        temp_plus_down = 0;        
        temp_neighborLength_down;
        gDistri_lengthDetail_noChoice_collapsed{target_seqLength}{tempi}{temp_neighborLength_down+1};
        
        if temp_neighborLength_down == 0
%             temp_plus_down = temp_plus_down + ...
%                 gDistri_lengthDetail_noChoice_collapsed{target_seqLength}{tempi}{temp_neighborLength_down+1}(1);            
        else
            for tempk=1:numSeq_extend(temp_neighborLength_down)
                temp_seq_2 = target_seqSet_extend{temp_neighborLength_down}{tempk};
                
                if sum(ismember(temp_seq, temp_seq_2)) == temp_neighborLength_down
                    temp_plus_down = temp_plus_down + ...
                        gDistri_lengthDetail_noChoice_collapsed{target_seqLength}{tempi}{temp_neighborLength_down+1}(tempk);                    
                end
            end            
        end

        temp_plus_up = 0;
        temp_neighborLength_up;
        gDistri_lengthDetail_noChoice_collapsed{target_seqLength}{tempi}{temp_neighborLength_up+1};
        if temp_neighborLength_up == (pointKindsNum+1) %5
            
        else
            for tempk=1:numSeq_extend(temp_neighborLength_up)
                temp_seq_2 = target_seqSet_extend{temp_neighborLength_up}{tempk};
                
                if sum(ismember(temp_seq, temp_seq_2)) == length(temp_seq)
                    temp_plus_up = temp_plus_up + ...
                        gDistri_lengthDetail_noChoice_collapsed{target_seqLength}{tempi}{temp_neighborLength_up+1}(tempk);
                end
            end
        end
        
        temp_plus_down;
        temp_plus_up;
        %gDistri_lengthDetail_noChoice_collapsed{target_seqLength}{tempi};
        gDistri_length_noChoice_collapsed{target_seqLength}{tempi}(temp_neighborLength_down+1);
        gDistri_length_noChoice_collapsed{target_seqLength}{tempi}(temp_neighborLength_up+1);
        
        gDistri_neighborLengthError_ismember_raw_noChoice_collapsed{target_seqLength}(tempi) = temp_plus_down+temp_plus_up;
        
        gDistri_neighborLengthError_ismember_noChoice_collapsed{target_seqLength}(tempi, 1) = ...
            (temp_plus_down+temp_plus_up) / ...
            ( ...
            gDistri_length_noChoice_collapsed{target_seqLength}{tempi}(temp_neighborLength_down+1)...
            +...
            gDistri_length_noChoice_collapsed{target_seqLength}{tempi}(temp_neighborLength_up+1)...
            );
        gDistri_neighborLengthError_ismember_noChoice_collapsed{target_seqLength}(tempi, 2) = ...
            (temp_plus_down) / ...
            ( ...
            gDistri_length_noChoice_collapsed{target_seqLength}{tempi}(temp_neighborLength_down+1)...
            );  
        gDistri_neighborLengthError_ismember_noChoice_collapsed{target_seqLength}(tempi, 3) = ...
            (temp_plus_up) / ...
            ( ...
            gDistri_length_noChoice_collapsed{target_seqLength}{tempi}(temp_neighborLength_up+1)...
            );         
        if target_seqLength == 3 && tempi == 18
            a = 1;
        end
    end    
end
gDistri_neighborLengthError_ismember_noChoice_collapsed;
% ignore length1 down & length4 up
gDistri_neighborLengthError_ismember_noChoice_collapsed{1}(:, 2) = NaN;
gDistri_neighborLengthError_ismember_noChoice_collapsed{1}(:, 1) = ...
    gDistri_neighborLengthError_ismember_noChoice_collapsed{1}(:, 3);
gDistri_neighborLengthError_ismember_noChoice_collapsed{4}(:, 3) = NaN;
gDistri_neighborLengthError_ismember_noChoice_collapsed{4}(:, 1) = ...
    gDistri_neighborLengthError_ismember_noChoice_collapsed{4}(:, 2);


gAcc_noChoice_inOrder;
gAcc_noChoice_inOrder_collapsed = cell(1, pointKindsNum);
for target_seqLength=1:pointKindsNum
    gAcc_noChoice_inOrder_collapsed{target_seqLength} = zeros(fileSize, target_seqLength);
    for tempk=1:target_seqLength
        temp = zeros(fileSize, 1);
        for file_index=1:fileSize
            temp(file_index) = gAcc_noChoice_inOrder{file_index, target_seqLength}(tempk);            
        end
        gAcc_noChoice_inOrder_collapsed{target_seqLength}(:, tempk) = temp;
    end
end

% gCorrect_trial_count_collapsed = target_seqSet;
% g_trial_count_collapsed = target_seqSet;
% gAcc_collapsed = target_seqSet;


%gError_merged
gError_merged = offloadingProb_merged;
gAcc_merged = offloadingProb_merged;
gCorrect_trial_count_merged = offloadingProb_merged;
g_trial_count_merged = offloadingProb_merged;
for fileMerged_index=1:fileMerged_size
    for target_seqLength=1:pointKindsNum
        current_numSeq = length(target_seqSet{target_seqLength}); 
        gCorrect_trial_count_merged{fileMerged_index, target_seqLength} = zeros(current_numSeq,1);
        g_trial_count_merged{fileMerged_index, target_seqLength} = zeros(current_numSeq,1);
        for temp_index=1:fileMergeNum
            
            gCorrect_trial_count_merged{fileMerged_index, target_seqLength} = gCorrect_trial_count_merged{fileMerged_index, target_seqLength} ...
                + gCorrect_trial_count{(fileMerged_index-1)*fileMergeNum+temp_index, target_seqLength};
            g_trial_count_merged{fileMerged_index, target_seqLength} = g_trial_count_merged{fileMerged_index, target_seqLength} ...
                + g_trial_count{(fileMerged_index-1)*fileMergeNum+temp_index, target_seqLength};

        end
        
        gAcc_merged{fileMerged_index, target_seqLength} = gCorrect_trial_count_merged{fileMerged_index, target_seqLength}./ ...
            g_trial_count_merged{fileMerged_index, target_seqLength};
        gError_merged{fileMerged_index, target_seqLength} = 1 - gAcc_merged{fileMerged_index, target_seqLength};
    end
end

%gError_choice_merged
gError_choice_merged = offloadingProb_merged;
gAcc_choice_merged = offloadingProb_merged;
gCorrect_choice_trial_count_merged = offloadingProb_merged;
g_choice_trial_count_merged = offloadingProb_merged;
for fileMerged_index=1:fileMerged_size
    for target_seqLength=1:pointKindsNum
        current_numSeq = length(target_seqSet{target_seqLength}); 
        gCorrect_choice_trial_count_merged{fileMerged_index, target_seqLength} = zeros(current_numSeq,1);
        g_choice_trial_count_merged{fileMerged_index, target_seqLength} = zeros(current_numSeq,1);
        for temp_index=1:fileMergeNum
            
            gCorrect_choice_trial_count_merged{fileMerged_index, target_seqLength} = gCorrect_choice_trial_count_merged{fileMerged_index, target_seqLength} ...
                + gCorrect_choice_trial_count{(fileMerged_index-1)*fileMergeNum+temp_index, target_seqLength};
            g_choice_trial_count_merged{fileMerged_index, target_seqLength} = g_choice_trial_count_merged{fileMerged_index, target_seqLength} ...
                + g_choice_trial_count{(fileMerged_index-1)*fileMergeNum+temp_index, target_seqLength};

        end
        
        gAcc_choice_merged{fileMerged_index, target_seqLength} = gCorrect_choice_trial_count_merged{fileMerged_index, target_seqLength}./ ...
            g_choice_trial_count_merged{fileMerged_index, target_seqLength};
        gError_choice_merged{fileMerged_index, target_seqLength} = 1 - gAcc_choice_merged{fileMerged_index, target_seqLength};
    end
end


%gError_noChoice_merged
gError_noChoice_merged = offloadingProb_merged;
gAcc_noChoice_merged = offloadingProb_merged;
gCorrect_noChoice_trial_count_merged = offloadingProb_merged;
g_noChoice_trial_count_merged = offloadingProb_merged;
for fileMerged_index=1:fileMerged_size
    for target_seqLength=1:pointKindsNum
        current_numSeq = length(target_seqSet{target_seqLength}); 
        gCorrect_noChoice_trial_count_merged{fileMerged_index, target_seqLength} = zeros(current_numSeq,1);
        g_noChoice_trial_count_merged{fileMerged_index, target_seqLength} = zeros(current_numSeq,1);
        for temp_index=1:fileMergeNum
            
            gCorrect_noChoice_trial_count_merged{fileMerged_index, target_seqLength} = gCorrect_noChoice_trial_count_merged{fileMerged_index, target_seqLength} ...
                + gCorrect_noChoice_trial_count{(fileMerged_index-1)*fileMergeNum+temp_index, target_seqLength};
            g_noChoice_trial_count_merged{fileMerged_index, target_seqLength} = g_noChoice_trial_count_merged{fileMerged_index, target_seqLength} ...
                + g_noChoice_trial_count{(fileMerged_index-1)*fileMergeNum+temp_index, target_seqLength};

        end
        
        gAcc_noChoice_merged{fileMerged_index, target_seqLength} = gCorrect_noChoice_trial_count_merged{fileMerged_index, target_seqLength}./ ...
            g_noChoice_trial_count_merged{fileMerged_index, target_seqLength};
        gError_noChoice_merged{fileMerged_index, target_seqLength} = 1 - gAcc_noChoice_merged{fileMerged_index, target_seqLength};
    end
end

gError_collapsed = gAcc_collapsed;
gError_noChoice_collapsed = gAcc_noChoice_collapsed;
gError_choice_collapsed = gAcc_choice_collapsed;
for target_seqLength=1:pointKindsNum
    gAcc_collapsed{target_seqLength} = gCorrect_trial_count_collapsed{target_seqLength}./...
        g_trial_count_collapsed{target_seqLength};    
    gError_collapsed{target_seqLength} = 1 - gAcc_collapsed{target_seqLength};
    
    gAcc_choice_collapsed{target_seqLength} = gCorrect_choice_trial_count_collapsed{target_seqLength}./...
        g_choice_trial_count_collapsed{target_seqLength}; 
    gError_choice_collapsed{target_seqLength} = 1 - gAcc_choice_collapsed{target_seqLength};
    
    gAcc_noChoice_collapsed{target_seqLength} = gCorrect_noChoice_trial_count_collapsed{target_seqLength}./...
        g_noChoice_trial_count_collapsed{target_seqLength}; 
    gError_noChoice_collapsed{target_seqLength} = 1 - gAcc_noChoice_collapsed{target_seqLength};
    
    
    current_numSeq = length(target_seqSet{target_seqLength});

    subplot(pointKindsNum, 1, target_seqLength);
    %nexttile
    
%     plot(1:current_numSeq, gError_collapsed{target_seqLength},'--s', 'LineWidth', 2,'MarkerSize', 10, 'MarkerFaceColor', 'k');
%     plot(1:current_numSeq, gError_noChoice_collapsed{target_seqLength},'--s', 'LineWidth', 2,'MarkerSize', 10, 'MarkerFaceColor', 'k');
    plot(1:current_numSeq, gError_noChoice_collapsed{target_seqLength},'-s', 'LineWidth', 2,'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.25 0.25 0.25]);

    hold on 
    
    for fileMerged_index=1:fileMerged_size
%         scatter(1:current_numSeq, gError_merged{fileMerged_index, target_seqLength}, 25, 'filled', ...
%             'MarkerFaceColor', scatter_rgbColor(fileMerged_index, :));       
        scatter(1:current_numSeq, gError_noChoice_merged{fileMerged_index, target_seqLength}, 10, 'filled', ...
            'MarkerFaceColor', scatter_rgbColor(fileMerged_index, :), ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);        

        hold on
    end           
    
    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
        
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength});%给坐标加标签
    ylabel('Error rate', 'FontSize', 12, 'FontWeight', 'bold');    
    title(sprintf('ErrorRate-noChoice-SeqLength=%d, All-%ddays', target_seqLength, fileSize), 'FontSize',15);    

    % 获取xticklabel的忿
    xtl=get(gca,'XTickLabel');
    % 获取xtick的忿
    xt=get(gca,'XTick');
    % 获取ytick的忿
    yt=get(gca,'YTick');
    % 设置text的x坐标位置仿
    xtextp=xt;
    % 设置text的y坐标位置仿
    ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.08;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',13,'FontWeight','bold');
    % 取消原始ticklabel
    set(gca,'xticklabel','');   
    % 取消右、上边框
    set(gca,'box','off');    
end
gAcc_noChoice_collapsed_inOne = [];
for target_seqLength=1:pointKindsNum  
    gAcc_noChoice_collapsed_inOne = ...
        [gAcc_noChoice_collapsed_inOne gAcc_noChoice_collapsed{target_seqLength}']; 
end

responseMatrix_noChoice = gDistri_lengthDetail_noChoice_trial_count_inOne_array;
responseMatrix_choice = gDistri_lengthDetail_choice_g_trial_count_inOne_array;
        
if ifSave_gAcc == 1
    currentName = ['gAcc', '_'];    
    fileName_gAcc = getFileNameMonkey_MAT_append(subject_name, currentName, monkey_name);
    %fileName_gAcc = getFileNameMonkey_MAT(subject_name, currentName, monkey_name);
    %fileName_gAcc(end-4) = '1';

    save(fileName_gAcc, 'gAcc_noChoice');     
    save(fileName_gAcc, 'gAcc_choice', '-append');   
    save(fileName_gAcc, 'gAcc_noChoice_merged', '-append');     
    save(fileName_gAcc, 'gAcc_choice_merged', '-append');     
    save(fileName_gAcc, 'gAcc_noChoice_collapsed', '-append'); 
    save(fileName_gAcc, 'gAcc_choice_collapsed', '-append');     
    save(fileName_gAcc, 'gAcc_noChoice_collapsed_inOne', '-append');   
    
    save(fileName_gAcc, 'gCorrect_choice_trial_count_collapsed', '-append');   
    save(fileName_gAcc, 'g_choice_trial_count_collapsed', '-append');   
    save(fileName_gAcc, 'gCorrect_noChoice_trial_count_collapsed', '-append');   
    save(fileName_gAcc, 'g_noChoice_trial_count_collapsed', '-append');   
    
    save(fileName_gAcc, 'gDistri_noChoice_collapsed', '-append');
    save(fileName_gAcc, 'gDistri_choice_collapsed', '-append');    
    save(fileName_gAcc, 'gDistri_length_noChoice_collapsed', '-append');
    save(fileName_gAcc, 'gDistri_length_choice_collapsed', '-append');    
    save(fileName_gAcc, 'gDistri_length_noChoice_collapsed_seqAVG', '-append');
    save(fileName_gAcc, 'gDistri_length_choice_collapsed_seqAVG', '-append');
    
    % stimuli_to_response_distri_noChoice = gDistri_lengthDetail_noChoice_trial_count_inOne_array;
    save(fileName_gAcc, 'gDistri_lengthDetail_noChoice_trial_count_inOne_array', '-append');
    % stimuli_to_response_distri_choice = gDistri_lengthDetail_choice_g_trial_count_inOne_array;
    save(fileName_gAcc, 'gDistri_lengthDetail_choice_g_trial_count_inOne_array', '-append');

    
    save(fileName_gAcc, 'seqPrecision_behavior', '-append');   
    
    if false
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\Z_responseMatrix_B2.mat';
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\Z_responseMatrix_2.mat'; %#ok<UNRCH>
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\Z_responseMatrix_Train2.mat';
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\Z_responseMatrix_6C.mat';
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\Z_responseMatrix_6D.mat';
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\Z_responseMatrix_early2p.mat';
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\Z_responseMatrix_all2p.mat';
        fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\Z_responseMatrix_almost2p.mat';
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\D_responseMatrix_2.mat'; %#ok<UNRCH>
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\D_responseMatrix_2p_A.mat';
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\D_responseMatrix_2p_B.mat';
        %fileName = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\D_responseMatrix_Train2.mat';
        save(fileName, 'responseMatrix_noChoice');
        save(fileName, 'responseMatrix_choice', '-append');
    end
end


currentFigName = ['acc', '_'];
% to generate a unique file name for saving figure
fileName_fig4 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(4) == 1
        exportgraphics(fig4, fileName_fig4, 'BackgroundColor', 'none');
end 

%% Plot rAcc, gAcc_choice, gAcc_noChoice, rProb in sequence level, instead of trial level
fig41 = figure(41);
set(gcf,'Position',[0 50 900 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

gAcc_noChoice_collapsed;
gAcc_choice_collapsed;
offloadingProb_all;

mean_gAcc_noChoice_collapsed = zeros(1, seq_length_rangeTail);
sem_gAcc_noChoice_collapsed = zeros(1, seq_length_rangeTail);
mean_gAcc_choice_collapsed = zeros(1, seq_length_rangeTail);
sem_gAcc_choice_collapsed = zeros(1, seq_length_rangeTail);
mean_offloadingProb_all = zeros(1, seq_length_rangeTail);
sem_offloadingProb_all = zeros(1, seq_length_rangeTail);

for tempi=1:seq_length_rangeTail
    mean_gAcc_noChoice_collapsed(tempi) = mean(gAcc_noChoice_collapsed{tempi},'omitnan');
    sem_gAcc_noChoice_collapsed(tempi) = std(gAcc_noChoice_collapsed{tempi},'omitnan')/sqrt(length(gAcc_noChoice_collapsed{tempi}));    
    mean_gAcc_choice_collapsed(tempi) = mean(gAcc_choice_collapsed{tempi},'omitnan');
    sem_gAcc_choice_collapsed(tempi) = std(gAcc_choice_collapsed{tempi},'omitnan')/sqrt(length(gAcc_choice_collapsed{tempi}));    
    mean_offloadingProb_all(tempi) = mean(offloadingProb_all{tempi},'omitnan');
    sem_offloadingProb_all(tempi) = std(offloadingProb_all{tempi},'omitnan')/sqrt(length(offloadingProb_all{tempi}));        
end

%% Plot accuracy
nexttile
errorbar(seq_length_rangeHead: seq_length_rangeTail, mean_gAcc_noChoice_collapsed,sem_gAcc_noChoice_collapsed, '-o', 'Color', [0.25 0.25 0.25 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
hold on
errorbar(seq_length_rangeHead: seq_length_rangeTail, mean_gAcc_choice_collapsed,sem_gAcc_choice_collapsed, '-o', 'Color', [0.25 0.75 0.25 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
hl = legend('internal|noChoice',...
    'internal|Choice',...
    'Location','southwest','fontsize',13);

set(hl,'box','off');
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
set(gca,'box','off');% 取消右、上边框
xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
title(sprintf('Accuracy in sequence level'));

%% Plot rProb
nexttile
errorbar(seq_length_rangeHead: seq_length_rangeTail, mean_offloadingProb_all,sem_offloadingProb_all,'-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'CapSize', 12, 'MarkerSize', 5);
legend('offloading rate','Location','northwest')

ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
set(gca,'box','off');% 取消右、上边框
xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');
title(sprintf('Offloading rate in sequence level'));



%% Plot gAcc_choice
fig44 = figure(44);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

gError_choice_collapsed;
gError_choice_merged;
for target_seqLength=1:pointKindsNum

    current_numSeq = length(target_seqSet{target_seqLength});

    subplot(pointKindsNum, 1, target_seqLength);
    
    plot(1:current_numSeq, gError_choice_collapsed{target_seqLength},'-s', 'LineWidth', 2,'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.25 0.75 0.25]);

    hold on 
    
    for fileMerged_index=1:fileMerged_size
        scatter(1:current_numSeq, gError_choice_merged{fileMerged_index, target_seqLength}, 10, 'filled', ...
            'MarkerFaceColor', scatter_rgbColor(fileMerged_index, :), ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);        

        hold on
    end           
    
    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
        
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength});%给坐标加标签
    ylabel('Error rate', 'FontSize', 12, 'FontWeight', 'bold');    
    title(sprintf('ErrorRate-choice-SeqLength=%d, All-%ddays', target_seqLength, fileSize), 'FontSize',15);    

    % 获取xticklabel的忿
    xtl=get(gca,'XTickLabel');
    % 获取xtick的忿
    xt=get(gca,'XTick');
    % 获取ytick的忿
    yt=get(gca,'YTick');
    % 设置text的x坐标位置仿
    xtextp=xt;
    % 设置text的y坐标位置仿
    ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.08;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',13,'FontWeight','bold');
    % 取消原始ticklabel
    set(gca,'xticklabel','');   
    % 取消右、上边框
    set(gca,'box','off');    
end
currentFigName = ['accChoice', '_'];
% to generate a unique file name for saving figure
fileName_fig44 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if false
        exportgraphics(fig44, fileName_fig44, 'BackgroundColor', 'none');
end 

%% %%%%%%%%%%%%%%%%%%%% g|choice - g|noChoice ,precise %%%%%%%%%%%%%%%%%%%
fig5 = figure(5);
set(gcf,'Position',[0 50 700 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

gCorrect_choice_trial_count;
g_choice_trial_count;
gAcc_choice;

gCorrect_noChoice_trial_count;
g_noChoice_trial_count;
gAcc_noChoice;

gAcc_choice_minus_noChoice = gAcc_choice;
for file_index=1:fileSize
    for target_seqLength=1:pointKindsNum
        gAcc_choice_minus_noChoice{file_index, target_seqLength} = gAcc_choice{file_index, target_seqLength} ...
            - gAcc_noChoice{file_index, target_seqLength};
    end
end

bound1 = 0.2;
bound2 = 0.8;

temp_RoiSeqSetLable = sorted_collapsed_seqSet_inOne(find(sorted_collapsed_offloadingProb>bound1 & sorted_collapsed_offloadingProb<bound2)); %#ok<*FNDSB>
temp_RoiSeqIndex = sorted_index_collapsed_offloadingProb(find(sorted_collapsed_offloadingProb>bound1 & sorted_collapsed_offloadingProb<bound2));

% gAcc_choice_minus_noChoice_RoiSeq = cell(fileSize, 1);
gAcc_choice_minus_noChoice_RoiSeq = zeros(fileSize, length(temp_RoiSeqIndex));

gAcc_choice_RoiSeq_collapsed= zeros(1, length(temp_RoiSeqIndex));
gAcc_noChoice_RoiSeq_collapsed = zeros(1, length(temp_RoiSeqIndex));

gCorrect_choice_RoiSeq_collapsed= zeros(1, length(temp_RoiSeqIndex));
g_choice_RoiSeq_collapsed= zeros(1, length(temp_RoiSeqIndex));

gCorrect_noChoice_RoiSeq_collapsed= zeros(1, length(temp_RoiSeqIndex));
g_noChoice_RoiSeq_collapsed= zeros(1, length(temp_RoiSeqIndex));


% current_numSeq = zeros(1, pointKindsNum);
accumulate_current_numSeq = zeros(1, pointKindsNum);
for target_seqLength=1:pointKindsNum
    if target_seqLength == 1
        accumulate_current_numSeq(target_seqLength) = length(target_seqSet{target_seqLength});
    else
    accumulate_current_numSeq(target_seqLength) = length(target_seqSet{target_seqLength}) + ...
        accumulate_current_numSeq(target_seqLength-1);
    end
    
end


temp_RoiSeqSetLable;
temp_target_seqLength = zeros(1, length(temp_RoiSeqSetLable));
for tempj=1:length(temp_RoiSeqIndex)
    if temp_RoiSeqIndex(tempj) <= accumulate_current_numSeq(1)
        temp_target_seqLength(tempj) = 1;
    elseif temp_RoiSeqIndex(tempj) <= accumulate_current_numSeq(2)
        temp_target_seqLength(tempj) = 2;
    elseif temp_RoiSeqIndex(tempj) <= accumulate_current_numSeq(3)
        temp_target_seqLength(tempj) = 3;
    elseif temp_RoiSeqIndex(tempj) <= accumulate_current_numSeq(4)
        temp_target_seqLength(tempj) = 4;
    end
    
    if temp_target_seqLength(tempj) == 1
        current_numSeq = temp_RoiSeqIndex(tempj);
    else
        current_numSeq = temp_RoiSeqIndex(tempj) - accumulate_current_numSeq(temp_target_seqLength(tempj)-1);
    end
    
    gAcc_choice_RoiSeq_collapsed(tempj) = gAcc_choice_collapsed{temp_target_seqLength(tempj)}(current_numSeq);
    gAcc_noChoice_RoiSeq_collapsed(tempj) = gAcc_noChoice_collapsed{temp_target_seqLength(tempj)}(current_numSeq);
    
    gCorrect_choice_RoiSeq_collapsed(tempj) = gCorrect_choice_trial_count_collapsed{temp_target_seqLength(tempj)}(current_numSeq);
    g_choice_RoiSeq_collapsed(tempj) = g_choice_trial_count_collapsed{temp_target_seqLength(tempj)}(current_numSeq);
    
    gCorrect_noChoice_RoiSeq_collapsed(tempj) = gCorrect_noChoice_trial_count_collapsed{temp_target_seqLength(tempj)}(current_numSeq);
    g_noChoice_RoiSeq_collapsed(tempj) = g_noChoice_trial_count_collapsed{temp_target_seqLength(tempj)}(current_numSeq);
    
    for file_index=1:fileSize
        gAcc_choice_minus_noChoice_RoiSeq(file_index, tempj) = ...
            gAcc_choice_minus_noChoice{file_index, temp_target_seqLength(tempj)}(current_numSeq);
    
    end
    
end
temp_delta_collapsed = gAcc_choice_RoiSeq_collapsed - gAcc_noChoice_RoiSeq_collapsed;
temp_delta = zeros(1, pointKindsNum);
temp_choice = zeros(1, pointKindsNum);
temp_noChoice = zeros(1, pointKindsNum);

mean(gAcc_choice_minus_noChoice_RoiSeq);
mean(gAcc_choice_RoiSeq_collapsed - gAcc_noChoice_RoiSeq_collapsed);

for target_seqLength=1:pointKindsNum
    temp_delta(target_seqLength) = mean(temp_delta_collapsed(find(temp_target_seqLength==target_seqLength)));    
    temp_choice(target_seqLength) = sum(gCorrect_choice_RoiSeq_collapsed(find(temp_target_seqLength==target_seqLength))) ...
        ./ sum(g_choice_RoiSeq_collapsed(find(temp_target_seqLength==target_seqLength))); 
    
    temp_noChoice(target_seqLength) = sum(gCorrect_noChoice_RoiSeq_collapsed(find(temp_target_seqLength==target_seqLength))) ...
        ./ sum(g_noChoice_RoiSeq_collapsed(find(temp_target_seqLength==target_seqLength))); 
    
end
temp_delta;
temp_choice - temp_noChoice;

plot(seq_length_rangeHead: seq_length_rangeTail, temp_choice - temp_noChoice, '-o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 5);
hold on 
plot([seq_length_rangeHead-1 seq_length_rangeHead: seq_length_rangeTail seq_length_rangeTail+1]...
    , zeros(1, pointKindsNum+2), '--');

ylim([-0.3 0.3]);
set(gca, 'FontSize', 20)
set(gca,'YTick',[-0.3:0.1:0.3]);%设置要显示坐标刻度的范围
set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
set(gca,'box','off');% 取消右、上边框
xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
title({'\fontsize{20}{\bf Select balanced sequence(offloading rate=0.2~0.8)}';
    '\fontsize{18}{\color[rgb]{0.25 0.75 0.25}internal|Choice} - internal|noChoice'});%,'FontWeight','bold');
% title({'\fontsize{20}{\bf Select balanced sequence(offloading rate=0.2~0.8)}';
%     '\fontsize{18}{\color{green}internal|Choice} - internal|noChoice'});%,'FontWeight','bold');
% title({'\fontsize{15}EndingHold RT-average';'All trials'});
% title('\fontsize{20}{\bf Select balanced sequence(offloading rate=0.2~0.8)} aa');%,'FontWeight','bold');

currentFigName = ['choiceMinusNoChoice_RoiSeq', '_'];
% to generate a unique file name for saving figure
fileName_fig5 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(5) == 1
        exportgraphics(fig5, fileName_fig5, 'BackgroundColor', 'none');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% scatter_rgbColor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig6 = figure(6);
rgbColor = zeros(1,fileMerged_size,3);
rgbColor(1,:,:) = scatter_rgbColor;
image(rgbColor);
set(gca,'position',[0 0 1 1]);
set(gca,'XTick',[0:0:0])
set(gca,'YTick',[0:0:0])

%% %%%%%%%%%%%%%%%%%%%% endingHold RT %%%%%%%%%%%%%%%%%%%
fig7 = figure(7);

for file_index=1:fileSize
    temp_load = eval(['MAT_file_load.file',sprintf('%d',file_index)]); 
    
    trial_num = temp_load.trial_para.trial_count;
    cumulativeError = temp_load.trial_para.cumulativeError;
    isCorrect = temp_load.trial_para.isCorrect;
    ifFreeTouch = temp_load.basic_para.ifFreeTouch;
    ifSelectOffloading_temp = temp_load.trial_para.ifSelectOffloading;
%     randArray2 = temp_load.basic_para.randArray2;
%     choiceCondition_boardLine1 = temp_load.basic_para.choiceCondition_boardLine1;
%     choiceCondition_boardLine2 = temp_load.basic_para.choiceCondition_boardLine2;
    
    selectFlag_newSeq = zeros(1, trial_num);
    %choiceCondition_flag = zeros(1, trial_num);
    
    for j=1:trial_num
        %to get selectFlag_newSeq        
        if isCorrect(j) == 1 && cumulativeError(j) == 0
            selectFlag_newSeq(j) = 1;
        elseif isCorrect(j) == -1 && cumulativeError(j) == 1
            selectFlag_newSeq(j) = 1;
        end
    end         
    
    choiceCondition_flag = temp_load.trial_para.choiceCondition_flag;
    
%     for j=1:trial_num
%         %to get choiceCondition_flag
%         if selectFlag_newSeq(j) == 1
%             if randArray2(j) <= choiceCondition_boardLine1
%                 choiceCondition_flag(j) = 0;
%             elseif randArray2(j) <= choiceCondition_boardLine2
%                 choiceCondition_flag(j) = 1;
%             else
%                 choiceCondition_flag(j) = 2;
%             end
%         elseif selectFlag_newSeq(j) == 0
%             choiceCondition_flag(j) = choiceCondition_flag(j-1); 
%         end
%     end    
    
    for target_seqLength=1:pointKindsNum
        tempIndex = find(temp_load.trial_para.seq_length(1: trial_num) == (target_seqLength+seq_length_rangeHead-1));
                        
        tempIndex_internalNoChoice = tempIndex(choiceCondition_flag(tempIndex) == 0);                
        tempIndex_allChoice = tempIndex(choiceCondition_flag(tempIndex) == 2);
        tempIndex_internalChoice = tempIndex_allChoice(ifSelectOffloading_temp(tempIndex_allChoice) == -1);                
        tempIndex_internalChoice_newSeq = tempIndex_internalChoice(selectFlag_newSeq(tempIndex_internalChoice) == 1);
        
        tempIndex_internal = sort([tempIndex_internalNoChoice tempIndex_internalChoice_newSeq]);
        
        %error stop correct g trials & free touch g trials (exclude error stop error trials)         
        temp1 = ifFreeTouch(tempIndex_internal) == 0 & isCorrect(tempIndex_internal) == -1;        
        tempIndex_internal_valid1 = tempIndex_internal(~temp1);
        
%         %exclude last item trials
%         temp_valid2 = ones(1, length(tempIndex_internal_valid1));
%         for tempi=1:length(tempIndex_internal_valid1)
%             if temp_load.trial_para.currentSequence{tempIndex_internal_valid1(tempi)}(end) == ...
%                     temp_load.basic_para.numFrames_rangeTail
%                 
%                 temp_valid2(tempi) = 0;
%             end
%         end
% %         tempIndex_internal_valid2 = tempIndex_internal_valid1(logical(temp_valid2));
% 
%         %only include last item trials
%         temp_valid3 = zeros(1, length(tempIndex_internal_valid1));
%         for tempi=1:length(tempIndex_internal_valid1)
%             if temp_load.trial_para.currentSequence{tempIndex_internal_valid1(tempi)}(end) == ...
%                     temp_load.basic_para.numFrames_rangeTail
%                 
%                 temp_valid3(tempi) = 1;
%             end
%         end
%         tempIndex_internal_valid3 = tempIndex_internal_valid1(logical(temp_valid3));        


        %exclude last item trials (from perspective of answer, or isFilled)
        temp_valid2 = ones(1, length(tempIndex_internal_valid1));
        for tempi=1:length(tempIndex_internal_valid1)
            if temp_load.trial_para.isFilled{tempIndex_internal_valid1(tempi)}(end) == 1
                
                temp_valid2(tempi) = 0;
            end
        end
        tempIndex_internal_valid2 = tempIndex_internal_valid1(logical(temp_valid2));

        
        %only include last item trials (from perspective of answer, or isFilled)
        temp_valid3 = zeros(1, length(tempIndex_internal_valid1));
        for tempi=1:length(tempIndex_internal_valid1)
            if temp_load.trial_para.isFilled{tempIndex_internal_valid1(tempi)}(end) == 1
                
                temp_valid3(tempi) = 1;
            end
        end
        tempIndex_internal_valid3 = tempIndex_internal_valid1(logical(temp_valid3));        
        
        
        
        tempIndex_internal_validFinal = tempIndex_internal_valid1;%all trials
%         tempIndex_internal_validFinal = tempIndex_internal_valid2;%exclude last item trials
%         tempIndex_internal_validFinal = tempIndex_internal_valid3;%only include last item trials
     
        tempIndex_internal_correct = tempIndex_internal_validFinal(isCorrect(tempIndex_internal_validFinal) == 1);
        tempIndex_internal_error = tempIndex_internal_validFinal(isCorrect(tempIndex_internal_validFinal) == -1);
        endingHold_correct{file_index, target_seqLength} = temp_load.trial_para.endingHold_RT(tempIndex_internal_correct);
        endingHold_error{file_index, target_seqLength} = temp_load.trial_para.endingHold_RT(tempIndex_internal_error);        
                
    end    
end

for file_index=1:fileSize
    temp_load = eval(['MAT_file_load.file',sprintf('%d',file_index)]); 
    
    trial_num = temp_load.trial_para.trial_count;
    currentSequence = temp_load.trial_para.currentSequence;
    isCorrect = temp_load.trial_para.isCorrect;
    choiceCondition_flag = temp_load.trial_para.choiceCondition_flag;
    ifSelectOffloading = temp_load.trial_para.ifSelectOffloading;
    endingHold_RT = temp_load.trial_para.endingHold_RT;
    target_seqSet;
    
    for tempi=1:trial_num
        temp_seq = currentSequence{tempi};
        temp_length = length(temp_seq);
        for tempj=1:numSeq(temp_length)
            temp_seq2 = target_seqSet{temp_length}{tempj};
            if sum(ismember(temp_seq, temp_seq2)) == temp_length
                temp_seq_index = tempj;
                break                
            end
        end
        temp_seq_index;
            
        if ifSelectOffloading(tempi) == 1
            temp_str = 'offload';
        else
            if choiceCondition_flag(tempi) == 2
                temp_str = 'choiceMemory';
            else
                temp_str = 'noChoiceMemory';
            end            
        end
        endingHold.(temp_str).all{temp_length}{temp_seq_index} = ...
            [endingHold.(temp_str).all{temp_length}{temp_seq_index} endingHold_RT(tempi)];
        
        if isCorrect(tempi) == 1
            % correct trial
            endingHold.(temp_str).correct{temp_length}{temp_seq_index} = ...
                [endingHold.(temp_str).correct{temp_length}{temp_seq_index} endingHold_RT(tempi)];
        else
            % error trial
            endingHold.(temp_str).error{temp_length}{temp_seq_index} = ...
                [endingHold.(temp_str).error{temp_length}{temp_seq_index} endingHold_RT(tempi)];           
        end
        
    end     
    
end
endingHold;
%if false
if ifSave_endingHold == 1
    %temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\endingHold_Z'; %#ok<*UNRCH>
    %temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\endingHold_Z_Train1'; %#ok<*UNRCH>
    %temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\endingHold_Z_Train2'; %#ok<*UNRCH>
    %temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\endingHold_Z_6C'; %#ok<*UNRCH>
    temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\endingHold_Z_6D'; %#ok<*UNRCH>    
    %temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\endingHold_D'; %#ok<*UNRCH>
    %temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\endingHold_D_Train2'; %#ok<*UNRCH>
    %temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\endingHold_D_2p_A'; %#ok<*UNRCH>
    %temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\endingHold_D_2p_B'; %#ok<*UNRCH>
    save(temp_str, 'endingHold');
end

endingHold_correct_collapsed = cell(1, pointKindsNum);
endingHold_error_collapsed = cell(1, pointKindsNum);

endingHold_correct_collapsedToOne = 0;
endingHold_error_collapsedToOne = 0;

endingHold_file = zeros(fileSize, 2);

for file_index=1:fileSize

    for target_seqLength=1:pointKindsNum
        
        endingHold_correct_collapsed{target_seqLength} = [endingHold_correct_collapsed{target_seqLength} ...
            endingHold_correct{file_index, target_seqLength}];
        endingHold_error_collapsed{target_seqLength} = [endingHold_error_collapsed{target_seqLength} ...
            endingHold_error{file_index, target_seqLength}];
        
        endingHold_correct_collapsedToOne = [endingHold_correct_collapsedToOne ...
            endingHold_correct{file_index, target_seqLength}]; %#ok<*AGROW>
        endingHold_error_collapsedToOne = [endingHold_error_collapsedToOne ...
            endingHold_error{file_index, target_seqLength}];        
    end
    
    subplot(1,fileSize+1,file_index);
    bar([0 1], [mean([endingHold_correct{file_index,:}]) mean([endingHold_error{file_index,:}])], ...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
    
    endingHold_file(file_index, 1) = mean([endingHold_correct{file_index,:}]);
    endingHold_file(file_index, 2) = mean([endingHold_error{file_index,:}]);

    ylim([0 1.3]);    
    set(gca,'XTickLabel', ["Correct"; "Error"]);%给坐标加标签
    %title(sprintf('EndingHold RT-%s',cell2mat(subject_name(file_index))));
    temp = cell2mat(subject_name(file_index));
    title(sprintf('EndingHold RT-%s',temp(end-4:end)));
    
end
subplot(1,fileSize+1,fileSize+1);
bar([0 1], [mean(endingHold_correct_collapsedToOne) mean(endingHold_error_collapsedToOne)], ...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
ylim([0 1.3]);
set(gca,'XTickLabel', ["Correct"; "Error"]);%给坐标加标签
title('EndingHold RT-average');

[h_endingHold_correct_error,p_endingHold_correct_error]=...
    ttest2(endingHold_file(:,1), endingHold_file(:,2));
           
%% %%%%%%%%%%%%%%%%%%%% endingHold RT average %%%%%%%%%%%%%%%%%%%
fig8 = figure(8);

bar([0 1], [mean(endingHold_file(:,1)) mean(endingHold_file(:,2))], ...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
hold on
errorbar([0 1], mean(endingHold_file),std(endingHold_file)./sqrt(length(endingHold_file)), '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
tempTxt = sprintf('');
if p_endingHold_correct_error < 0.001
    tempTxt = sprintf('***');
elseif p_endingHold_correct_error < 0.01
    tempTxt = sprintf('**');
elseif p_endingHold_correct_error < 0.05
    tempTxt = sprintf('*');
end
text(0.5,max(mean(endingHold_file))+0.1,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
    'HorizontalAlignment','center');

ylim([0 max(mean(endingHold_file))+0.2]);
set(gca,'YTick',[0:0.2:1.3]);%设置要显示坐标刻度的范围
set(gca, 'FontSize', 20)
set(gca,'XTickLabel', ["Correct"; "Error"]);%给坐标加标签
set(gca,'box','off');% 取消右、上边框
ylabel('Reaction time (s)', 'FontSize', 18, 'FontWeight', 'bold');
title({'\fontsize{15}EndingHold RT, internal memory';'Each session'});
% title({'\fontsize{15}EndingHold RT-average';'Exclude last item trials (answer)'});
% title({'\fontsize{15}EndingHold RT-average';'Only include last item trials (answer)'});
% title('EndingHold RT-average');

currentFigName = ['endingHold', '_'];
% to generate a unique file name for saving figure
fileName_fig8 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(8) == 1
        exportgraphics(fig8, fileName_fig8, 'BackgroundColor', 'none');
end


%% Save data of paper fig3
if false == 1
    currentName = ['fig3_data', '_'];
    fileName_figData = getFileNameMonkey_MAT(subject_name, currentName, monkey_name);
    
    fig3_panelA_data = struct;
    fig3_panelA_data.endingHold_file = endingHold_file;
    fig3_panelA_data.p_endingHold_correct_error = p_endingHold_correct_error;
    
    save(fileName_figData, 'fig3_panelA_data');
end
%% %%%%%%%%%%%% correlation between rProb and error rate %%%%%%%%%%%%%%%
fig9 = figure(9);

offloadingProb_all;
gError_noChoice_collapsed;
offloadingProb_all_inOne = zeros(1, length(collapsed_seqSet_inOne));
gError_noChoice_collapsed_inOne = zeros(1, length(collapsed_seqSet_inOne));
offloadingAcc_all_inOne = zeros(1, length(collapsed_seqSet_inOne));

tempj = 1;
for target_seqLength=1:pointKindsNum
    current_numSeq = length(target_seqSet{target_seqLength}); 
    for tempi=1:current_numSeq
    % for tempi=1:length(collapsed_seqSet_inOne)
        offloadingProb_all_inOne(tempj) = offloadingProb_all{target_seqLength}(tempi);
        gError_noChoice_collapsed_inOne(tempj) = gError_noChoice_collapsed{target_seqLength}(tempi);
        offloadingAcc_all_inOne(tempj) = offloadingAcc_all{target_seqLength}(tempi);        
        tempj = tempj + 1;
    end
end

isvalid = (~isnan(gError_noChoice_collapsed_inOne)) & (~isnan(offloadingProb_all_inOne));

x = gError_noChoice_collapsed_inOne(isvalid);
y = offloadingProb_all_inOne(isvalid);
% x = offloadingProb_all_inOne;
% y = gError_noChoice_collapsed_inOne;
n = 1;
[p_mapping,S] = polyfit(x,y,n);
% [p,S] = polyfit(offloadingProb_all_inOne,gError_noChoice_collapsed_inOne,1);
r2 = 1 - (S.normr/norm(y - mean(y)))^2;
[r, p_corr] = corr(x',y');
tempTxt = 'p>0.05';
if p_corr < 0.05    
    tempTxt = 'p<0.05';
end
if p_corr < 0.01
    tempTxt = 'p<0.01';
end
if p_corr < 0.001
    tempTxt = 'p<0.001';
end

[p_mapping_poly,S_poly] = polyfit(x,y,3);
r2_poly = 1 - (S_poly.normr/norm(y - mean(y)))^2;


if ifSaveMapping == 1
    currentName = ['mappingParam', '_'];    
    fileName_mappingParam = getFileNameMonkey_MAT(subject_name, currentName, monkey_name);
    save(fileName_mappingParam, 'p_mapping');
    save(fileName_mappingParam, 'p_mapping_poly','-append');   
    save(fileName_mappingParam, 'target_seqSet','-append');        
    save(fileName_mappingParam, 'seqSet_inOne','-append'); 
end

% y_hat = polyval(p,x);
x_fit = 0:0.1:1;
y_fit = polyval(p_mapping,x_fit);

% temp_color = [0.25 0.75 0.75;
%             0.75 0.25 0.75;
%             0.75 0.75 0.25;
%             0.25 0.25 0.25];
% temp_color = [0.25 0.25 0.25;
%             0.25 0.25 0.25;
%             0.25 0.25 0.25;
%             0.25 0.25 0.25];
% temp_color = [0.85 0.85 0.85;
%             0.6 0.6 0.6;
%             0.25 0.25 0.25;
%             0 0 0];
% temp_color = [[1 0 0];
%             [1 0 0]*0.75;
%             [1 0 0]*0.55;
%             [1 0 0]*0.3];
% my_hot = hot;
% temp_color = [my_hot(floor(256*0.7),:);
%             my_hot(floor(256*0.5),:);
%             my_hot(floor(256*0.25),:);
%             my_hot(floor(256*0.05),:)];
% temp_color = [[1 0 0];
%             [1 0 0];
%             [1 0 0];
%             [1 0 0]];        

temp_color = [[166,97,26]/255;
            [223,194,125]/255;
            [128,205,193]/255;
            [1,133,113]/255];



% temp_faceAlpha = [0.15 0.35 0.7 1];
% temp_faceAlpha = [0.75 0.75 0.75 0.75];
temp_faceAlpha = 0.5;
% temp_size = [10, 30, 50, 70];

for target_seqLength=1:pointKindsNum
%     scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
%         temp_size(target_seqLength), 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
%         'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); 
    
%     scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
%         25, 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
%         'MarkerFaceAlpha', temp_faceAlpha, 'MarkerEdgeAlpha', 0.7);   
    
    temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(gError_noChoice_collapsed{target_seqLength}));    
    scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', temp_faceAlpha, 'MarkerEdgeAlpha', 0.7);       

%     scatter(offloadingProb_all{target_seqLength}, gError_noChoice_collapsed{target_seqLength}, ...
%         25, 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
%         'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
end
% plot(x_fit, y_fit, '-', 'LineWidth', 3, 'Color', [0.25 0.25 0.75]);
plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);

% legend('1 item','2 item','3 item','4 item','fit','Location','southeast','fontsize',12)
% text(0.4,0.4,sprintf('r2 = %.3f',r2), 'fontsize',20,'FontWeight','bold');
% legend('1 item','2 item','3 item','4 item','fit','Location','southeast','fontsize',8)
legend('length 1','length 2','length 3','length 4','fit','Location','southeast','fontsize',12)
% text(0.1,0.9,sprintf('r2 = %.3f',r2), 'fontsize',16,'FontWeight','bold');
text(0.1,0.9,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',14,'FontWeight','bold');


xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Error rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
% title('Correlation between offloading rate & error rate', 'FontSize', 20, 'FontWeight', 'bold'); 
title('Correlation, offloading rate & error rate', 'FontSize', 18, 'FontWeight', 'bold'); 

currentFigName = ['correlation_rProbandError', '_'];
% to generate a unique file name for saving figure
fileName_fig9 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(9) == 1
        exportgraphics(fig9, fileName_fig9, 'BackgroundColor', 'none');
end

if ifSave_paperFig2Data == 1
    currentName = ['fig2_data', '_'];    
    fileName_figData = getFileNameMonkey_MAT_append(subject_name, currentName, monkey_name);
    fig2_panelC_data = struct;
    %fig2_panelC_data.errorRate = x;
    %fig2_panelC_data.offloadingRate = y;
    fig2_panelC_data.errorRate = gError_noChoice_collapsed;
    fig2_panelC_data.offloadingRate = offloadingProb_all;
    fig2_panelC_data.maxLength = pointKindsNum;
    fig2_panelC_data.x_fit = x_fit;
    fig2_panelC_data.y_fit = y_fit;
    fig2_panelC_data.r2 = r2;
    fig2_panelC_data.r = r;
    fig2_panelC_data.p_corr = p_corr;
    save(fileName_figData, 'fig2_panelC_data','-append'); 
end



fig91 = figure(91);
set(gcf,'Position',[0 50 1800 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(1,pointKindsNum,'TileSpacing','Compact','Padding','Compact');

errorRate = gError_noChoice_merged;
offloadingRate = offloadingProb_merged;
maxLength = pointKindsNum;

temp_size = 15;
temp_faceAlpha = 1;%0.7
temp_edgeAlpha = 1;%0.7
for target_seqLength=1:maxLength
    nexttile
    temp_errorRate = [];
    temp_offloadingRate = [];
    
    for fileMerged_index=1:fileMerged_size
        if fileMerged_index == 1 || fileMerged_index == fileMerged_size
            temp_size2 = 1.5*temp_size;
        else
            temp_size2 = temp_size;
        end
        
        scatter(errorRate{fileMerged_index,target_seqLength}, offloadingRate{fileMerged_index,target_seqLength}, ...
            temp_size2, 'filled', 'MarkerFaceColor', scatter_rgbColor(fileMerged_index, :), ...
            'MarkerFaceAlpha', temp_faceAlpha, 'MarkerEdgeAlpha', temp_edgeAlpha);
        hold on
        temp_errorRate = [temp_errorRate errorRate{fileMerged_index,target_seqLength}]; %#ok<*AGROW>
        temp_offloadingRate = [temp_offloadingRate offloadingRate{fileMerged_index,target_seqLength}];
    end
    
    temp_numSeq = numSeq(target_seqLength);
    for tempi=1:temp_numSeq
        %plot(temp_errorRate(tempi,:), temp_offloadingRate(tempi,:),'Color',[0 0 0]);
        plot(temp_errorRate(tempi,:), temp_offloadingRate(tempi,:));
        hold on
    end
    
    xlim([0 1]);
    ylim([0 1]);
    set(gca, 'FontSize', 14)
    set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    set(gca,'box','off');% 取消右、上边框
    xlabel('Error rate', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Offloading rate', 'FontSize', 14, 'FontWeight', 'bold');
    
    temp_str1 = 'Correlation, offloading rate & error rate';
    temp_str2 = sprintf('length = %d', target_seqLength);
    [t,s] = title(temp_str1,temp_str2);
    t.FontSize = 14;
    t.FontWeight = 'bold';
    s.FontSize = 12;
end
currentFigName = ['correlation_rProbandError_eachMerge', '_'];
% to generate a unique file name for saving figure
fileName_fig91 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if false
        exportgraphics(fig91, fileName_fig91, 'BackgroundColor', 'none');
end

%% %%%%%%%%%%%% correlation between rProb and error rate, low & high error rate %%%%%%%%%%%%%%%
fig10 = figure(10);

offloadingProb_all_inOne;
gError_noChoice_collapsed_inOne;
gError_noChoice_collapsed;

temp1 = isnan(gError_noChoice_collapsed_inOne);
tempValidSeq_sorted = sort(gError_noChoice_collapsed_inOne(~temp1));
errorRateThreshold_big = tempValidSeq_sorted(floor(length(tempValidSeq_sorted)/2))+0.0001;
errorRateThreshold_small = tempValidSeq_sorted(floor(length(tempValidSeq_sorted)/2))-0.0001;
tempResidual_big = abs(sum(tempValidSeq_sorted<errorRateThreshold_big) - floor(length(tempValidSeq_sorted)/2));
tempResidual_small = abs(sum(tempValidSeq_sorted<errorRateThreshold_small) - floor(length(tempValidSeq_sorted)/2));
if tempResidual_big < tempResidual_small
    errorRateThreshold = errorRateThreshold_big;
else
    errorRateThreshold = errorRateThreshold_small;
end

temp_color = [0.25 0.75 0.75;
            0.75 0.25 0.75;
            0.75 0.75 0.25;
            0.25 0.25 0.25];
for target_seqLength=1:pointKindsNum
    scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
        25, 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
end
% plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
hold on
plot([errorRateThreshold errorRateThreshold], [0 1], '-', 'LineWidth', 3, 'Color', [0.75 0.25 0.25]);
text(0.1,0.9,sprintf('50%%Threshold = %.4f',errorRateThreshold), 'fontsize',16,'FontWeight','bold', 'Color', [0.75 0.25 0.25]);

if pointKindsNum == 4
    legend('1 item','2 item','3 item','4 item','Threshold','Location','southeast','fontsize',8)
elseif pointKindsNum == 3
    legend('1 item','2 item','3 item','Threshold','Location','southeast','fontsize',8)
end
% text(0.1,0.9,sprintf('r2 = %.3f',r2), 'fontsize',16,'FontWeight','bold');

xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Error rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
title('Correlation, offloading rate & error rate', 'FontSize', 18, 'FontWeight', 'bold'); 


currentFigName = ['errorRate_threshold', '_'];
% to generate a unique file name for saving figure
fileName_fig10 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(10) == 1
        exportgraphics(fig10, fileName_fig10, 'BackgroundColor', 'none');
end


fig11 = figure(11);

% for low & high error rate
seqSetIndex_lowErrorRate = cell(1, pointKindsNum);
seqSetIndex_highErrorRate = cell(1, pointKindsNum);
for target_seqLength=1:pointKindsNum
    seqSetIndex_lowErrorRate{target_seqLength} = gError_noChoice_collapsed{target_seqLength} < errorRateThreshold;
    seqSetIndex_highErrorRate{target_seqLength} = gError_noChoice_collapsed{target_seqLength} >= errorRateThreshold;
end
seqSetIndex_lowErrorRate_collapsed = gError_noChoice_collapsed_inOne < errorRateThreshold;
seqSetIndex_highErrorRate_collapsed = gError_noChoice_collapsed_inOne >= errorRateThreshold;
  


% % for seq_length=1,2 & seq_length=3,4
% seqSetIndex_lowErrorRate = cell(1, pointKindsNum);
% seqSetIndex_highErrorRate = cell(1, pointKindsNum);
% for target_seqLength=1:pointKindsNum
%     if target_seqLength <= 2
%         seqSetIndex_lowErrorRate{target_seqLength} = true(1, length(gError_noChoice_collapsed{target_seqLength}));
%         seqSetIndex_highErrorRate{target_seqLength} = ~seqSetIndex_lowErrorRate{target_seqLength};
%     else
%         seqSetIndex_lowErrorRate{target_seqLength} = false(1, length(gError_noChoice_collapsed{target_seqLength}));
%         seqSetIndex_highErrorRate{target_seqLength} = ~seqSetIndex_lowErrorRate{target_seqLength};
%     end
% end
% seqSetIndex_lowErrorRate_collapsed = false(1, length(offloadingProb_all_inOne));
% seqSetIndex_lowErrorRate_collapsed(1:length(offloadingProb_all{1}) + length(offloadingProb_all{2})) = 1;
% seqSetIndex_highErrorRate_collapsed = ~seqSetIndex_lowErrorRate_collapsed;


x = gError_noChoice_collapsed_inOne(seqSetIndex_lowErrorRate_collapsed);
y = offloadingProb_all_inOne(seqSetIndex_lowErrorRate_collapsed);
n = 1;
[p,S] = polyfit(x,y,n);
r2 = 1 - (S.normr/norm(y - mean(y)))^2;

r = corr(x',y');

x_fit_low = 0:0.1:1;
y_fit_low = polyval(p,x_fit_low);


temp_color = [0.25 0.75 0.75;
            0.75 0.25 0.75;
            0.75 0.75 0.25;
            0.25 0.25 0.25];
        
for target_seqLength=1:pointKindsNum
    scatter(gError_noChoice_collapsed{target_seqLength}(seqSetIndex_lowErrorRate{target_seqLength})...
        , offloadingProb_all{target_seqLength}(seqSetIndex_lowErrorRate{target_seqLength})...
        , ...
        25, 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
end
               
% plot(x_fit_low, y_fit_low, '-', 'LineWidth', 3, 'Color', [0.25 0.25 0.75]);
plot(x_fit_low, y_fit_low, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);

if pointKindsNum == 4
    legend('1 item','2 item','3 item','4 item','fit','Location','southeast','fontsize',8)
elseif pointKindsNum == 3
    legend('1 item','2 item','3 item','fit','Location','southeast','fontsize',8)
end

text(0.1,0.9,sprintf(' r2 = %.3f\n r = %.3f',r2,r), 'fontsize',16,'FontWeight','bold');



xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Error rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
title('Correlation, low error rate', 'FontSize', 20, 'FontWeight', 'bold'); 
% title('Correlation, seqLength = 1,2', 'FontSize', 20, 'FontWeight', 'bold'); 

currentFigName = ['correlation_rProbandError_lowErrorRate', '_'];
% currentFigName = ['correlation_rProbandError_length12', '_'];
% to generate a unique file name for saving figure
fileName_fig11 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(11) == 1
        exportgraphics(fig11, fileName_fig11, 'BackgroundColor', 'none');
end


fig12 = figure(12);

seqSetIndex_highErrorRate_collapsed;

x = gError_noChoice_collapsed_inOne(seqSetIndex_highErrorRate_collapsed);
y = offloadingProb_all_inOne(seqSetIndex_highErrorRate_collapsed);
n = 1;
[p,S] = polyfit(x,y,n);
r2 = 1 - (S.normr/norm(y - mean(y)))^2;

r = corr(x',y');

x_fit_high = 0:0.1:1;
y_fit_high = polyval(p,x_fit_high);

temp_color = [0.25 0.75 0.75;
            0.75 0.25 0.75;
            0.75 0.75 0.25;
            0.25 0.25 0.25];
        
for target_seqLength=1:pointKindsNum
    scatter(gError_noChoice_collapsed{target_seqLength}(seqSetIndex_highErrorRate{target_seqLength})...
        , offloadingProb_all{target_seqLength}(seqSetIndex_highErrorRate{target_seqLength})...
        , ...
        25, 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
end
               
% plot(x_fit_high, y_fit_high, '-', 'LineWidth', 3, 'Color', [0.25 0.25 0.75]);
plot(x_fit_high, y_fit_high, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
if pointKindsNum == 4
    legend('1 item','2 item','3 item','4 item','fit','Location','southeast','fontsize',8)
elseif pointKindsNum == 3
    legend('1 item','2 item','3 item','fit','Location','southeast','fontsize',8)
end
text(0.1,0.9,sprintf(' r2 = %.3f\n r = %.3f',r2,r), 'fontsize',16,'FontWeight','bold');

xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Error rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
title('Correlation, high error rate', 'FontSize', 20, 'FontWeight', 'bold'); 
% title('Correlation, seqLength = 3,4', 'FontSize', 20, 'FontWeight', 'bold'); 

currentFigName = ['correlation_rProbandError_highErrorRate', '_'];
% currentFigName = ['correlation_rProbandError_length34', '_'];
% to generate a unique file name for saving figure
fileName_fig12 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(12) == 1
        exportgraphics(fig12, fileName_fig12, 'BackgroundColor', 'none');
end

%% %%%%%%%%%%%% rProb, low & high error rate %%%%%%%%%%%%%%%
fig13 = figure(13);

offloadingProb_all_inOne;
gError_noChoice_collapsed_inOne;
gError_noChoice_collapsed;

% for low & high error rate
% seqSetIndex_lowErrorRate = cell(1, pointKindsNum);
% seqSetIndex_highErrorRate = cell(1, pointKindsNum);
% for target_seqLength=1:pointKindsNum
%     seqSetIndex_lowErrorRate{target_seqLength} = gError_noChoice_collapsed{target_seqLength} < 0.5;
%     seqSetIndex_highErrorRate{target_seqLength} = ~seqSetIndex_lowErrorRate{target_seqLength};
% end
seqSetIndex_lowErrorRate_collapsed = gError_noChoice_collapsed_inOne < errorRateThreshold;
seqSetIndex_highErrorRate_collapsed = ~seqSetIndex_lowErrorRate_collapsed;


% xbin = 0:0.1:(1-0.1);
% ybin_count_low = zeros(1, length(xbin));
% tempOffload = offloadingProb_all_inOne(seqSetIndex_lowErrorRate_collapsed);
% for tempi=1:length(tempOffload)
%     temp0 = ceil(tempOffload(tempi) * 10);
%     if temp0 == 0
%         temp0 = 1;
%     end    
%     ybin_count_low(temp0) = ybin_count_low(temp0) + 1;    
% end
% bar(xbin+0.05, ybin_count_low, 0.8);

xbin = 0:0.2:(1-0.1);
ybin_count_low = zeros(1, length(xbin));
tempOffload = offloadingProb_all_inOne(seqSetIndex_lowErrorRate_collapsed);
for tempi=1:length(tempOffload)
    temp0 = ceil(tempOffload(tempi) * 5);
    if temp0 == 0
        temp0 = 1;
    end
    ybin_count_low(temp0) = ybin_count_low(temp0) + 1;    
end
bar(xbin+0.1, ybin_count_low, 0.8, ...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);


xlim([0 1]);
% ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
% set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Count', 'FontSize', 18, 'FontWeight', 'bold');  

title(sprintf('offloading rate distribution, low error rate\n total %d sequences'...
    ,length(tempOffload)), 'FontSize', 20, 'FontWeight', 'bold');

currentFigName = ['rProbDistri_lowErrorRate', '_'];
% to generate a unique file name for saving figure
fileName_fig13 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(13) == 1
        exportgraphics(fig13, fileName_fig13, 'BackgroundColor', 'none');
end


fig14 = figure(14);

% xbin = 0:0.1:(1-0.1);
% ybin_count_high = zeros(1, length(xbin));
% tempOffload = offloadingProb_all_inOne(seqSetIndex_highErrorRate_collapsed);
% for tempi=1:length(tempOffload)
%     temp0 = ceil(tempOffload(tempi) * 10);
%     if temp0 == 0
%         temp0 = 1;
%     end    
%     ybin_count_high(temp0) = ybin_count_high(temp0) + 1;    
% end
% bar(xbin+0.05, ybin_count_high, 0.8);

xbin = 0:0.2:(1-0.1);
ybin_count_high = zeros(1, length(xbin));
tempOffload = offloadingProb_all_inOne(seqSetIndex_highErrorRate_collapsed);
for tempi=1:length(tempOffload)
    temp0 = ceil(tempOffload(tempi) * 5);
    if temp0 == 0
        temp0 = 1;
    end
    ybin_count_high(temp0) = ybin_count_high(temp0) + 1;    
end
bar(xbin+0.1, ybin_count_high, 0.8, ...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);



xlim([0 1]);
% ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
% set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Count', 'FontSize', 18, 'FontWeight', 'bold');  
% title('offloading rate distribution, high error rate', 'FontSize', 20, 'FontWeight', 'bold'); 
title(sprintf('offloading rate distribution, high error rate\n total %d sequences'...
    ,length(tempOffload)), 'FontSize', 20, 'FontWeight', 'bold');

currentFigName = ['rProbDistri_highErrorRate', '_'];
% to generate a unique file name for saving figure
fileName_fig14 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(14) == 1
        exportgraphics(fig14, fileName_fig14, 'BackgroundColor', 'none');
end

%% %%%%%%%%%%%% lightChunking & heavyChunking, Initialization %%%%%%%%%%%%%%%
numChunk_seqSet = target_seqSet;
for target_seqLength=1:pointKindsNum
    numChunk_seqSet(target_seqLength) = {ones(length(numChunk_seqSet{target_seqLength}), 1)};
end

for target_seqLength=1:pointKindsNum
    for tempi=1:length(target_seqSet{target_seqLength})
        if target_seqLength > 1
            for temp_seqLength=1:target_seqLength-1
                if target_seqSet{target_seqLength}{tempi}(temp_seqLength)  ~= ...
                       target_seqSet{target_seqLength}{tempi}(temp_seqLength+1)-1
                   
                   numChunk_seqSet{target_seqLength}(tempi) = numChunk_seqSet{target_seqLength}(tempi) + 1;
                end
            end
        end
    end
end
numChunk_seqSet_collapsed = [];
for target_seqLength=1:pointKindsNum
    numChunk_seqSet_collapsed = [numChunk_seqSet_collapsed numChunk_seqSet{target_seqLength}'];
end

numChunk_seqSet;
numChunk_seqSet_collapsed;

chunkDegree_seqSet_0light1heavy = cell(1, pointKindsNum);
for target_seqLength=1:pointKindsNum
    chunkDegree_seqSet_0light1heavy{target_seqLength} = -1*ones(length(numChunk_seqSet{target_seqLength}), 1);
    if target_seqLength > 1
        for tempi=1:length(numChunk_seqSet{target_seqLength})
            if numChunk_seqSet{target_seqLength}(tempi) == max(numChunk_seqSet{target_seqLength})
                chunkDegree_seqSet_0light1heavy{target_seqLength}(tempi) = 0;
            %elseif numChunk_seqSet{target_seqLength}(tempi) == min(numChunk_seqSet{target_seqLength})
            elseif numChunk_seqSet{target_seqLength}(tempi) < max(numChunk_seqSet{target_seqLength})
                chunkDegree_seqSet_0light1heavy{target_seqLength}(tempi) = 1;
            end
        end
    end
end
chunkDegree_seqSet_0light1heavy_collapsed = [];
for target_seqLength=1:pointKindsNum
    chunkDegree_seqSet_0light1heavy_collapsed = [chunkDegree_seqSet_0light1heavy_collapsed...
        chunkDegree_seqSet_0light1heavy{target_seqLength}'];
end

%% %%%%%%%%%%%% lightChunking & heavyChunking, Correlation of error rate & rProb %%%%%%%%%%%%%%%
fig15 = figure(15);

chunkDegree_seqSet_0light1heavy;
chunkDegree_seqSet_0light1heavy_collapsed;

seqSetIndex_lightChunk = cell(1, pointKindsNum);
seqSetIndex_heavyChunk = cell(1, pointKindsNum);
for target_seqLength=1:pointKindsNum
    seqSetIndex_lightChunk{target_seqLength} = chunkDegree_seqSet_0light1heavy{target_seqLength}==0;
    seqSetIndex_heavyChunk{target_seqLength} = chunkDegree_seqSet_0light1heavy{target_seqLength}==1;
end

seqSetIndex_lightChunk_collapsed = chunkDegree_seqSet_0light1heavy_collapsed==0;
seqSetIndex_heavyChunk_collapsed = chunkDegree_seqSet_0light1heavy_collapsed==1;

x = gError_noChoice_collapsed_inOne(seqSetIndex_lightChunk_collapsed);
y = offloadingProb_all_inOne(seqSetIndex_lightChunk_collapsed);
n = 1;
[p,S] = polyfit(x,y,n);
r2 = 1 - (S.normr/norm(y - mean(y)))^2;

x_fit_low = 0:0.1:1;
y_fit_low = polyval(p,x_fit_low);


temp_color = [0.25 0.75 0.75;
            0.75 0.25 0.75;
            0.75 0.75 0.25;
            0.25 0.25 0.25];

handle = zeros(1, pointKindsNum+2);        
for target_seqLength=1:pointKindsNum
    handle(target_seqLength) = scatter(gError_noChoice_collapsed{target_seqLength}(seqSetIndex_lightChunk{target_seqLength})...
        , offloadingProb_all{target_seqLength}(seqSetIndex_lightChunk{target_seqLength})...
        , ...
        25, 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
end

mean(offloadingProb_all_inOne(seqSetIndex_lightChunk_collapsed));
mean(gError_noChoice_collapsed_inOne(seqSetIndex_lightChunk_collapsed));

handle(pointKindsNum+1) = scatter(mean(gError_noChoice_collapsed_inOne(seqSetIndex_lightChunk_collapsed))...
    , mean(offloadingProb_all_inOne(seqSetIndex_lightChunk_collapsed))...
    , 200, 'filled', 'square', ...
    'MarkerFaceColor', [0.25 0.75 0.25], ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
hold on
plot([mean(gError_noChoice_collapsed_inOne(seqSetIndex_lightChunk_collapsed)) ...
    mean(gError_noChoice_collapsed_inOne(seqSetIndex_lightChunk_collapsed))],...
    [0 mean(offloadingProb_all_inOne(seqSetIndex_lightChunk_collapsed))]...        
    , '--', 'LineWidth', 2, 'Color', [0.25 0.75 0.25 0.7]);
hold on
plot([0 mean(gError_noChoice_collapsed_inOne(seqSetIndex_lightChunk_collapsed))],...
    [mean(offloadingProb_all_inOne(seqSetIndex_lightChunk_collapsed))...
    mean(offloadingProb_all_inOne(seqSetIndex_lightChunk_collapsed))]...
    , '--', 'LineWidth', 2, 'Color', [0.25 0.75 0.25 0.7]);
hold on

% handle(6) = plot(x_fit_low, y_fit_low, '-', 'LineWidth', 3, 'Color', [0.25 0.25 0.75 0.7]);
handle(pointKindsNum+2) = plot(x_fit_low, y_fit_low, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
if pointKindsNum == 4
    legend(handle, '1 item','2 item','3 item','4 item','centroid','fit','Location','southeast','fontsize',8)
elseif pointKindsNum == 3
    legend(handle, '1 item','2 item','3 item','centroid','fit','Location','southeast','fontsize',8)
end
text(0.1,0.9,sprintf('r2 = %.3f',r2), 'fontsize',16,'FontWeight','bold');

xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Error rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
title('Correlation, light degree of chunking', 'FontSize', 20, 'FontWeight', 'bold'); 
% title('Correlation, no chunking sequences', 'FontSize', 20, 'FontWeight', 'bold'); 

currentFigName = ['correlation_rProbandError_lightChunk', '_'];
% to generate a unique file name for saving figure
fileName_fig15 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(15) == 1
        exportgraphics(fig15, fileName_fig15, 'BackgroundColor', 'none');
end

fig16 = figure(16);

seqSetIndex_heavyChunk_collapsed;

x = gError_noChoice_collapsed_inOne(seqSetIndex_heavyChunk_collapsed);
y = offloadingProb_all_inOne(seqSetIndex_heavyChunk_collapsed);
n = 1;
[p,S] = polyfit(x,y,n);
r2 = 1 - (S.normr/norm(y - mean(y)))^2;

x_fit_high = 0:0.1:1;
y_fit_high = polyval(p,x_fit_high);

temp_color = [0.25 0.75 0.75;
            0.75 0.25 0.75;
            0.75 0.75 0.25;
            0.25 0.25 0.25];

handle = zeros(1, pointKindsNum+2);
for target_seqLength=1:pointKindsNum
    handle(target_seqLength) = scatter(gError_noChoice_collapsed{target_seqLength}(seqSetIndex_heavyChunk{target_seqLength})...
        , offloadingProb_all{target_seqLength}(seqSetIndex_heavyChunk{target_seqLength})...
        , ...
        25, 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
end
               

mean(offloadingProb_all_inOne(seqSetIndex_heavyChunk_collapsed));
mean(gError_noChoice_collapsed_inOne(seqSetIndex_heavyChunk_collapsed));

handle(pointKindsNum+1) = scatter(mean(gError_noChoice_collapsed_inOne(seqSetIndex_heavyChunk_collapsed))...
    , mean(offloadingProb_all_inOne(seqSetIndex_heavyChunk_collapsed))...
    , 200, 'filled', 'square', ...
    'MarkerFaceColor', [0.25 0.75 0.25], ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
hold on
plot([mean(gError_noChoice_collapsed_inOne(seqSetIndex_heavyChunk_collapsed)) ...
    mean(gError_noChoice_collapsed_inOne(seqSetIndex_heavyChunk_collapsed))],...
    [0 mean(offloadingProb_all_inOne(seqSetIndex_heavyChunk_collapsed))]...    
    , '--', 'LineWidth', 2, 'Color', [0.25 0.75 0.25 0.7]);
hold on
plot([0 mean(gError_noChoice_collapsed_inOne(seqSetIndex_heavyChunk_collapsed))],...
    [mean(offloadingProb_all_inOne(seqSetIndex_heavyChunk_collapsed))...
    mean(offloadingProb_all_inOne(seqSetIndex_heavyChunk_collapsed))]...    
    , '--', 'LineWidth', 2, 'Color', [0.25 0.75 0.25 0.7]);
hold on

% handle(6) = plot(x_fit_low, y_fit_low, '-', 'LineWidth', 3, 'Color', [0.25 0.25 0.75 0.7]);
handle(pointKindsNum+2) = plot(x_fit_high, y_fit_high, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);

if pointKindsNum == 4
    legend(handle, '1 item','2 item','3 item','4 item','centroid','fit','Location','southeast','fontsize',8)
elseif pointKindsNum == 3
    legend(handle, '1 item','2 item','3 item','centroid','fit','Location','southeast','fontsize',8)
end
text(0.1,0.9,sprintf('r2 = %.3f',r2), 'fontsize',16,'FontWeight','bold');




xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Error rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
title('Correlation, heavy degree of chunking', 'FontSize', 20, 'FontWeight', 'bold'); 
% title('Correlation, chunking sequences', 'FontSize', 20, 'FontWeight', 'bold'); 

collapsed_seqSet_inOne;
heavyChunk_seqCluster1_index = ((offloadingProb_all_inOne<0.5)+seqSetIndex_heavyChunk_collapsed)==2;
heavyChunk_seqCluster2_index = ((offloadingProb_all_inOne>=0.5)+seqSetIndex_heavyChunk_collapsed)==2;
heavyChunk_seqCluster1 = collapsed_seqSet_inOne(heavyChunk_seqCluster1_index);
heavyChunk_seqCluster2 = collapsed_seqSet_inOne(heavyChunk_seqCluster2_index);

currentFigName = ['correlation_rProbandError_heavyChunk', '_'];
% to generate a unique file name for saving figure
fileName_fig16 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(16) == 1
        exportgraphics(fig16, fileName_fig16, 'BackgroundColor', 'none');
end

%% %%%%%%%%%%%%%%% ProbOffload & Accuracy, in one figure %%%%%%%%%%%%%%%%%
fig17 = figure(17);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point


for target_seqLength=1:pointKindsNum
    current_numSeq = length(target_seqSet{target_seqLength});

    subplot(pointKindsNum, 1, target_seqLength);
    plot(1:current_numSeq, offloadingProb_all{target_seqLength},'-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.75 0.25 0.25]);    
    hold on
    plot(1:current_numSeq, gError_noChoice_collapsed{target_seqLength},'-s', 'LineWidth', 2,'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.25 0.25 0.25]);
    hold on     
    
    if target_seqLength == 1
       legend('offloading rate', 'error rate', 'Location','northwest','fontsize',10) 
    end
    
    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
    
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength});%给坐标加标签
    ylabel(sprintf('Offloading rate \n or Error rate'), 'FontSize', 12, 'FontWeight', 'bold');    
    title(sprintf('Offloading rate & Error rate, length=%d, All-%ddays', target_seqLength, fileSize), 'FontSize',15);    
    % 获取xticklabel的忿
    xtl=get(gca,'XTickLabel');
    % 获取xtick的忿
    xt=get(gca,'XTick');
    % 获取ytick的忿
    yt=get(gca,'YTick');
    % 设置text的x坐标位置仿
    xtextp=xt;
    % 设置text的y坐标位置仿
    ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.08;
    %ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))*0.4-0.05;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    
    temp_fontsize = 13;%11-->13
    if numFrames == 8 && target_seqLength >= 3
        temp_fontsize = 9;
    end
    text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',temp_fontsize,'FontWeight','bold');
    % 取消原始ticklabel
    set(gca,'xticklabel','');    
    % 取消右、上边框
    set(gca,'box','off');
    
end

currentFigName = ['rProbANDacc', '_'];
% to generate a unique file name for saving figure
fileName_fig17 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(17) == 1
    exportgraphics(fig17, fileName_fig17, 'BackgroundColor', 'none');
end

%% %%%%%%%%%%%%%%%%%%%% RT %%%%%%%%%%%%%%%%%%%
% close all
fig18 = figure (18);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

offloadingProb;
selecting_RT;
selecting_RT_seqSet = offloadingProb;
selecting_RT_seqSet_all = offloadingProb_all;
selecting_RT_seqSet_collapsed = offloadingProb_all;

% ifSelectOffloading = selecting_RT;
ifSelectOffloading_seqSet = offloadingProb;
% ifSelectOffloading_seqSet_all = offloadingProb_all;
ifSelectOffloading_seqSet_collapsed = offloadingProb_all;

for file_index=1:fileSize
    temp_load = eval(['MAT_file_load.file',sprintf('%d',file_index)]);    
    for target_seqLength=1:pointKindsNum
        current_numSeq = length(target_seqSet{target_seqLength});        
        %offloading_trial_count_newSeq{file_index, target_seqLength} = zeros(current_numSeq, 1);
        selecting_RT_seqSet{file_index, target_seqLength} = cell(current_numSeq, 1); 
        ifSelectOffloading_seqSet{file_index, target_seqLength} = cell(current_numSeq, 1);
        for tempi=1:current_numSeq
            tempValid = 1;
            for tempj=1:length(target_seqSetIndex{file_index, target_seqLength}{tempi})                
                trial_count = target_seqSetIndex{file_index, target_seqLength}{tempi}(tempj);                
                if temp_load.trial_para.choiceCondition_flag(trial_count) == 2
                    selecting_RT_seqSet{file_index, target_seqLength}{tempi}(tempValid) = selecting_RT{file_index}(trial_count);
                    ifSelectOffloading_seqSet{file_index, target_seqLength}{tempi}(tempValid) =...
                        temp_load.trial_para.ifSelectOffloading(trial_count);                    
                    tempValid = tempValid + 1;
                end
            end
        end
        
    end
end

for file_index=1:fileSize
% for file_index=1:15
% for file_index=16:31
    if file_index == 1
        for target_seqLength=1:pointKindsNum
            current_numSeq = length(target_seqSet{target_seqLength});
            selecting_RT_seqSet_collapsed{target_seqLength} = cell(current_numSeq, 1);
            ifSelectOffloading_seqSet_collapsed{target_seqLength} = cell(current_numSeq, 1);
        end
    end
    for target_seqLength=1:pointKindsNum
        current_numSeq = length(target_seqSet{target_seqLength});
        for tempi=1:current_numSeq
            selecting_RT_seqSet_collapsed{target_seqLength}{tempi} = ...
                [ selecting_RT_seqSet_collapsed{target_seqLength}{tempi} ...
            selecting_RT_seqSet{file_index, target_seqLength}{tempi} ];   
        
            ifSelectOffloading_seqSet_collapsed{target_seqLength}{tempi} = ...
                [ ifSelectOffloading_seqSet_collapsed{target_seqLength}{tempi} ...
            ifSelectOffloading_seqSet{file_index, target_seqLength}{tempi} ];          
        end
    end
end

% get selecting_RT_seqSet_all
for target_seqLength=1:pointKindsNum
    current_numSeq = length(target_seqSet{target_seqLength});
    for tempi=1:current_numSeq
        selecting_RT_seqSet_all{target_seqLength}(tempi) = ...
            mean(selecting_RT_seqSet_collapsed{target_seqLength}{tempi});
    end
end 

% get selecting_RT_seqSet_all_g & selecting_RT_seqSet_all_r
selecting_RT_seqSet_all_g = selecting_RT_seqSet_all;
var_RT_seqSet_all_g = selecting_RT_seqSet_all;
selecting_RT_seqSet_all_r = selecting_RT_seqSet_all;
var_RT_seqSet_all_r = selecting_RT_seqSet_all;
for target_seqLength=1:pointKindsNum
    current_numSeq = length(target_seqSet{target_seqLength});
    for tempi=1:current_numSeq
        temp_g_index = ifSelectOffloading_seqSet_collapsed{target_seqLength}{tempi}==-1;
%         %use raw data
%         selecting_RT_seqSet_all_g{target_seqLength}(tempi) = ...
%             mean(selecting_RT_seqSet_collapsed{target_seqLength}{tempi}(temp_g_index));
%         var_RT_seqSet_all_g{target_seqLength}(tempi) = ...
%             var(selecting_RT_seqSet_collapsed{target_seqLength}{tempi}(temp_g_index)); 

        % wash data
        temp_RT_g = selecting_RT_seqSet_collapsed{target_seqLength}{tempi}(temp_g_index);
        temp_mean_g = mean(temp_RT_g);
        temp_var_g = var(temp_RT_g);
        temp_bool_g = temp_RT_g < (temp_mean_g-temp_var_g*3) | temp_RT_g > (temp_mean_g+temp_var_g*3);
        temp_RT_g_valid = temp_RT_g(~temp_bool_g);
        selecting_RT_seqSet_all_g{target_seqLength}(tempi) = mean(temp_RT_g_valid);
        var_RT_seqSet_all_g{target_seqLength}(tempi) = var(temp_RT_g_valid);
        
        temp_r_index = ifSelectOffloading_seqSet_collapsed{target_seqLength}{tempi}==1;
%         %use raw data        
%         selecting_RT_seqSet_all_r{target_seqLength}(tempi) = ...
%             mean(selecting_RT_seqSet_collapsed{target_seqLength}{tempi}(temp_r_index));
%         var_RT_seqSet_all_r{target_seqLength}(tempi) = ...
%             var(selecting_RT_seqSet_collapsed{target_seqLength}{tempi}(temp_r_index));
        
        % wash data
        temp_RT_r = selecting_RT_seqSet_collapsed{target_seqLength}{tempi}(temp_r_index);
        temp_mean_r = mean(temp_RT_r);
        temp_var_r = var(temp_RT_r);
        temp_bool_r = temp_RT_r < (temp_mean_r-temp_var_r*3) | temp_RT_r > (temp_mean_r+temp_var_r*3);
        temp_RT_r_valid = temp_RT_r(~temp_bool_r);
        selecting_RT_seqSet_all_r{target_seqLength}(tempi) = mean(temp_RT_r_valid);
        var_RT_seqSet_all_r{target_seqLength}(tempi) = var(temp_RT_r_valid);

        
    end        
end 

% ifSelectOffloading_seqSet_collapsed{target_seqLength}{tempi}

for target_seqLength=1:pointKindsNum
    current_numSeq = length(target_seqSet{target_seqLength});
       
    subplot(pointKindsNum, 1, target_seqLength);
    plot(1:current_numSeq, selecting_RT_seqSet_all{target_seqLength},'-s', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.25 0.25 0.75]);    
    hold on
    plot(1:current_numSeq, selecting_RT_seqSet_all_g{target_seqLength},'-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.25 0.75 0.25]);    
    hold on
    plot(1:current_numSeq, selecting_RT_seqSet_all_r{target_seqLength},'-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.75 0.25 0.25]);    
    hold on
    
    
    
    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
        
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength});%给坐标加标签
    ylabel('Reaction time (s)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('RT-SeqLength=%d, All-%ddays', target_seqLength, fileSize), 'FontSize',15);    
    % 获取xticklabel的忿
    xtl=get(gca,'XTickLabel');
    % 获取xtick的忿
    xt=get(gca,'XTick');
    % 获取ytick的忿
    yt=get(gca,'YTick');
    % 设置text的x坐标位置仿
    xtextp=xt;
    % 设置text的y坐标位置仿
    ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.08;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    text(xtextp,0.92*ones(1,current_numSeq),num2str(round(var_RT_seqSet_all_g{target_seqLength}*1000)),'HorizontalAlignment','center','rotation',0,'fontsize',13,'Color',[0.25 0.75 0.25]);    
    text(xtextp,0.4*ones(1,current_numSeq),num2str(round(var_RT_seqSet_all_r{target_seqLength}*1000)),'HorizontalAlignment','center','rotation',0,'fontsize',13,'Color',[0.75 0.25 0.25]);    
    text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',13,'FontWeight','bold');
    % 取消原始ticklabel
    set(gca,'xticklabel','');    
    % 取消右、上边框
    set(gca,'box','off');
    
end
if pointKindsNum == 4
    var_RT_g_inOne = [var_RT_seqSet_all_g{1}' var_RT_seqSet_all_g{2}' var_RT_seqSet_all_g{3}' var_RT_seqSet_all_g{4}'];
    var_RT_r_inOne = [var_RT_seqSet_all_r{1}' var_RT_seqSet_all_r{2}' var_RT_seqSet_all_r{3}' var_RT_seqSet_all_r{4}'];
elseif pointKindsNum == 3
    var_RT_g_inOne = [var_RT_seqSet_all_g{1}' var_RT_seqSet_all_g{2}' var_RT_seqSet_all_g{3}'];
    var_RT_r_inOne = [var_RT_seqSet_all_r{1}' var_RT_seqSet_all_r{2}' var_RT_seqSet_all_r{3}'];
end

var_RT_rMinusG_inOne = var_RT_r_inOne - var_RT_g_inOne;
var_RT_rMinusG = var_RT_seqSet_all_r;


for j=1:pointKindsNum  
    var_RT_rMinusG{j} = (var_RT_seqSet_all_r{j} - var_RT_seqSet_all_g{j})*1000;
    var_RT_rMinusG{j}(var_RT_rMinusG{j} < -50) = 0;
    var_RT_rMinusG{j}(isnan(var_RT_rMinusG{j})) = [];
    tempZero = zeros(1, length(var_RT_rMinusG{j}))';
    %paired-sample t-test
    [h_t(j),p_t(j)]=ttest(var_RT_rMinusG{j}, tempZero);  %#ok<*SAGROW>
    p(j) = p_t(j); 
    %Bonferroni Correction
%     if p(j) < 0.05/pointKindsNum(1)
%         tempTxt = sprintf('*');
%         if p(j) < 0.01/pointKindsNum(1)
%             tempTxt = sprintf('**');
%         elseif p(j) < 0.001/pointKindsNum(1)
%             tempTxt = sprintf('***');
%         end
% %         text(j+seq_length_rangeHead-1,mean(choiceMinusNoChoice(:,j))-0.08,tempTxt,'Color','black','FontSize',30,'FontWeight','bold',...
% %             'HorizontalAlignment','center');
%     end
end


currentFigName = ['RT_distri', '_'];
% to generate a unique file name for saving figure
fileName_fig18 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(18) == 1
        exportgraphics(fig18, fileName_fig18, 'BackgroundColor', 'none');
end

%% %%%%%%% correlation between rProb and error rate, low & high offloading rate %%%%%%%%%%%%
fig19 = figure(19);

offloadingProb_all_inOne;

temp1 = isnan(offloadingProb_all_inOne);
tempValidSeq_sorted = sort(offloadingProb_all_inOne(~temp1));
OffloadingRateThreshold_big = tempValidSeq_sorted(floor(length(tempValidSeq_sorted)*0.45))+0.0001;
offloadingRateThreshold_small = tempValidSeq_sorted(floor(length(tempValidSeq_sorted)*0.45))-0.0001;
tempResidual_big = abs(sum(tempValidSeq_sorted<OffloadingRateThreshold_big) - floor(length(tempValidSeq_sorted)/2));
tempResidual_small = abs(sum(tempValidSeq_sorted<offloadingRateThreshold_small) - floor(length(tempValidSeq_sorted)/2));
if tempResidual_big < tempResidual_small
    offloadingRateThreshold = OffloadingRateThreshold_big;
else
    offloadingRateThreshold = offloadingRateThreshold_small;
end

temp_color = [0.25 0.75 0.75;
            0.75 0.25 0.75;
            0.75 0.75 0.25;
            0.25 0.25 0.25];
for target_seqLength=1:pointKindsNum
    scatter(gError_noChoice_collapsed{target_seqLength}, offloadingProb_all{target_seqLength}, ...
        25, 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
end

hold on
plot([0 1], [offloadingRateThreshold offloadingRateThreshold], '-', 'LineWidth', 3, 'Color', [0.75 0.25 0.25]);
text(0.1,0.9,sprintf('45%%Threshold = %.4f',offloadingRateThreshold), 'fontsize',16,'FontWeight','bold', 'Color', [0.75 0.25 0.25]);

if pointKindsNum == 4
    legend('1 item','2 item','3 item','4 item','Threshold','Location','southeast','fontsize',8)
elseif pointKindsNum == 3
    legend('1 item','2 item','3 item','Threshold','Location','southeast','fontsize',8)
end

% text(0.1,0.9,sprintf('r2 = %.3f',r2), 'fontsize',16,'FontWeight','bold');

xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Error rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
title('Correlation, offloading rate & error rate', 'FontSize', 18, 'FontWeight', 'bold'); 


currentFigName = ['offloadingRate_threshold', '_'];
% to generate a unique file name for saving figure
fileName_fig19 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(19) == 1
        exportgraphics(fig19, fileName_fig19, 'BackgroundColor', 'none');
end


% offloadingRateThreshold = offloadingRateThreshold - 0.1;

fig20 = figure(20);

% for low & high offloading rate
seqSetIndex_lowOffloadingRate = cell(1, pointKindsNum);
seqSetIndex_highOffloadingRate = cell(1, pointKindsNum);
for target_seqLength=1:pointKindsNum
    seqSetIndex_lowOffloadingRate{target_seqLength} = offloadingProb_all{target_seqLength} < offloadingRateThreshold;
    seqSetIndex_highOffloadingRate{target_seqLength} = offloadingProb_all{target_seqLength} >= offloadingRateThreshold;
end
seqSetIndex_lowOffloadingRate_collapsed = offloadingProb_all_inOne < offloadingRateThreshold;
seqSetIndex_highOffloadingRate_collapsed = offloadingProb_all_inOne >= offloadingRateThreshold;
  

x = gError_noChoice_collapsed_inOne(seqSetIndex_lowOffloadingRate_collapsed);
y = offloadingProb_all_inOne(seqSetIndex_lowOffloadingRate_collapsed);
n = 1;
[p,S] = polyfit(x,y,n);
r2 = 1 - (S.normr/norm(y - mean(y)))^2;

x_fit_low = 0:0.1:1;
y_fit_low = polyval(p,x_fit_low);


temp_color = [0.25 0.75 0.75;
            0.75 0.25 0.75;
            0.75 0.75 0.25;
            0.25 0.25 0.25];
        
for target_seqLength=1:pointKindsNum
    scatter(gError_noChoice_collapsed{target_seqLength}(seqSetIndex_lowOffloadingRate{target_seqLength})...
        , offloadingProb_all{target_seqLength}(seqSetIndex_lowOffloadingRate{target_seqLength})...
        , ...
        25, 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
end
               
plot(x_fit_low, y_fit_low, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);

if pointKindsNum == 4
    legend('1 item','2 item','3 item','4 item','fit','Location','southeast','fontsize',8)
elseif pointKindsNum == 3
    legend('1 item','2 item','3 item','fit','Location','southeast','fontsize',8)
end

text(0.1,0.9,sprintf('r2 = %.3f',r2), 'fontsize',16,'FontWeight','bold');



xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Error rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
title('Correlation, low Offloading rate', 'FontSize', 20, 'FontWeight', 'bold'); 

currentFigName = ['correlation_rProbandError_lowOffloadingRate', '_'];
% to generate a unique file name for saving figure
fileName_fig20 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(20) == 1
        exportgraphics(fig20, fileName_fig20, 'BackgroundColor', 'none');
end


fig21 = figure(21);

seqSetIndex_highOffloadingRate_collapsed;

x = gError_noChoice_collapsed_inOne(seqSetIndex_highOffloadingRate_collapsed);
y = offloadingProb_all_inOne(seqSetIndex_highOffloadingRate_collapsed);
n = 1;
[p,S] = polyfit(x,y,n);
r2 = 1 - (S.normr/norm(y - mean(y)))^2;

x_fit_low = 0:0.1:1;
y_fit_low = polyval(p,x_fit_low);


temp_color = [0.25 0.75 0.75;
            0.75 0.25 0.75;
            0.75 0.75 0.25;
            0.25 0.25 0.25];
        
for target_seqLength=1:pointKindsNum
    scatter(gError_noChoice_collapsed{target_seqLength}(seqSetIndex_highOffloadingRate{target_seqLength})...
        , offloadingProb_all{target_seqLength}(seqSetIndex_highOffloadingRate{target_seqLength})...
        , ...
        25, 'filled', 'MarkerFaceColor', temp_color(target_seqLength, :), ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
end
               
plot(x_fit_low, y_fit_low, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);

if pointKindsNum == 4
    legend('1 item','2 item','3 item','4 item','fit','Location','southeast','fontsize',8)
elseif pointKindsNum == 3
    legend('1 item','2 item','3 item','fit','Location','southeast','fontsize',8)
end
text(0.1,0.9,sprintf('r2 = %.3f',r2), 'fontsize',16,'FontWeight','bold');



xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Error rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
title('Correlation, high Offloading rate', 'FontSize', 20, 'FontWeight', 'bold'); 


currentFigName = ['correlation_rProbandError_highOffloadingRate', '_'];
% to generate a unique file name for saving figure
fileName_fig21 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(21) == 1
        exportgraphics(fig21, fileName_fig21, 'BackgroundColor', 'none');
end

%% %%%%%%%%%%%%%%% choice and noChoice gAcc, in one figure %%%%%%%%%%%%%%%%%
fig22 = figure(22);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point


for target_seqLength=1:pointKindsNum
    current_numSeq = length(target_seqSet{target_seqLength});

    subplot(pointKindsNum, 1, target_seqLength);
    plot(1:current_numSeq, gAcc_choice_collapsed{target_seqLength},'-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.25 0.75 0.25]);    
    hold on
    plot(1:current_numSeq, gAcc_noChoice_collapsed{target_seqLength},'-s', 'LineWidth', 2,'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.25 0.25 0.25]);
    hold on     
    
    if target_seqLength == 1
       legend('choice accuracy', 'no choice accuracy', 'Location','southwest','fontsize',10) 
    end
    
    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
    
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength});%给坐标加标签
    ylabel(sprintf('Accuracy'), 'FontSize', 12, 'FontWeight', 'bold');    
    title(sprintf('Choice and no choice accuracy, length=%d, All-%ddays', target_seqLength, fileSize), 'FontSize',15);    
    % 获取xticklabel的忿
    xtl=get(gca,'XTickLabel');
    % 获取xtick的忿
    xt=get(gca,'XTick');
    % 获取ytick的忿
    yt=get(gca,'YTick');
    % 设置text的x坐标位置仿
    xtextp=xt;
    % 设置text的y坐标位置仿
    ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.08;
    %ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))*0.4-0.05;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    temp_fontsize = 13;%11-->13
    if numFrames == 8 && target_seqLength >= 3
        temp_fontsize = 9;
    end    
    text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',temp_fontsize,'FontWeight','bold');
    % 取消原始ticklabel
    set(gca,'xticklabel','');    
    % 取消右、上边框
    set(gca,'box','off');
    
end

currentFigName = ['choiceANDnoChoice', '_'];
% to generate a unique file name for saving figure
fileName_fig22 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(22) == 1
    exportgraphics(fig22, fileName_fig22, 'BackgroundColor', 'none');
end


%% %%%%%%%%%%%%%%%%%%%% gAcc_noChoice_inOrder %%%%%%%%%%%%%%%%%%%
fig23 = figure(23);
set(gcf,'Position',[0 50 700 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

gAcc_noChoice_inOrder_collapsed;
temp_color = [0.00 0.00 0.00;
              0.28 0.28 0.28;
              0.56 0.56 0.56;
              0.84 0.84 0.84];
temp_mean = cell(1, pointKindsNum);
temp_SEM = cell(1, pointKindsNum);
temp_mean_collaped = [];
%for target_seqLength=1:pointKindsNum
for target_seqLength=pointKindsNum:-1:1
    %gAcc_noChoice_inOrder_collapsed{target_seqLength};
    temp_mean{target_seqLength} = mean(gAcc_noChoice_inOrder_collapsed{target_seqLength}, 1);
    temp_SEM{target_seqLength} = std(gAcc_noChoice_inOrder_collapsed{target_seqLength}, 1)...
        ./sqrt( size(gAcc_noChoice_inOrder_collapsed{target_seqLength}, 1) );
    errorbar(seq_length_rangeHead:target_seqLength,temp_mean{target_seqLength},temp_SEM{target_seqLength},...
        '-o','Color',temp_color(target_seqLength,:),'LineWidth',1,'CapSize',12,'MarkerSize',5);
    hold on
    
    temp_mean_collaped = [temp_mean_collaped temp_mean{target_seqLength}];
end

legend('length=4','length=3','length=2','length=1',...
    'Location','southwest','fontsize',13)


% plot([seq_length_rangeHead-1 seq_length_rangeHead: seq_length_rangeTail seq_length_rangeTail+1]...
%     , zeros(1, pointKindsNum+2), '--');

ylim([min(temp_mean_collaped)-0.1 1]);
set(gca, 'FontSize', 20)
% set(gca,'YTick',[-0.3:0.1:0.3]);%设置要显示坐标刻度的范围
set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
set(gca,'box','off');% 取消右、上边框
xlabel('Order', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
title({'\fontsize{20}{\bf Accuracy of no choice condition, in order}'});

currentFigName = ['gAcc_noChoice_inOrder', '_'];
% to generate a unique file name for saving figure
fileName_fig23 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(23) == 1
        exportgraphics(fig23, fileName_fig23, 'BackgroundColor', 'none');
end
%% %%%%%%%%%%%%%%%%%%%% gDistri_length %%%%%%%%%%%%%%%%%%%
fig24 = figure(24);
set(gcf,'Position',[0 50 700 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

gDistri_length_noChoice_collapsed_seqAVG;
gDistri_length_noChoice_collapsed_std;
gDistri_length_noChoice_collapsed;
length_rangeHead = 0;
length_rangeTail = numFrames;% 6

temp_color = [0.00 0.00 0.00;
              0.28 0.28 0.28;
              0.56 0.56 0.56;
              0.84 0.84 0.84];
temp_mean = cell(1, pointKindsNum);
temp_std = cell(1, pointKindsNum);
for target_seqLength=1:pointKindsNum
    temp_mean{target_seqLength} = gDistri_length_noChoice_collapsed_seqAVG{target_seqLength};
    temp_std{target_seqLength} = gDistri_length_noChoice_collapsed_std{target_seqLength};    
    plot(length_rangeHead:length_rangeTail,temp_mean{target_seqLength},...
        '-o','Color',temp_color(target_seqLength,:),'LineWidth',1,'MarkerSize',5);
%     errorbar(length_rangeHead:length_rangeTail,temp_mean{target_seqLength},temp_std{target_seqLength},...
%         '-o','Color',temp_color(target_seqLength,:),'LineWidth',1,'CapSize',12,'MarkerSize',5);
    hold on    
end

legend('length=1','length=2','length=3','length=4',...
    'Location','northeast','fontsize',13)

% ylim([min(temp_mean_collaped)-0.1 1]);
set(gca, 'FontSize', 20)
% set(gca,'YTick',[-0.3:0.1:0.3]);%设置要显示坐标刻度的范围
set(gca,'XLim',[length_rangeHead-1 length_rangeTail+1]);%X轴的数据显示范围
set(gca,'XTick',[length_rangeHead:1:length_rangeTail]);%设置要显示坐标刻度的范围
set(gca,'XTickLabel',[length_rangeHead:1:length_rangeTail]);%给坐标加标签
set(gca,'box','off');% 取消右、上边框
xlabel('Response length', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 18, 'FontWeight', 'bold');
title({'\fontsize{20}{\bf Distribution of response length}'});

%% %%%%%%%%%%%%%%%%%%%% gDistri_neighborLengthError_ismember %%%%%%%%%%%%%%%%%%%
fig25 = figure(25);
set(gcf,'Position',[0 50 700 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

gDistri_neighborLengthError_ismember_noChoice_collapsed;
length_rangeHead = seq_length_rangeHead;
length_rangeTail = seq_length_rangeTail;

temp_color = [0.0 0.0 0.0;
              0.4 0.4 0.4;
              0.8 0.8 0.8];

temp_mean = zeros(pointKindsNum, 3);
temp_SEM = zeros(pointKindsNum, 3);
for tempi=1:3
    for target_seqLength=1:pointKindsNum
        temp_valid = ~isnan(gDistri_neighborLengthError_ismember_noChoice_collapsed{target_seqLength}(:, tempi));
        
        temp_mean(target_seqLength, tempi) = mean(gDistri_neighborLengthError_ismember_noChoice_collapsed{target_seqLength}(temp_valid, tempi));
        temp_SEM(target_seqLength, tempi) = std(gDistri_neighborLengthError_ismember_noChoice_collapsed{target_seqLength}(temp_valid, tempi))...
            / sqrt(numSeq(target_seqLength));        
    end      
    % plot(length_rangeHead:length_rangeTail,temp_mean,...
    %     '-o','Color', temp_color(tempi, :)'LineWidth',1,'MarkerSize',5);
    errorbar(length_rangeHead:length_rangeTail,temp_mean(:, tempi),temp_SEM(:, tempi),...
        '-o','Color', temp_color(tempi, :),'LineWidth',1,'CapSize',12,'MarkerSize',5);
    hold on
end
legend('down and up length','down length','up length',...
    'Location','southeast','fontsize',13);


ylim([0 1]);
set(gca, 'FontSize', 20)
% set(gca,'YTick',[-0.3:0.1:0.3]);%设置要显示坐标刻度的范围
set(gca,'XLim',[length_rangeHead-1 length_rangeTail+1]);%X轴的数据显示范围
set(gca,'XTick',[length_rangeHead:1:length_rangeTail]);%设置要显示坐标刻度的范围
set(gca,'XTickLabel',[length_rangeHead:1:length_rangeTail]);%给坐标加标签
set(gca,'box','off');% 取消右、上边框
xlabel('Length', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 18, 'FontWeight', 'bold');
title({'\fontsize{20}{\bf Neighbor length error, delet one or add one}'});

%% Get gAcc_noChoice_2, gAcc_choice_2
gAcc_noChoice;
gAcc_choice;
gAcc_noChoice_2 = cell(1, pointKindsNum);
gAcc_choice_2 = cell(1, pointKindsNum);
isoutlier_noChoice = cell(1, pointKindsNum);
isoutlier_choice = cell(1, pointKindsNum);
for tempi=1:pointKindsNum
    gAcc_noChoice_2{tempi} = zeros(fileSize, numSeq(tempi));
    gAcc_choice_2{tempi} = zeros(fileSize, numSeq(tempi));
    for tempj=1:numSeq(tempi)
        for tempk=1:fileSize
            gAcc_noChoice_2{tempi}(tempk, tempj) = gAcc_noChoice{tempk, tempi}(tempj);
            gAcc_choice_2{tempi}(tempk, tempj) = gAcc_choice{tempk, tempi}(tempj);
        end
    end
    
    isoutlier_noChoice{tempi} = isoutlier(gAcc_noChoice_2{tempi}, 'quartiles');
    isoutlier_choice{tempi} = isoutlier(gAcc_choice_2{tempi}, 'quartiles');
end
gAcc_noChoice_2;
gAcc_choice_2;


%% Plot gAcc_noChoice with violinplot
fig26 = figure(26);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

for target_seqLength=1:pointKindsNum

    current_numSeq = length(target_seqSet{target_seqLength});

    subplot(pointKindsNum, 1, target_seqLength);
    
    
    %boxplot(gAcc_noChoice_2{target_seqLength});
    temp_color1 = [0.25 0.25 0.25];
    temp_color2 = repmat(temp_color1, numSeq(target_seqLength), 1);
    temp_plot = gAcc_noChoice_2{target_seqLength};
    for tempi=1:numSeq(target_seqLength)
        if sum(isnan(temp_plot(:, tempi))) == size(temp_plot, 1)
            temp_plot(1, tempi) = 0;
            
        end
    end
    violinplot(temp_plot,[],'ViolinColor',temp_color2);         
    %violinplot(gAcc_noChoice_2{target_seqLength},[],'ViolinColor',temp_color2);
    
    hold on 
    plot(0:current_numSeq+1, 0.05*ones(1,current_numSeq+2),'--','Color',[0.75,0.25,0.25]);

    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
        
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength});%给坐标加标签
    ylabel('Accuracy', 'FontSize', 12, 'FontWeight', 'bold');    
    title(sprintf('Accuracy-noChoice-SeqLength=%d, All-%ddays', target_seqLength, fileSize), 'FontSize',15);    

    % 获取xticklabel的忿
    xtl=get(gca,'XTickLabel');
    % 获取xtick的忿
    xt=get(gca,'XTick');
    % 获取ytick的忿
    yt=get(gca,'YTick');
    % 设置text的x坐标位置仿
    xtextp=xt;
    % 设置text的y坐标位置仿
    ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.08;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',13,'FontWeight','bold');
    % 取消原始ticklabel
    set(gca,'xticklabel','');   
    % 取消右、上边框
    set(gca,'box','off');    
end
currentFigName = ['gAcc_noChoice_violinplot', '_'];
% to generate a unique file name for saving figure
fileName_fig26 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(26) == 1
        exportgraphics(fig26, fileName_fig26, 'BackgroundColor', 'none');
end

%% Plot gAcc_choice with violinplot
fig27 = figure(27);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

for target_seqLength=1:pointKindsNum

    current_numSeq = length(target_seqSet{target_seqLength});

    subplot(pointKindsNum, 1, target_seqLength);
    
    
    %boxplot(gAcc_noChoice_2{target_seqLength});
    temp_color1 = [0.25 0.75 0.25];
    temp_color2 = repmat(temp_color1, numSeq(target_seqLength), 1);
    temp_plot = gAcc_choice_2{target_seqLength};
    for tempi=1:numSeq(target_seqLength)
        if sum(isnan(temp_plot(:, tempi))) == size(temp_plot, 1)
            temp_plot(1, tempi) = 0;
            
        end
    end
    violinplot(temp_plot,[],'ViolinColor',temp_color2); 
        
    hold on 
    plot(0:current_numSeq+1, 0.05*ones(1,current_numSeq+2),'--','Color',[0.75,0.25,0.25]);
    
    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
        
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength});%给坐标加标签
    ylabel('Accuracy', 'FontSize', 12, 'FontWeight', 'bold');    
    title(sprintf('Accuracy-choice-SeqLength=%d, All-%ddays', target_seqLength, fileSize), 'FontSize',15);    

    % 获取xticklabel的忿
    xtl=get(gca,'XTickLabel');
    % 获取xtick的忿
    xt=get(gca,'XTick');
    % 获取ytick的忿
    yt=get(gca,'YTick');
    % 设置text的x坐标位置仿
    xtextp=xt;
    % 设置text的y坐标位置仿
    ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.08;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    text(xtextp,ytextp,xtl,'HorizontalAlignment','center','rotation',45,'fontsize',13,'FontWeight','bold');
    % 取消原始ticklabel
    set(gca,'xticklabel','');   
    % 取消右、上边框
    set(gca,'box','off');    
end
currentFigName = ['gAcc_choice_violinplot', '_'];
% to generate a unique file name for saving figure
fileName_fig27 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(27) == 1
        exportgraphics(fig27, fileName_fig27, 'BackgroundColor', 'none');
end

%% Selection RT
fig92 = figure(92);

temp_RT_offset = -0.345;
% temp_RT_offset = 0;

selecting_RT_seqSet_all;

selecting_RT_seqSet_all_inOne = [];
for tempi=1:4
    selecting_RT_seqSet_all_inOne = [selecting_RT_seqSet_all_inOne selecting_RT_seqSet_all{tempi}'];
end

selecting_RT_seqSet_all_inOne;

selecting_RT_seqSet_all_inOne_adjusted = selecting_RT_seqSet_all_inOne + temp_RT_offset;

offloadingProb_all_inOne;

x = offloadingProb_all_inOne;
y = selecting_RT_seqSet_all_inOne_adjusted;

% [rA,pA] = corr(offloadingProb_all_inOne',selecting_RT_seqSet_all_inOne_adjusted');

% [rB,pB] = corr(abs(offloadingProb_all_inOne-0.5)',selecting_RT_seqSet_all_inOne_adjusted');

% temp_coeffA = polyfit(offloadingProb_all_inOne,selecting_RT_seqSet_all_inOne_adjusted,1);

[ymin,ymax] = bounds(y);

[p_mapping,S] = polyfit(x,y,n);
r2 = 1 - (S.normr/norm(y - mean(y)))^2;
[r, p_corr] = corr(x',y');
tempTxt = 'p>0.05';
if p_corr < 0.05    
    tempTxt = 'p<0.05';
end
if p_corr < 0.01
    tempTxt = 'p<0.01';
end
if p_corr < 0.001
    tempTxt = 'p<0.001';
end

x_fit = 0:0.1:1;
y_fit = polyval(p_mapping,x_fit);

temp_color = [[166,97,26]/255;
            [223,194,125]/255;
            [128,205,193]/255;
            [1,133,113]/255];

temp_faceAlpha = 0.5;

for target_seqLength=1:pointKindsNum    
    temp_size = ((target_seqLength.^3)*2 + 3) .* ones(1, length(offloadingProb_all{target_seqLength}));    
    scatter(offloadingProb_all{target_seqLength}, selecting_RT_seqSet_all{target_seqLength}+temp_RT_offset, ...
        temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', temp_faceAlpha, 'MarkerEdgeAlpha', 0.7);       
    hold on
end
plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
% legend('length 1','length 2','length 3','length 4','fit','Location','southeast','fontsize',12)
text(0.1,ymax+(ymax-ymin)*0.01,sprintf('r=%.3f, %s',r,tempTxt), 'fontsize',14,'FontWeight','bold');

xlim([0 1]);
% ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
% set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');  
ylabel('Selection RT (s)', 'FontSize', 18, 'FontWeight', 'bold');  
title('Correlation, Selection RT & offloading rate', 'FontSize', 18, 'FontWeight', 'bold'); 


%% End
cd 'C:\ASDROOT\STUDY\MonkeyDataAnalysis'
