% This script: To preprocess behavioral data
%% Initialization
clear all %#ok<CLALL>
close all

%% Some flags to control the following script
if_monkey_D0_Z1 = 1;% To decide whether dealing with Ding's data or Zelku's data
if_ruleOutRToutlier = 1;
if_isCorrect_modified = 1;
if_ruleOut_someResponse = 1;
if_ruleOut_errorStopTrial = 1;%1
if_ruleOut_outlierSeqSession = 0;
% if_backTouch_modified = 1;

path_results = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\data';
if if_monkey_D0_Z1 == 0
    % searchName = 'DEYECopy1';
    %searchName = 'DnumFrames6_A';
    %searchName = 'DnumFrames8_A';
    %searchName = 'DnumFrames7_A';
    %searchName = 'DnumFrames6_Train1';
    %searchName = 'DnumFrames6_Train2';
    searchName = 'DnumFrames6_2p';
    if if_ruleOut_outlierSeqSession == 1
        searchName_isoutlier = 'D_isoutlier_1';
        searchName_mappingParam = 'from21-12-06to21-12-13_D_mappingParam_1';
    end
elseif if_monkey_D0_Z1 == 1
    %searchName = 'ZnumFrames6_A';
    %searchName = 'ZnumFrames6_B';
    %searchName = 'ZnumFrames6_B2';
    %searchName = 'ZnumFrames7_A';
    %searchName = 'ZnumFrames6_Train1';
    %searchName = 'ZnumFrames6_Train2';
    %searchName = 'ZnumFrames6_Train3';
    %searchName = 'ZnumFrames6_C';
    %searchName = 'ZnumFrames6_D';
    %searchName = 'ZnumFrames6_early2p';
    %searchName = 'ZnumFrames6_all2p';
    searchName = 'ZnumFrames6_almost2p';
    if if_ruleOut_outlierSeqSession == 1
        searchName_isoutlier = 'Z_isoutlier_2';
        searchName_mappingParam = 'from21-06-16to21-08-04_Z_mappingParam_1.mat';
    end
end

if if_ruleOut_outlierSeqSession == 1
    load_isoutlierSeqSession = loadMat_single(searchName_isoutlier, path_results);
    isoutlier_noChoice = load_isoutlierSeqSession.isoutlier_noChoice;
    isoutlier_choice = load_isoutlierSeqSession.isoutlier_choice;
    
    isoutlier_noChoice_inOne = [];
    isoutlier_choice_inOne = [];
    for tempi=1:length(isoutlier_noChoice)
        isoutlier_noChoice_inOne = [isoutlier_noChoice_inOne isoutlier_noChoice{tempi}]; %#ok<*AGROW>
        isoutlier_choice_inOne = [isoutlier_choice_inOne isoutlier_choice{tempi}];      
    end
    
    % To load some existing parameters
    load_mapping = loadMat_single(searchName_mappingParam, path_results);
    target_seqSet = load_mapping.target_seqSet;
    pointKindsNum = length(target_seqSet);
    target_seqSet_inOne = [];
    for tempi=1:pointKindsNum
        target_seqSet_inOne = [target_seqSet_inOne; target_seqSet{tempi}];     %#ok<*AGROW>
    end
    numSeq = zeros(1, pointKindsNum);
    for target_seqLength=1:pointKindsNum
        numSeq(target_seqLength) = size(target_seqSet{target_seqLength}, 1);
    end
end

% targetPATH = 'C:\ASDROOT\STUDY\backups\MatlabScript_226MIDDLE_20211014\MATDATA';
% targetPATH = 'C:\ASDROOT\STUDY\backups\MATDATA_20220630';
targetPATH = 'C:\ASDROOT\STUDY\backups\MATDATA_master';
targetPATH = [targetPATH, '\', searchName];

cd(targetPATH)


MAT_file=dir([searchName, '*']);
MAT_file_name = cell(1, max(size(MAT_file)));
MAT_file_load_cell = cell(1, max(size(MAT_file)));
rowHeadings = cell(1, max(size(MAT_file)));
for i=1:size(MAT_file)
    MAT_file_name(i) = {MAT_file(i).name};
    MAT_file_load_cell(i) = {load(cell2mat(MAT_file_name(i)))};    
    rowHeadings(i) = {sprintf('file%d', i)};   
end
MAT_file_load = cell2struct(MAT_file_load_cell,rowHeadings,2);
subject_name = cell(1, max(size(MAT_file)));



initialization = 0;
for file_index=1:size(MAT_file)
    temp_load = eval(['MAT_file_load.file',sprintf('%d',file_index)]);    
    
    seq_length_rangeHead = temp_load.basic_para.seq_length_rangeHead;
    seq_length_rangeTail = temp_load.basic_para.seq_length_rangeTail;           
    pointKindsNum = seq_length_rangeTail-seq_length_rangeHead+1;
    
    
    if initialization == 0
        
        RT_outlierTrial = cell(max(size(MAT_file)), pointKindsNum);
        sumRT_outlierTrial = zeros(max(size(MAT_file)), pointKindsNum);
        RT_outlierTrial_collapsed = cell(max(size(MAT_file)), 1);
        
        selecting_RT = cell(max(size(MAT_file)), pointKindsNum);
        meanSelecting_RT = zeros(max(size(MAT_file)), pointKindsNum);
        stdSelecting_RT = zeros(max(size(MAT_file)), pointKindsNum);
        
        response_outlierTrial = cell(max(size(MAT_file)), 1);
        errorStop_trial = cell(max(size(MAT_file)), 1);
        xSeqxSession_outlierTrial = cell(max(size(MAT_file)), 1);
        all_outlierTrial = cell(max(size(MAT_file)), 1);        
        
        initialization = 1;
    end
    
    
    trial_num = temp_load.trial_para.trial_count;
    cumulativeError = temp_load.trial_para.cumulativeError;
    
    if if_isCorrect_modified == 1
        for tempi=1:trial_num
            %temp_load.trial_para.isCorrect(tempi);
            %temp_load.trial_para.currentSequence(tempi);
            %temp_load.trial_para.isFilled(tempi);
            current_seq = temp_load.trial_para.currentSequence{tempi};
            current_isFilled = temp_load.trial_para.isFilled{tempi};
            current_trial_response = find(current_isFilled==0);
            if length(current_trial_response) == length(current_seq) 
                if temp_load.trial_para.isCorrect(tempi) == -1 && ...
                        sum(current_trial_response == current_seq) == length(current_seq)
                    temp_load.trial_para.isCorrect(tempi) = 1;
                end
            end
        end
    end
    
    
    isCorrect = temp_load.trial_para.isCorrect;

    
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
    
%     randArray2 = temp_load.basic_para.randArray2;
%     choiceCondition_boardLine1 = temp_load.basic_para.choiceCondition_boardLine1;
%     choiceCondition_boardLine2 = temp_load.basic_para.choiceCondition_boardLine2;
%     
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
    
    
    RT_outlierTrial_collapsed{file_index} = false(1, trial_num);
    for j=1:pointKindsNum
        tempIndex = find(temp_load.trial_para.seq_length(1: trial_num) == (j+seq_length_rangeHead-1));
                                        
        tempIndex_allChoice = tempIndex(choiceCondition_flag(tempIndex) == 2);
                
        tempIndex_allChoice_newSeq = tempIndex_allChoice(selectFlag_newSeq(tempIndex_allChoice) == 1);

        if isempty(tempIndex_allChoice_newSeq) == 0
            if tempIndex_allChoice_newSeq(1) == 1
               tempIndex_allChoice_newSeq = tempIndex_allChoice_newSeq(2:end); 
            end
        end
        
                     
        selecting_RT{file_index, j} = temp_load.trial_para.selecting_RT(tempIndex_allChoice_newSeq);
        meanSelecting_RT(file_index, j) = sum(selecting_RT{file_index, j})/length(selecting_RT{file_index, j});
        stdSelecting_RT(file_index, j) = std(selecting_RT{file_index, j});
        RT_outlierTrial{file_index, j} = false(1, length(selecting_RT{file_index, j}));
        if if_ruleOutRToutlier == 1
            for tempi=1:length(selecting_RT{file_index, j})
                if abs(selecting_RT{file_index, j}(tempi)-meanSelecting_RT(file_index, j)) > 3*stdSelecting_RT(file_index, j)
                    RT_outlierTrial{file_index, j}(tempi) = true;
                end
            end
        end
        sumRT_outlierTrial(file_index, j) = sum(RT_outlierTrial{file_index, j});
        
        RT_outlierTrial_collapsed{file_index}(tempIndex_allChoice_newSeq(RT_outlierTrial{file_index, j})) = true;               
        
    end
    %tempa(file_index)=sum(RT_outlierTrial_collapsed{file_index}) - sum(sumRT_outlierTrial(file_index, :));

    RT_outlierTrial_collapsed{file_index};
    
    % rule out specific response, isFilled are zeros or ones or
    % response_length >= 6
    response_outlierTrial{file_index} = false(1, trial_num);
    numFrames = temp_load.basic_para.numSquares;
    if if_ruleOut_someResponse == 1
        for trial_count=1:trial_num
            trial_count;
            current_isFilled = temp_load.trial_para.isFilled{trial_count};
            if sum(current_isFilled) == 0 || sum(current_isFilled) == numFrames
                response_outlierTrial{file_index}(trial_count) = true;
            end
            current_trial_response = find(current_isFilled==0);
            if length(current_trial_response) >= 6
                response_outlierTrial{file_index}(trial_count) = true;
            end
        end
    end
    
    % rule out error stop trials, only use free touch trials.
    errorStop_trial{file_index} = false(1, trial_num);
    ifFreeTouch = temp_load.basic_para.ifFreeTouch(1:trial_num);
    if if_ruleOut_errorStopTrial == 1        
        errorStop_trial{file_index} = ifFreeTouch==0;
    end
    
    % rule out specific seqs in specific sessions, which are outliers (95% CI)
    ifSelectOffloading = temp_load.trial_para.ifSelectOffloading;
    xSeqxSession_outlierTrial{file_index} = false(1, trial_num);
    if if_ruleOut_outlierSeqSession == 1
        temp_isoutlier_noChoice = isoutlier_noChoice_inOne(file_index, :);
        temp_isoutlier_choice = isoutlier_choice_inOne(file_index, :);
        target_seqSet_inOne;
        
        for trial_count=1:trial_num
            temp_seq = temp_load.trial_para.currentSequence{trial_count};
            temp_length = length(temp_seq);
            temp_seqSet_index = 0;
            for tempi=1:sum(numSeq)
                if temp_length == length(target_seqSet_inOne{tempi})
                    if sum(ismember(temp_seq, target_seqSet_inOne{tempi})) == temp_length
                        temp_seqSet_index = tempi;
                        break
                    end
                end
            end
            temp_seqSet_index;
            choiceCondition_flag(trial_count); %#ok<*NOEFF>
            if choiceCondition_flag(trial_count) == 0
                % no choice
                if temp_isoutlier_noChoice(temp_seqSet_index) == 1
                    xSeqxSession_outlierTrial{file_index}(trial_count) = 1;                    
                end
                
%             elseif choiceCondition_flag(trial_count) == 2 && ...
%                     ifSelectOffloading(trial_count) == -1
            elseif choiceCondition_flag(trial_count) == 2
                % choice
                if temp_isoutlier_choice(temp_seqSet_index) == 1
                    xSeqxSession_outlierTrial{file_index}(trial_count) = 1;                    
                end            
                
            end
        end
    end
    
    all_outlierTrial{file_index} = RT_outlierTrial_collapsed{file_index} | ...
        response_outlierTrial{file_index} | errorStop_trial{file_index} | ...
        xSeqxSession_outlierTrial{file_index};
    

    
    % basic_para
    temp_load.basic_para.randArray1 = temp_load.basic_para.randArray1(all_outlierTrial{file_index}==0);
    temp_load.basic_para.randArray2 = temp_load.basic_para.randArray2(all_outlierTrial{file_index}==0);
    temp_load.basic_para.randArray3 = temp_load.basic_para.randArray3(all_outlierTrial{file_index}==0);
    temp_load.basic_para.randArray4 = temp_load.basic_para.randArray4(all_outlierTrial{file_index}==0);
    temp_load.basic_para.ifFreeTouch = temp_load.basic_para.ifFreeTouch(all_outlierTrial{file_index}==0);
    temp_load.basic_para.fixationTime = temp_load.basic_para.fixationTime(all_outlierTrial{file_index}==0);
    temp_load.basic_para.twoFixationTime = temp_load.basic_para.twoFixationTime(all_outlierTrial{file_index}==0);
    temp_load.basic_para.cumulativeErrorLimit = temp_load.basic_para.cumulativeErrorLimit(all_outlierTrial{file_index}==0);
    temp_load.basic_para.if_outlierTrial = all_outlierTrial{file_index};
    
    % trial_para
    raw_trial_count = temp_load.trial_para.trial_count;
    temp_load.trial_para.trial_count = sum(all_outlierTrial{file_index}==0);
    temp_load.trial_para.uniqueSeq_trial_count = temp_load.trial_para.uniqueSeq_trial_count - ...
        (raw_trial_count - temp_load.trial_para.trial_count);
    temp_load.trial_para.isGreenRectLeft = temp_load.trial_para.isGreenRectLeft(all_outlierTrial{file_index}==0);
    temp_load.trial_para.seq_length = temp_load.trial_para.seq_length(all_outlierTrial{file_index}==0);
    temp_load.trial_para.numFrames = temp_load.trial_para.numFrames(all_outlierTrial{file_index}==0);
    
    if isfield(temp_load.trial_para,'isWaitGocueBreak') == 1    
        temp_min_length = min(length(temp_load.trial_para.isWaitGocueBreak),length(all_outlierTrial{file_index}));
        temp1 = all_outlierTrial{file_index}(1:temp_min_length);
        temp_load.trial_para.isWaitGocueBreak = temp_load.trial_para.isWaitGocueBreak(temp1==0);
    else
        temp_load.trial_para.isWaitGocueBreak = zeros(1, sum(all_outlierTrial{file_index}==0));
    end
    
    temp_load.trial_para.selecting_RT = temp_load.trial_para.selecting_RT(all_outlierTrial{file_index}==0);
    temp_load.trial_para.selecting_point = temp_load.trial_para.selecting_point(all_outlierTrial{file_index}==0);
    temp_load.trial_para.ifSelectOffloading = temp_load.trial_para.ifSelectOffloading(all_outlierTrial{file_index}==0);
    temp_load.trial_para.endingHold_RT = temp_load.trial_para.endingHold_RT(...
        all_outlierTrial{file_index}(1:length(temp_load.trial_para.endingHold_RT))==0);
    
    temp_load.trial_para.isCorrect = temp_load.trial_para.isCorrect(all_outlierTrial{file_index}==0);
    temp_load.trial_para.isIppatsuCorrect = temp_load.trial_para.isIppatsuCorrect(all_outlierTrial{file_index}==0);
    temp_load.trial_para.isTouchFrameCorrect = temp_load.trial_para.isTouchFrameCorrect(all_outlierTrial{file_index}==0);
    temp_load.trial_para.isTouchingTimeOut = temp_load.trial_para.isTouchingTimeOut(all_outlierTrial{file_index}==0);
    temp_load.trial_para.cumulativeCorrect_greenRedSwitch = temp_load.trial_para.cumulativeCorrect_greenRedSwitch(all_outlierTrial{file_index}==0);
    temp_load.trial_para.consecutiveCorrect = temp_load.trial_para.consecutiveCorrect(all_outlierTrial{file_index}==0);
    temp_load.trial_para.cumulativeError = temp_load.trial_para.cumulativeError(all_outlierTrial{file_index}==0);
    temp_load.trial_para.currentFrames = temp_load.trial_para.currentFrames(all_outlierTrial{file_index}==0);
    temp_load.trial_para.currentSequence = temp_load.trial_para.currentSequence(all_outlierTrial{file_index}==0);
    temp_load.trial_para.isFilled = temp_load.trial_para.isFilled(all_outlierTrial{file_index}==0);
    temp_load.trial_para.taskMode_noChoice_Green0_Red1 = temp_load.trial_para.taskMode_noChoice_Green0_Red1(all_outlierTrial{file_index}==0);
    temp_load.trial_para.selectFlag_newSeq = temp_load.trial_para.selectFlag_newSeq(all_outlierTrial{file_index}==0);
    temp_load.trial_para.choiceMode = temp_load.trial_para.choiceMode(all_outlierTrial{file_index}==0);
    temp_load.trial_para.choiceCondition_flag = temp_load.trial_para.choiceCondition_flag(all_outlierTrial{file_index}==0);
    
    if isfield(temp_load.trial_para,'reCenterGazeBias_x')
        temp_load.trial_para.reCenterGazeBias_x = temp_load.trial_para.reCenterGazeBias_x(all_outlierTrial{file_index}==0);
        temp_load.trial_para.reCenterGazeBias_y = temp_load.trial_para.reCenterGazeBias_y(all_outlierTrial{file_index}==0);
    end
    
    %save single session in a mat file
    cd 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\preprocessed'
    file_name = ['preprocessed_' cell2mat(MAT_file_name(file_index))];    
    
    basic_para = temp_load.basic_para;
    trial_para = temp_load.trial_para;
    save(file_name, 'basic_para');
    save(file_name, 'trial_para','-append'); 

end



cd 'C:\ASDROOT\STUDY\MonkeyDataAnalysis'