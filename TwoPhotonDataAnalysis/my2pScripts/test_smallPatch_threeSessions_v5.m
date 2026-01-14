%% Initialization
% clear
close all

if_load = 1;

if_save = 1;


% temp_fileNameAB1 = '20230510A_and_20230511A';
% temp_fileNameAB2 = '20230510B_and_20230511A';
% temp_fileNameAB12_save = '20230510A_and_20230510B_and_20230511A';

temp_fileNameAB1 = '20230504A_and_20230509A';
temp_fileNameAB2 = '20230508A_and_20230509A';
temp_fileNameAB12_save = '20230504A_and_20230508A_and_20230509A';


output_shortPath = 'D:\twoPhotonData_motionCorrected';

output_pathAB1 = [output_shortPath '\twoSessions\' temp_fileNameAB1];
output_pathAB2 = [output_shortPath '\twoSessions\' temp_fileNameAB2];

if if_load == 1
    temp_load = load(output_pathAB1,'decodingDataSimplified_AB');
    decodingDataSimplified_AB1 = temp_load.decodingDataSimplified_AB;
    temp_load = load(output_pathAB2,'decodingDataSimplified_AB');
    decodingDataSimplified_AB2 = temp_load.decodingDataSimplified_AB;
end

a = 1;

tempMappingCellIndex_suite2p_AB1 = decodingDataSimplified_AB1.extraForMerged.tempMappingCellIndex_suite2p;
tempMappingCellIndex_suite2p_AB2 = decodingDataSimplified_AB2.extraForMerged.tempMappingCellIndex_suite2p;

tempMappingCellIndex_suite2p_AB1_B = decodingDataSimplified_AB1.extraForMerged.tempMappingCellIndex_suite2p(:,2);
tempMappingCellIndex_suite2p_AB2_B = decodingDataSimplified_AB2.extraForMerged.tempMappingCellIndex_suite2p(:,2);

temptempBoolIndex = ismember(tempMappingCellIndex_suite2p_AB1_B,tempMappingCellIndex_suite2p_AB2_B);

tempMappingCellIndex_suite2p_AB12_B_sorted = sort(tempMappingCellIndex_suite2p_AB1_B(temptempBoolIndex));
roiNum_AB12 = size(tempMappingCellIndex_suite2p_AB12_B_sorted,1);

tempMappingCellIndex_suite2p_AB12 = zeros(roiNum_AB12,3);
for tempi=1:roiNum_AB12
    temp_suite2p_B = tempMappingCellIndex_suite2p_AB12_B_sorted(tempi);
    
    temptempIndex_1 = find(tempMappingCellIndex_suite2p_AB1_B==temp_suite2p_B);
    temptempIndex_2 = find(tempMappingCellIndex_suite2p_AB2_B==temp_suite2p_B);
    
    tempMappingCellIndex_suite2p_AB12(tempi,1) = tempMappingCellIndex_suite2p_AB1(temptempIndex_1,1);
    tempMappingCellIndex_suite2p_AB12(tempi,2) = tempMappingCellIndex_suite2p_AB2(temptempIndex_2,1);
    tempMappingCellIndex_suite2p_AB12(tempi,3) = temp_suite2p_B;
end
tempMappingCellIndex_suite2p_AB12;

tempBoolIndex_AB1_B = ismember(tempMappingCellIndex_suite2p_AB1(:,2),tempMappingCellIndex_suite2p_AB12_B_sorted);
tempBoolIndex_AB2_B = ismember(tempMappingCellIndex_suite2p_AB2(:,2),tempMappingCellIndex_suite2p_AB12_B_sorted);


trialNum_A1 = decodingDataSimplified_AB1.extraForMerged.trialNum_A;
trialNum_A2 = decodingDataSimplified_AB2.extraForMerged.trialNum_A;
trialNum_B = decodingDataSimplified_AB2.extraForMerged.trialNum_B;
trialNum_AB12 = trialNum_A1 + trialNum_A2 + trialNum_B;

roiNum_AB1 = decodingDataSimplified_AB1.extraForMerged.roiNum_AB;
roiNum_AB2 = decodingDataSimplified_AB2.extraForMerged.roiNum_AB;

roiNum_A1 = decodingDataSimplified_AB1.extraForMerged.roiNum_A;
roiNum_A2 = decodingDataSimplified_AB2.extraForMerged.roiNum_A;

%% Get neuron data from A1, A2, B

% F_dff_decisionPeriodA
temp1 = decodingDataSimplified_AB1.F_dff_decisionPeriodA;
temp2 = decodingDataSimplified_AB2.F_dff_decisionPeriodA;
temp3 = zeros(roiNum_AB12,trialNum_A1+trialNum_A2+trialNum_B,size(temp1,3));
for tempi=1:roiNum_AB12   
    temp_suite2p_B = tempMappingCellIndex_suite2p_AB12(tempi,3);    
    temptempIndex_1 = find(tempMappingCellIndex_suite2p_AB1_B==temp_suite2p_B);
    temptempIndex_2 = find(tempMappingCellIndex_suite2p_AB2_B==temp_suite2p_B);
    
    temptemp_dff_A1 = squeeze(temp1(temptempIndex_1,1:trialNum_A1,:));
    temptemp_dff_A2 = squeeze(temp2(temptempIndex_2,1:trialNum_A2,:));        
    temptemp_dff_B = squeeze(temp2(temptempIndex_2,(trialNum_A2+1):end,:));
    
    temptemp_dff = [temptemp_dff_A1;temptemp_dff_A2;temptemp_dff_B];        
    temp3(tempi,:,:) = temptemp_dff;
end
F_dff_decisionPeriodA = temp3;


% F_dff_decisionPeriod
temp1 = decodingDataSimplified_AB1.F_dff_decisionPeriod;
temp2 = decodingDataSimplified_AB2.F_dff_decisionPeriod;
temp3 = zeros(roiNum_AB12,trialNum_A1+trialNum_A2+trialNum_B,size(temp1,3));
for tempi=1:roiNum_AB12   
    temp_suite2p_B = tempMappingCellIndex_suite2p_AB12(tempi,3);    
    temptempIndex_1 = find(tempMappingCellIndex_suite2p_AB1_B==temp_suite2p_B);
    temptempIndex_2 = find(tempMappingCellIndex_suite2p_AB2_B==temp_suite2p_B);
    
    temptemp_dff_A1 = squeeze(temp1(temptempIndex_1,1:trialNum_A1,:));
    temptemp_dff_A2 = squeeze(temp2(temptempIndex_2,1:trialNum_A2,:));        
    temptemp_dff_B = squeeze(temp2(temptempIndex_2,(trialNum_A2+1):end,:));
    
    temptemp_dff = [temptemp_dff_A1;temptemp_dff_A2;temptemp_dff_B];        
    temp3(tempi,:,:) = temptemp_dff;
end
F_dff_decisionPeriod = temp3;

a = 1;

% F_dff_baselinePeriod
temp1 = decodingDataSimplified_AB1.F_dff_baselinePeriod;
temp2 = decodingDataSimplified_AB2.F_dff_baselinePeriod;
temp3 = zeros(roiNum_AB12,trialNum_A1+trialNum_A2+trialNum_B,size(temp1,3));
for tempi=1:roiNum_AB12   
    temp_suite2p_B = tempMappingCellIndex_suite2p_AB12(tempi,3);    
    temptempIndex_1 = find(tempMappingCellIndex_suite2p_AB1_B==temp_suite2p_B);
    temptempIndex_2 = find(tempMappingCellIndex_suite2p_AB2_B==temp_suite2p_B);
    
    temptemp_dff_A1 = squeeze(temp1(temptempIndex_1,1:trialNum_A1,:));
    temptemp_dff_A2 = squeeze(temp2(temptempIndex_2,1:trialNum_A2,:));        
    temptemp_dff_B = squeeze(temp2(temptempIndex_2,(trialNum_A2+1):end,:));
    
    temptemp_dff = [temptemp_dff_A1;temptemp_dff_A2;temptemp_dff_B];        
    temp3(tempi,:,:) = temptemp_dff;
end
F_dff_baselinePeriod = temp3;


% F_dff_length1_sample
temp1 = decodingDataSimplified_AB1.F_dff_length1_sample;
temp2 = decodingDataSimplified_AB2.F_dff_length1_sample;
temp3 = zeros(roiNum_AB12,trialNum_A1+trialNum_A2+trialNum_B,size(temp1,3));
for tempi=1:roiNum_AB12   
    temp_suite2p_B = tempMappingCellIndex_suite2p_AB12(tempi,3);    
    temptempIndex_1 = find(tempMappingCellIndex_suite2p_AB1_B==temp_suite2p_B);
    temptempIndex_2 = find(tempMappingCellIndex_suite2p_AB2_B==temp_suite2p_B);
    
    temptemp_dff_A1 = squeeze(temp1(temptempIndex_1,1:trialNum_A1,:));
    temptemp_dff_A2 = squeeze(temp2(temptempIndex_2,1:trialNum_A2,:));        
    temptemp_dff_B = squeeze(temp2(temptempIndex_2,(trialNum_A2+1):end,:));
    
    temptemp_dff = [temptemp_dff_A1;temptemp_dff_A2;temptemp_dff_B];        
    temp3(tempi,:,:) = temptemp_dff;
end
F_dff_length1_sample = temp3;


% F_dff_length2_sample
temp1 = decodingDataSimplified_AB1.F_dff_length2_sample;
temp2 = decodingDataSimplified_AB2.F_dff_length2_sample;
temp3 = zeros(roiNum_AB12,trialNum_A1+trialNum_A2+trialNum_B,size(temp1,3));
for tempi=1:roiNum_AB12   
    temp_suite2p_B = tempMappingCellIndex_suite2p_AB12(tempi,3);    
    temptempIndex_1 = find(tempMappingCellIndex_suite2p_AB1_B==temp_suite2p_B);
    temptempIndex_2 = find(tempMappingCellIndex_suite2p_AB2_B==temp_suite2p_B);
    
    temptemp_dff_A1 = squeeze(temp1(temptempIndex_1,1:trialNum_A1,:));
    temptemp_dff_A2 = squeeze(temp2(temptempIndex_2,1:trialNum_A2,:));        
    temptemp_dff_B = squeeze(temp2(temptempIndex_2,(trialNum_A2+1):end,:));
    
    temptemp_dff = [temptemp_dff_A1;temptemp_dff_A2;temptemp_dff_B];        
    temp3(tempi,:,:) = temptemp_dff;
end
F_dff_length2_sample = temp3;


% F_dff_length3_sample
temp1 = decodingDataSimplified_AB1.F_dff_length3_sample;
temp2 = decodingDataSimplified_AB2.F_dff_length3_sample;
temp3 = zeros(roiNum_AB12,trialNum_A1+trialNum_A2+trialNum_B,size(temp1,3));
for tempi=1:roiNum_AB12   
    temp_suite2p_B = tempMappingCellIndex_suite2p_AB12(tempi,3);    
    temptempIndex_1 = find(tempMappingCellIndex_suite2p_AB1_B==temp_suite2p_B);
    temptempIndex_2 = find(tempMappingCellIndex_suite2p_AB2_B==temp_suite2p_B);
    
    temptemp_dff_A1 = squeeze(temp1(temptempIndex_1,1:trialNum_A1,:));
    temptemp_dff_A2 = squeeze(temp2(temptempIndex_2,1:trialNum_A2,:));        
    temptemp_dff_B = squeeze(temp2(temptempIndex_2,(trialNum_A2+1):end,:));
    
    temptemp_dff = [temptemp_dff_A1;temptemp_dff_A2;temptemp_dff_B];        
    temp3(tempi,:,:) = temptemp_dff;
end
F_dff_length3_sample = temp3;


%% Get trial data from A and B
seqIndex = [decodingDataSimplified_AB1.seqIndex(1:trialNum_A1),decodingDataSimplified_AB2.seqIndex(1:trialNum_A2),decodingDataSimplified_AB2.seqIndex(trialNum_A2+1:end)];
trial_para_currentSequence = [decodingDataSimplified_AB1.trial_para_currentSequence(1:trialNum_A1),decodingDataSimplified_AB2.trial_para_currentSequence(1:trialNum_A2),decodingDataSimplified_AB2.trial_para_currentSequence(trialNum_A2+1:end)];
trial_para_ifSelectOffloading = [decodingDataSimplified_AB1.trial_para_ifSelectOffloading(1:trialNum_A1),decodingDataSimplified_AB2.trial_para_ifSelectOffloading(1:trialNum_A2),decodingDataSimplified_AB2.trial_para_ifSelectOffloading(trialNum_A2+1:end)];
trial_para_isFilled = [decodingDataSimplified_AB1.trial_para_isFilled(1:trialNum_A1),decodingDataSimplified_AB2.trial_para_isFilled(1:trialNum_A2),decodingDataSimplified_AB2.trial_para_isFilled(trialNum_A2+1:end)];
trialIndex_bool_memoryCorrect = [decodingDataSimplified_AB1.trialIndex_bool_memoryCorrect(1:trialNum_A1),decodingDataSimplified_AB2.trialIndex_bool_memoryCorrect(1:trialNum_A2),decodingDataSimplified_AB2.trialIndex_bool_memoryCorrect(trialNum_A2+1:end)];

trial_para_choiceCondition_flag = [decodingDataSimplified_AB1.trial_para_choiceCondition_flag(1:trialNum_A1),decodingDataSimplified_AB2.trial_para_choiceCondition_flag(1:trialNum_A2),decodingDataSimplified_AB2.trial_para_choiceCondition_flag(trialNum_A2+1:end)];
trial_para_isCorrect = [decodingDataSimplified_AB1.trial_para_isCorrect(1:trialNum_A1),decodingDataSimplified_AB2.trial_para_isCorrect(1:trialNum_A2),decodingDataSimplified_AB2.trial_para_isCorrect(trialNum_A2+1:end)];


%% Generate decodingDataSimplified_AB12
decodingDataSimplified_AB12 = struct;
decodingDataSimplified_AB12.F_dff_baselinePeriod = F_dff_baselinePeriod;
decodingDataSimplified_AB12.baselinePeriod_interval = decodingDataSimplified_AB1.baselinePeriod_interval;
decodingDataSimplified_AB12.F_dff_length1_sample = F_dff_length1_sample;
decodingDataSimplified_AB12.length1_sample_interval = decodingDataSimplified_AB1.length1_sample_interval;
decodingDataSimplified_AB12.F_dff_length2_sample = F_dff_length2_sample;
decodingDataSimplified_AB12.length2_sample_interval = decodingDataSimplified_AB1.length2_sample_interval;
decodingDataSimplified_AB12.F_dff_length3_sample = F_dff_length3_sample;
decodingDataSimplified_AB12.length3_sample_interval = decodingDataSimplified_AB1.length3_sample_interval;

decodingDataSimplified_AB12.F_dff_decisionPeriodA = F_dff_decisionPeriodA;
decodingDataSimplified_AB12.decisionPeriodA_interval = decodingDataSimplified_AB1.decisionPeriodA_interval;
decodingDataSimplified_AB12.F_dff_decisionPeriod = F_dff_decisionPeriod;
decodingDataSimplified_AB12.decisionPeriod_interval = decodingDataSimplified_AB1.decisionPeriod_interval;
decodingDataSimplified_AB12.seqIndex = seqIndex;
decodingDataSimplified_AB12.trial_para_currentSequence = trial_para_currentSequence;
decodingDataSimplified_AB12.trialIndex_bool_memoryCorrect = trialIndex_bool_memoryCorrect;
decodingDataSimplified_AB12.target_seqSet_inOne = decodingDataSimplified_AB1.target_seqSet_inOne;
decodingDataSimplified_AB12.target_seqSet = decodingDataSimplified_AB1.target_seqSet;
decodingDataSimplified_AB12.trial_para_ifSelectOffloading = trial_para_ifSelectOffloading;
decodingDataSimplified_AB12.trial_para_isFilled = trial_para_isFilled;
decodingDataSimplified_AB12.trial_para_choiceCondition_flag = trial_para_choiceCondition_flag;
decodingDataSimplified_AB12.trial_para_isCorrect = trial_para_isCorrect;

% extra variable for merged sessions
extraForMerged = struct;
extraForMerged.roiNum_A1 = roiNum_A1;
extraForMerged.roiNum_A2 = roiNum_A2;
extraForMerged.roiNum_B = roiNum_B;
extraForMerged.roiNum_AB12 = roiNum_AB12;
extraForMerged.trialNum_A1 = trialNum_A1;
extraForMerged.trialNum_A2 = trialNum_A2;
extraForMerged.trialNum_B = trialNum_B;
extraForMerged.trialNum_AB12 = trialNum_AB12;
% extraForMerged.tempMappingCellIndex = tempMappingCellIndex_highCorr;
extraForMerged.tempMappingCellIndex_suite2p = tempMappingCellIndex_suite2p_AB12;
% extraForMerged.tempMappingCellIndex_suite2p_withCorr_sorted = tempMappingCellIndex_suite2p_withCorr_sorted;

decodingDataSimplified_AB12.extraForMerged = extraForMerged;

fprintf('roiNum_A1=%d,roiNum_A2=%d,roiNum_B=%d,roiNum_AB12=%d.\n',roiNum_A1,roiNum_A2,roiNum_B,roiNum_AB12);
fprintf('trialNum_A1=%d,trialNum_A2=%d,trialNum_B=%d,trialNum_AB12=%d.\n',trialNum_A1,trialNum_A2,trialNum_B,trialNum_AB12);


decodingDataSimplified_AB = decodingDataSimplified_AB12;

%% Save
if if_save == 1    
    output_path = [output_shortPath '\twoSessions\' temp_fileNameAB12_save];
    save(output_path,'decodingDataSimplified_AB');
end
%% End