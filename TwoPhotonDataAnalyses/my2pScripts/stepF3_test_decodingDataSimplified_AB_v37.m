% Chuan's 29th script (20251214)
% This script: Neuron alignment of two-session recording of the same FOV.
%% Initialization
% clear
% clear decodingDataSimplified_A decodingDataSimplified_B
close all

if_monkey_D0_Z1 = 1;

if_compute = 1;
if_plot = 1;
plotRoiNum = 100;

if_save = 0;



% % % 0122A-0123A (Z)
% currentSessionIndex_A = 7;
% currentSessionIndex_B = 8;

% % 0124A-0126A (Z)
% currentSessionIndex_A = 9;
% currentSessionIndex_B = 10;

% % 0129A-0131A (Z)
% currentSessionIndex_A = 11;
% currentSessionIndex_B = 12;

% % 0202A-0203A (Z)
% currentSessionIndex_A = 13;
% currentSessionIndex_B = 14;

% % 0207A-0208A (Z)
% currentSessionIndex_A = 15;
% currentSessionIndex_B = 16;

% % 0210A-0211A (Z)
% currentSessionIndex_A = 17;
% currentSessionIndex_B = 18;

% % 0216A-0218A (Z)
% currentSessionIndex_A = 19;
% currentSessionIndex_B = 20;

% % 0220A-0221A (Z)
% currentSessionIndex_A = 21;
% currentSessionIndex_B = 22;

% % 0226A-0227A (Z)
% currentSessionIndex_A = 23;
% currentSessionIndex_B = 24;

% % 0229A-0301A (Z)
% currentSessionIndex_A = 25;
% currentSessionIndex_B = 26;

% % 0304A-0305A (Z)
% currentSessionIndex_A = 27;
% currentSessionIndex_B = 28;

% % 0307A-0308A (Z)
% currentSessionIndex_A = 29;
% currentSessionIndex_B = 30;

% % 0312A-0315A (Z)
currentSessionIndex_A = 31;
currentSessionIndex_B = 32;

% % 0319A-0320A (Z)
% currentSessionIndex_A = 33;
% currentSessionIndex_B = 34;

% % 0322A-0323A (Z)
% currentSessionIndex_A = 35;
% currentSessionIndex_B = 36;

% % 0329A-0330A (Z)
% currentSessionIndex_A = 37;
% currentSessionIndex_B = 38;

% % 0402A-0403A (Z)
% currentSessionIndex_A = 39;
% currentSessionIndex_B = 40;

% % 0410A-0411A (Z)
% currentSessionIndex_A = 41;
% currentSessionIndex_B = 42;


% % 0527A-0528A
% currentSessionIndex_A = 19;%19
% currentSessionIndex_B = 20;%20

% % 0512A-0513A
% currentSessionIndex_A = 11;
% currentSessionIndex_B = 12;

% % 0530A-0531A
% currentSessionIndex_A = 21;
% currentSessionIndex_B = 22;


% % 0426A-0427A
% currentSessionIndex_A = 1;
% currentSessionIndex_B = 2;

% % 0502A-0503A
% currentSessionIndex_A = 3;
% currentSessionIndex_B = 4;

% % 0504A-0509A
% currentSessionIndex_A = 5;
% currentSessionIndex_B = 7;

% % 0508A-0509A
% currentSessionIndex_A = 6;
% currentSessionIndex_B = 7;

% % 0515A-0516A
% currentSessionIndex_A = 13;
% currentSessionIndex_B = 14;

% % 0525A-0526A
% currentSessionIndex_A = 17;
% currentSessionIndex_B = 18;

% % 0601A-0602A
% currentSessionIndex_A = 23;
% currentSessionIndex_B = 24;

% % 0604A-0605A
% currentSessionIndex_A = 25;
% currentSessionIndex_B = 26;

% % 0614A-0615A
% currentSessionIndex_A = 28;
% currentSessionIndex_B = 29;

% % 0510A-0511A
% currentSessionIndex_A = 8;
% currentSessionIndex_B = 10;

% % 0510B-0511A
% currentSessionIndex_A = 9;
% currentSessionIndex_B = 10;




r_corr_AB_threshold = -1;







currentSession_multi = string;

if if_monkey_D0_Z1 == 0
    currentSession_multi = [currentSession_multi; '113Recording_20230426A_Ding_Site16'];%1
    currentSession_multi = [currentSession_multi; '113Recording_20230427A_Ding_Site16_sameFOV0426'];%2
    currentSession_multi = [currentSession_multi; '113Recording_20230502A_Ding_Site13'];%3
    currentSession_multi = [currentSession_multi; '113Recording_20230503A_Ding_Site13_sameFOV0502'];%4
    currentSession_multi = [currentSession_multi; '113Recording_20230504A_Ding_Site02'];%5
    currentSession_multi = [currentSession_multi; '113Recording_20230508A_Ding_Site02_sameFOV0509'];%6
    currentSession_multi = [currentSession_multi; '113Recording_20230509A_Ding_Site02'];%7, 660000 frames, easy to crash
    
    currentSession_multi = [currentSession_multi; '113Recording_20230510A_Ding_Site05_sameFOV0511'];%8
    currentSession_multi = [currentSession_multi; '113Recording_20230510B_Ding_Site05_sameFOV0511'];%9
    currentSession_multi = [currentSession_multi; '113Recording_20230511A_Ding_Site05'];%10
    currentSession_multi = [currentSession_multi; '113Recording_20230512A_Ding_Site09'];%11
    currentSession_multi = [currentSession_multi; '113Recording_20230513A_Ding_Site09_sameFOV0512'];%12
    
    currentSession_multi = [currentSession_multi; '113Recording_20230515A_Ding_Site24_sameFOV0516'];%13
    currentSession_multi = [currentSession_multi; '113Recording_20230516A_Ding_Site24'];%14
    currentSession_multi = [currentSession_multi; '113Recording_20230517A_Ding_Site16B'];%15
    currentSession_multi = [currentSession_multi; '113Recording_20230522A_Ding_Site05B'];%16
    currentSession_multi = [currentSession_multi; '113Recording_20230525A_Ding_Site09B'];%17
    currentSession_multi = [currentSession_multi; '113Recording_20230526A_Ding_Site09B_sameFOV0525'];%18
    
    currentSession_multi = [currentSession_multi; '113Recording_20230527A_Ding_Site02B'];%19
    currentSession_multi = [currentSession_multi; '113Recording_20230528A_Ding_Site02B_sameFOV0527'];%20
    currentSession_multi = [currentSession_multi; '113Recording_20230530A_Ding_Site05C'];%21
    currentSession_multi = [currentSession_multi; '113Recording_20230531A_Ding_Site05C_sameFOV0530'];%22
    currentSession_multi = [currentSession_multi; '113Recording_20230601A_Ding_Site13B'];%23
    currentSession_multi = [currentSession_multi; '113Recording_20230602A_Ding_Site13B_sameFOV0601'];%24
    
    currentSession_multi = [currentSession_multi; '113Recording_20230604A_Ding_Site07'];%25
    currentSession_multi = [currentSession_multi; '113Recording_20230605A_Ding_Site07_sameFOV0604'];%26
    currentSession_multi = [currentSession_multi; '113Recording_20230612A_Ding_Site14'];%27
    currentSession_multi = [currentSession_multi; '113Recording_20230614A_Ding_Site05D'];%28
    currentSession_multi = [currentSession_multi; '113Recording_20230615A_Ding_Site05D_sameFOV0614'];%29
    currentSession_multi = [currentSession_multi; '113Recording_20230619A_Ding_Site02C'];%30
    currentSession_multi = [currentSession_multi; '113Recording_20230620A_Ding_Site05E'];%31
    
elseif if_monkey_D0_Z1 == 1
    currentSession_multi = [currentSession_multi; '113Recording_20240111A_Zelku_Site09A'];%1
    currentSession_multi = [currentSession_multi; '113Recording_20240112A_Zelku_Site06A'];%2
    currentSession_multi = [currentSession_multi; '113Recording_20240115A_Zelku_Site06A'];%3
    currentSession_multi = [currentSession_multi; '113Recording_20240117A_Zelku_Site14A'];%4
    currentSession_multi = [currentSession_multi; '113Recording_20240118A_Zelku_Site18A'];%5
    currentSession_multi = [currentSession_multi; '113Recording_20240119A_Zelku_Site17A'];%6
    
    currentSession_multi = [currentSession_multi; '113Recording_20240122A_Zelku_Site09B'];%7
    currentSession_multi = [currentSession_multi; '113Recording_20240123A_Zelku_Site09B_sameFOV0122'];%8
    
    currentSession_multi = [currentSession_multi; '113Recording_20240124A_Zelku_Site06B'];%9
    currentSession_multi = [currentSession_multi; '113Recording_20240126A_Zelku_Site06B_sameFOV0124'];%10
    
    currentSession_multi = [currentSession_multi; '113Recording_20240129A_Zelku_Site07A'];%11
    currentSession_multi = [currentSession_multi; '113Recording_20240131A_Zelku_Site07A_sameFOV0129'];%12
    
    currentSession_multi = [currentSession_multi; '113Recording_20240202A_Zelku_Site06XA'];%13
    currentSession_multi = [currentSession_multi; '113Recording_20240203A_Zelku_Site06XA_sameFOV0202'];%14
    
    currentSession_multi = [currentSession_multi; '113Recording_20240207A_Zelku_Site05A'];%15
    currentSession_multi = [currentSession_multi; '113Recording_20240208A_Zelku_Site05A_sameFOV0207'];%16

    currentSession_multi = [currentSession_multi; '113Recording_20240210A_Zelku_Site10A'];%17
    currentSession_multi = [currentSession_multi; '113Recording_20240211A_Zelku_Site10A_sameFOV0210'];%18
    
    currentSession_multi = [currentSession_multi; '113Recording_20240216A_Zelku_Site09C'];%19
    currentSession_multi = [currentSession_multi; '113Recording_20240218A_Zelku_Site09C_sameFOV0216'];%20
    
    currentSession_multi = [currentSession_multi; '113Recording_20240220A_Zelku_Site06XB'];%21
    currentSession_multi = [currentSession_multi; '113Recording_20240221A_Zelku_Site06XB_sameFOV0220'];%22
    
    currentSession_multi = [currentSession_multi; '113Recording_20240226A_Zelku_Site10B'];%23
    currentSession_multi = [currentSession_multi; '113Recording_20240227A_Zelku_Site10B_sameFOV0226'];%24
    
    currentSession_multi = [currentSession_multi; '113Recording_20240229A_Zelku_Site06C'];%25
    currentSession_multi = [currentSession_multi; '113Recording_20240301A_Zelku_Site06C_sameFOV0229'];%26
    
    currentSession_multi = [currentSession_multi; '113Recording_20240304A_Zelku_Site09D'];%27
    currentSession_multi = [currentSession_multi; '113Recording_20240305A_Zelku_Site09D_sameFOV0304'];%28
    
    currentSession_multi = [currentSession_multi; '113Recording_20240307A_Zelku_Site10C'];%29
    currentSession_multi = [currentSession_multi; '113Recording_20240308A_Zelku_Site10C_sameFOV0307'];%30
    
    currentSession_multi = [currentSession_multi; '113Recording_20240312A_Zelku_Site06RA'];%31
    currentSession_multi = [currentSession_multi; '113Recording_20240315A_Zelku_Site06RA_sameFOV0312'];%32
    
    currentSession_multi = [currentSession_multi; '113Recording_20240319A_Zelku_Site09E'];%33
    currentSession_multi = [currentSession_multi; '113Recording_20240320A_Zelku_Site09E_sameFOV0319'];%34
    
    currentSession_multi = [currentSession_multi; '113Recording_20240322A_Zelku_Site07B'];%35
    currentSession_multi = [currentSession_multi; '113Recording_20240323A_Zelku_Site07B_sameFOV0322'];%36
    
    currentSession_multi = [currentSession_multi; '113Recording_20240329A_Zelku_Site05B'];%37
    currentSession_multi = [currentSession_multi; '113Recording_20240330A_Zelku_Site05B_sameFOV0329'];%38
    
    currentSession_multi = [currentSession_multi; '113Recording_20240402A_Zelku_Site14B'];%39
    currentSession_multi = [currentSession_multi; '113Recording_20240403A_Zelku_Site14B_sameFOV0402'];%40
    
    currentSession_multi = [currentSession_multi; '113Recording_20240410A_Zelku_Site17B'];%41
    currentSession_multi = [currentSession_multi; '113Recording_20240411A_Zelku_Site17B_sameFOV0410'];%42
    
end

currentSession_multi(1) = [];
num_FOV = length(currentSession_multi);


% currentSession_A = '113Recording_20230527A_Ding_Site02B';
% currentSession_B = '113Recording_20230528A_Ding_Site02B_sameFOV0527';

currentSession_A = currentSession_multi{currentSessionIndex_A};
currentSession_B = currentSession_multi{currentSessionIndex_B};

fprintf('currentSession_A = %s.\n',currentSession_A);
fprintf('currentSession_B = %s.\n',currentSession_B);


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)

output_shortPath = 'D:\twoPhotonData_motionCorrected';
temp_if_max0_min1 = 0;

%% Load
currentSession = currentSession_A;
temp_currentSession_path = [output_shortPath '\' currentSession];
input_path_A = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
if exist('decodingDataSimplified_A','var') == 0
    temp_load = load([input_path_A,'\','decodingDataSimplified.mat'],'decodingDataSimplified');
    decodingDataSimplified_A = temp_load.decodingDataSimplified;
end

currentSession = currentSession_B;
temp_currentSession_path = [output_shortPath '\' currentSession];
input_path_B = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
if exist('decodingDataSimplified_B','var') == 0
    temp_load = load([input_path_B,'\','decodingDataSimplified.mat'],'decodingDataSimplified');
    decodingDataSimplified_B = temp_load.decodingDataSimplified;
end



%% Correct trial_para_ifSelectOffloading
decodingDataSimplified_A.trial_para_ifSelectOffloading = ...
    decodingDataSimplified_A.trial_para_ifSelectOffloading(1:length(decodingDataSimplified_A.trial_para_currentSequence));
decodingDataSimplified_B.trial_para_ifSelectOffloading = ...
    decodingDataSimplified_B.trial_para_ifSelectOffloading(1:length(decodingDataSimplified_B.trial_para_currentSequence));



%% Get trial data from A and B
seqIndex = [decodingDataSimplified_A.seqIndex,decodingDataSimplified_B.seqIndex];
trial_para_currentSequence = [decodingDataSimplified_A.trial_para_currentSequence,decodingDataSimplified_B.trial_para_currentSequence];
trial_para_ifSelectOffloading = [decodingDataSimplified_A.trial_para_ifSelectOffloading,decodingDataSimplified_B.trial_para_ifSelectOffloading];
trial_para_isFilled = [decodingDataSimplified_A.trial_para_isFilled,decodingDataSimplified_B.trial_para_isFilled];
trialIndex_bool_memoryCorrect = [decodingDataSimplified_A.trialIndex_bool_memoryCorrect,decodingDataSimplified_B.trialIndex_bool_memoryCorrect];
trial_para_choiceCondition_flag = [decodingDataSimplified_A.trial_para_choiceCondition_flag,decodingDataSimplified_B.trial_para_choiceCondition_flag];
trial_para_isCorrect = [decodingDataSimplified_A.trial_para_isCorrect,decodingDataSimplified_B.trial_para_isCorrect];

trialNum_A = length(decodingDataSimplified_A.trial_para_currentSequence);
trialNum_B = length(decodingDataSimplified_B.trial_para_currentSequence);
trialNum_AB = trialNum_A + trialNum_B;

a = 1; %#ok<*NASGU>

% %% Get neuron data from A and B
roiNum_A = size(decodingDataSimplified_A.F_dff_decisionPeriodA,1);
roiNum_B = size(decodingDataSimplified_B.F_dff_decisionPeriodA,1);

tempSessionIndex = currentSessionIndex_A;
tempCellIndex_currentFOV = find(selevtivity_multiFOV.cell_FOVIndex_multiFOV == tempSessionIndex);
dff_seqMean_delay1_currentFOV_A = selevtivity_multiFOV.dff_seqMean_delay1_multiFOV(tempCellIndex_currentFOV,:);
dff_seqMean_delay2_currentFOV_A = selevtivity_multiFOV.dff_seqMean_delay2_multiFOV(tempCellIndex_currentFOV,:);
dff_seqFrame_delay1_currentFOV_A = selevtivity_multiFOV.dff_seqFrame_delay1_multiFOV(tempCellIndex_currentFOV,:,:);
dff_seqFrame_delay2_currentFOV_A = selevtivity_multiFOV.dff_seqFrame_delay2_multiFOV(tempCellIndex_currentFOV,:,:);

tempSessionIndex = currentSessionIndex_B;
tempCellIndex_currentFOV = find(selevtivity_multiFOV.cell_FOVIndex_multiFOV == tempSessionIndex);
dff_seqMean_delay1_currentFOV_B = selevtivity_multiFOV.dff_seqMean_delay1_multiFOV(tempCellIndex_currentFOV,:);
dff_seqMean_delay2_currentFOV_B = selevtivity_multiFOV.dff_seqMean_delay2_multiFOV(tempCellIndex_currentFOV,:);
dff_seqFrame_delay1_currentFOV_B = selevtivity_multiFOV.dff_seqFrame_delay1_multiFOV(tempCellIndex_currentFOV,:,:);
dff_seqFrame_delay2_currentFOV_B = selevtivity_multiFOV.dff_seqFrame_delay2_multiFOV(tempCellIndex_currentFOV,:,:);


%% Get covMatrix_AB
if if_compute == 1
    covMatrix_AB_delay1 = zeros(roiNum_A,roiNum_B);
    covMatrix_AB_delay2 = zeros(roiNum_A,roiNum_B);
    p_covMatrix_AB_delay1 = zeros(roiNum_A,roiNum_B);
    p_covMatrix_AB_delay2 = zeros(roiNum_A,roiNum_B);
    % for tempi=1:roiNum_A
    parfor tempi=1:roiNum_A
        for tempj=1:roiNum_B
            [covMatrix_AB_delay1(tempi,tempj),p_covMatrix_AB_delay1(tempi,tempj)] = ...
                corr(dff_seqMean_delay1_currentFOV_A(tempi,:)',dff_seqMean_delay1_currentFOV_B(tempj,:)'); %#ok<*PFBNS>
            [covMatrix_AB_delay2(tempi,tempj),p_covMatrix_AB_delay2(tempi,tempj)] = ...
                corr(dff_seqMean_delay2_currentFOV_A(tempi,:)',dff_seqMean_delay2_currentFOV_B(tempj,:)'); %#ok<*PFBNS>
        end
    end
    covMatrix_AB_delay1_raw = covMatrix_AB_delay1;
    covMatrix_AB_delay2_raw = covMatrix_AB_delay2;
    covMatrix_AB_delay1(p_covMatrix_AB_delay1>0.05) = 0;
    covMatrix_AB_delay2(p_covMatrix_AB_delay2>0.05) = 0;
    
    covMatrix_AB = zeros(roiNum_A,roiNum_B);
    for tempi=1:roiNum_A
        for tempj=1:roiNum_B
            %covMatrix_AB(tempi,tempj) = max(covMatrix_AB_delay1(tempi,tempj),covMatrix_AB_delay2(tempi,tempj));
            covMatrix_AB(tempi,tempj) = mean([covMatrix_AB_delay1(tempi,tempj),covMatrix_AB_delay2(tempi,tempj)]);
        end
    end
    covMatrix_AB_raw = covMatrix_AB;
    covMatrix_AB(covMatrix_AB<0) = 0;
    %covMatrix_AB(covMatrix_AB<0.25) = 0;
end
covMatrix_AB;
if_compute;

%% Get meanSSEMatrix_AB, sse is sum of squares error of seqMean
meanSSEMatrix_AB_delay1 = zeros(roiNum_A,roiNum_B);
meanSSEMatrix_AB_delay2 = zeros(roiNum_A,roiNum_B);
for tempi=1:roiNum_A
    for tempj=1:roiNum_B
        temp_error = dff_seqMean_delay1_currentFOV_A(tempi,:)-dff_seqMean_delay1_currentFOV_B(tempj,:);
        meanSSEMatrix_AB_delay1(tempi,tempj) = sum(temp_error.^2);
        temp_error = dff_seqMean_delay2_currentFOV_A(tempi,:)-dff_seqMean_delay2_currentFOV_B(tempj,:);
        meanSSEMatrix_AB_delay2(tempi,tempj) = sum(temp_error.^2);
    end
end
meanSSEMatrix_AB = meanSSEMatrix_AB_delay1 + meanSSEMatrix_AB_delay2;

% meanSSEMatrix_AB_delay1(1,[2,993,446])
% meanSSEMatrix_AB_delay2(1,[2,993,446])
%
% meanSSEMatrix_AB(1,[2,993,446])



std_seqFrame_delay1_currentFOV_A = std(dff_seqFrame_delay1_currentFOV_A,0,3);
std_seqFrame_delay2_currentFOV_A = std(dff_seqFrame_delay2_currentFOV_A,0,3);
std_seqFrame_delay1_currentFOV_B = std(dff_seqFrame_delay1_currentFOV_B,0,3);
std_seqFrame_delay2_currentFOV_B = std(dff_seqFrame_delay2_currentFOV_B,0,3);


stdSSEMatrix_AB_delay1 = zeros(roiNum_A,roiNum_B);
stdSSEMatrix_AB_delay2 = zeros(roiNum_A,roiNum_B);
for tempi=1:roiNum_A
    for tempj=1:roiNum_B
        temptemp_stdA = std_seqFrame_delay1_currentFOV_A(tempi,:);
        temptemp_stdB = std_seqFrame_delay1_currentFOV_B(tempj,:);
        temp_error = temptemp_stdA - temptemp_stdB;
        stdSSEMatrix_AB_delay1(tempi,tempj) = sum(temp_error.^2);
        
        temptemp_stdA = std_seqFrame_delay2_currentFOV_A(tempi,:);
        temptemp_stdB = std_seqFrame_delay2_currentFOV_B(tempj,:);
        temp_error = temptemp_stdA - temptemp_stdB;
        stdSSEMatrix_AB_delay2(tempi,tempj) = sum(temp_error.^2);
    end
end

stdSSEMatrix_AB = stdSSEMatrix_AB_delay1 + stdSSEMatrix_AB_delay2;


% stdSSEMatrix_AB_delay1(1,[2,993,446])
% stdSSEMatrix_AB_delay2(1,[2,993,446])
%
% stdSSEMatrix_AB(1,[2,993,446])
%
%
% meanSSEMatrix_AB_delay1(1,[2,993,446]) .* stdSSEMatrix_AB_delay1(1,[2,993,446])
%
% meanSSEMatrix_AB_delay2(1,[2,993,446]) .* stdSSEMatrix_AB_delay2(1,[2,993,446])
%
% meanSSEMatrix_AB(1,[2,993,446]) .* stdSSEMatrix_AB(1,[2,993,446])

%% Get tempMappingCellIndex_suite2p

tempCellIndex_currentFOV= find(selevtivity_multiFOV.cell_FOVIndex_multiFOV == currentSessionIndex_A);
tempExactCellIndex_suite2p_currentFOV_A = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(tempCellIndex_currentFOV);

tempCellIndex_currentFOV= find(selevtivity_multiFOV.cell_FOVIndex_multiFOV == currentSessionIndex_B);
tempExactCellIndex_suite2p_currentFOV_B = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(tempCellIndex_currentFOV);


%% Get mapping matrix
script_mappingTwoSessions_afterPipeline_Name_v = autoGetFunName('test_mappingTwoSessions_afterPipeline',targetPATH);
script_mappingTwoSessions_afterPipeline= str2func(script_mappingTwoSessions_afterPipeline_Name_v);
if if_compute == 1
    script_mappingTwoSessions_afterPipeline();
end
valid_pairIndex_AandB;
tempMappingCellIndex = valid_pairIndex_AandB;

if false
    tempMappingCellIndex = temp_extra.tempMappingCellIndex;
end

roiNum_AB_shared = size(tempMappingCellIndex,1);

%% Get neuron data from A and B
a = 1;
% F_dff_baselinePeriod
temp1 = decodingDataSimplified_A.F_dff_baselinePeriod(tempMappingCellIndex(:,1),:,:);
temp2 = decodingDataSimplified_B.F_dff_baselinePeriod(tempMappingCellIndex(:,2),:,:);
temp3 = zeros(roiNum_AB_shared,trialNum_A+trialNum_B,size(decodingDataSimplified_A.F_dff_baselinePeriod,3));
temp3(:,1:trialNum_A,:) = temp1;
temp3(:,trialNum_A+(1:trialNum_B),:) = temp2;
F_dff_baselinePeriod = temp3;

a = 1;

% F_dff_length1_sample

if false
    a1 = size(decodingDataSimplified_A.F_dff_length1_sample,3); %#ok<*UNRCH>
    b1 = size(decodingDataSimplified_B.F_dff_length1_sample,3);
    tempValidSize = min(a1,b1);
    decodingDataSimplified_A.F_dff_length1_sample = decodingDataSimplified_A.F_dff_length1_sample(:,:,1:tempValidSize);
    decodingDataSimplified_A.length1_sample_interval = [1,tempValidSize];
end

temp1 = decodingDataSimplified_A.F_dff_length1_sample(tempMappingCellIndex(:,1),:,:);
temp2 = decodingDataSimplified_B.F_dff_length1_sample(tempMappingCellIndex(:,2),:,:);
temp3 = zeros(roiNum_AB_shared,trialNum_A+trialNum_B,size(decodingDataSimplified_A.F_dff_length1_sample,3));
temp3(:,1:trialNum_A,:) = temp1;
temp3(:,trialNum_A+(1:trialNum_B),:) = temp2;
F_dff_length1_sample = temp3;

% F_dff_length2_sample
temp1 = decodingDataSimplified_A.F_dff_length2_sample(tempMappingCellIndex(:,1),:,:);
temp2 = decodingDataSimplified_B.F_dff_length2_sample(tempMappingCellIndex(:,2),:,:);
temp3 = zeros(roiNum_AB_shared,trialNum_A+trialNum_B,size(decodingDataSimplified_A.F_dff_length2_sample,3));
temp3(:,1:trialNum_A,:) = temp1;
temp3(:,trialNum_A+(1:trialNum_B),:) = temp2;
F_dff_length2_sample = temp3;

% F_dff_length3_sample
temp1 = decodingDataSimplified_A.F_dff_length3_sample(tempMappingCellIndex(:,1),:,:);
temp2 = decodingDataSimplified_B.F_dff_length3_sample(tempMappingCellIndex(:,2),:,:);
temp3 = zeros(roiNum_AB_shared,trialNum_A+trialNum_B,size(decodingDataSimplified_A.F_dff_length3_sample,3));
temp3(:,1:trialNum_A,:) = temp1;
temp3(:,trialNum_A+(1:trialNum_B),:) = temp2;
F_dff_length3_sample = temp3;

% F_dff_decisionPeriodB
temp1 = decodingDataSimplified_A.F_dff_decisionPeriodB(tempMappingCellIndex(:,1),:,:);
temp2 = decodingDataSimplified_B.F_dff_decisionPeriodB(tempMappingCellIndex(:,2),:,:);
temp3 = zeros(roiNum_AB_shared,trialNum_A+trialNum_B,size(decodingDataSimplified_A.F_dff_decisionPeriodB,3));
temp3(:,1:trialNum_A,:) = temp1;
temp3(:,trialNum_A+(1:trialNum_B),:) = temp2;
F_dff_decisionPeriodB = temp3;


% F_dff_decisionPeriodA
temp1 = decodingDataSimplified_A.F_dff_decisionPeriodA(tempMappingCellIndex(:,1),:,:);
temp2 = decodingDataSimplified_B.F_dff_decisionPeriodA(tempMappingCellIndex(:,2),:,:);
temp3 = zeros(roiNum_AB_shared,trialNum_A+trialNum_B,size(decodingDataSimplified_A.F_dff_decisionPeriodA,3));
temp3(:,1:trialNum_A,:) = temp1;
temp3(:,trialNum_A+(1:trialNum_B),:) = temp2;
F_dff_decisionPeriodA = temp3;

% F_dff_decisionPeriod
temp1 = decodingDataSimplified_A.F_dff_decisionPeriod(tempMappingCellIndex(:,1),:,:);
temp2 = decodingDataSimplified_B.F_dff_decisionPeriod(tempMappingCellIndex(:,2),:,:);
temp3 = zeros(roiNum_AB_shared,trialNum_A+trialNum_B,size(decodingDataSimplified_A.F_dff_decisionPeriod,3));
temp3(:,1:trialNum_A,:) = temp1;
temp3(:,trialNum_A+(1:trialNum_B),:) = temp2;
F_dff_decisionPeriod = temp3;


%% Evaluation from function, Similarity of dff (delay 1) in each seqMean between A and B
% tempSessionIndex = currentSessionIndex_A;
% tempCellIndex_currentFOV = find(selevtivity_multiFOV.cell_FOVIndex_multiFOV == tempSessionIndex);
% dff_seqMean_delay1_currentFOV_A = selevtivity_multiFOV.dff_seqMean_delay1_multiFOV(tempCellIndex_currentFOV,:);
% dff_seqMean_delay2_currentFOV_A = selevtivity_multiFOV.dff_seqMean_delay2_multiFOV(tempCellIndex_currentFOV,:);
%
% tempSessionIndex = currentSessionIndex_B;
% tempCellIndex_currentFOV = find(selevtivity_multiFOV.cell_FOVIndex_multiFOV == tempSessionIndex);
% dff_seqMean_delay1_currentFOV_B = selevtivity_multiFOV.dff_seqMean_delay1_multiFOV(tempCellIndex_currentFOV,:);
% dff_seqMean_delay2_currentFOV_B = selevtivity_multiFOV.dff_seqMean_delay2_multiFOV(tempCellIndex_currentFOV,:);

dff_seqMean_delay1_currentFOV_A_shared = dff_seqMean_delay1_currentFOV_A(tempMappingCellIndex(:,1),:);
dff_seqMean_delay1_currentFOV_B_shared = dff_seqMean_delay1_currentFOV_B(tempMappingCellIndex(:,2),:);
dff_seqMean_delay2_currentFOV_A_shared = dff_seqMean_delay2_currentFOV_A(tempMappingCellIndex(:,1),:);
dff_seqMean_delay2_currentFOV_B_shared = dff_seqMean_delay2_currentFOV_B(tempMappingCellIndex(:,2),:);

r_corr_AB_delay1 = zeros(roiNum_AB_shared,1);
p_corr_AB_delay1 = zeros(roiNum_AB_shared,1);
for tempi=1:roiNum_AB_shared
    %[r_corr_AB_delay1(tempi),p_corr_AB_delay1(tempi)] = ...
    %    corr(dff_seqMean_delay1_currentFOV_A_shared(tempi,:)',dff_seqMean_delay1_currentFOV_B_shared(tempi,:)');
    
    temp1 = dff_seqMean_delay1_currentFOV_A_shared(tempi,:)';
    temp2 = dff_seqMean_delay1_currentFOV_B_shared(tempi,:)';
    
    tempValidIndex = (~isnan(temp1)) & (~isnan(temp2));
    
    [r_corr_AB_delay1(tempi),p_corr_AB_delay1(tempi)] = ...
        corr(temp1(tempValidIndex),temp2(tempValidIndex));
end
% r_corr_AB_delay1(p_corr_AB_delay1>0.05) = 0;

% a1 = sum(r_corr_AB_delay1<0 & p_corr_AB_delay1<0.05);

r_corr_AB_delay2 = zeros(roiNum_AB_shared,1);
p_corr_AB_delay2 = zeros(roiNum_AB_shared,1);
for tempi=1:roiNum_AB_shared
    [r_corr_AB_delay2(tempi),p_corr_AB_delay2(tempi)] = ...
        corr(dff_seqMean_delay2_currentFOV_A_shared(tempi,:)',dff_seqMean_delay2_currentFOV_B_shared(tempi,:)');
end
% r_corr_AB_delay2(p_corr_AB_delay2>0.05) = 0;

a1 = sum(r_corr_AB_delay1>r_corr_AB_threshold);
a2 = sum(r_corr_AB_delay2>r_corr_AB_threshold);

highCorrBoolIndex_AB_shared_delay1 = r_corr_AB_delay1>r_corr_AB_threshold;
highCorrBoolIndex_AB_shared_delay2 = r_corr_AB_delay2>r_corr_AB_threshold;

% highCorrBoolIndex_AB_shared = r_corr_AB>r_corr_AB_threshold;
highCorrBoolIndex_AB_shared = highCorrBoolIndex_AB_shared_delay1 | highCorrBoolIndex_AB_shared_delay2;
% fprintf('num(highCorrBoolIndex_AB_shared) = %d.\n',sum(highCorrBoolIndex_AB_shared));
a = 1;

% r_corr_AB = max([r_corr_AB_delay1,r_corr_AB_delay2],[],2);
r_corr_AB = mean([r_corr_AB_delay1,r_corr_AB_delay2],2);


% %%
% tempMappingCellIndex;
%% Get tempMappingCellIndex_suite2p
%
% tempCellIndex_currentFOV= find(selevtivity_multiFOV.cell_FOVIndex_multiFOV == currentSessionIndex_A);
% tempExactCellIndex_suite2p_currentFOV_A = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(tempCellIndex_currentFOV);
%
% tempCellIndex_currentFOV= find(selevtivity_multiFOV.cell_FOVIndex_multiFOV == currentSessionIndex_B);
% tempExactCellIndex_suite2p_currentFOV_B = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(tempCellIndex_currentFOV);
%
tempExactCellIndex_suite2p_currentFOV_A;
tempExactCellIndex_suite2p_currentFOV_B;

tempMappingCellIndex_suite2p = zeros(roiNum_AB_shared,2);
tempMappingCellIndex_suite2p(:,1) = tempExactCellIndex_suite2p_currentFOV_A(tempMappingCellIndex(:,1));
tempMappingCellIndex_suite2p(:,2) = tempExactCellIndex_suite2p_currentFOV_B(tempMappingCellIndex(:,2));



%% Apply highCorrBoolIndex_AB_shared
F_dff_baselinePeriod_highCorr = F_dff_baselinePeriod(highCorrBoolIndex_AB_shared,:,:);
F_dff_length1_sample_highCorr = F_dff_length1_sample(highCorrBoolIndex_AB_shared,:,:);
F_dff_length2_sample_highCorr = F_dff_length2_sample(highCorrBoolIndex_AB_shared,:,:);
F_dff_length3_sample_highCorr = F_dff_length3_sample(highCorrBoolIndex_AB_shared,:,:);

F_dff_decisionPeriodB_highCorr = F_dff_decisionPeriodB(highCorrBoolIndex_AB_shared,:,:);
F_dff_decisionPeriodA_highCorr = F_dff_decisionPeriodA(highCorrBoolIndex_AB_shared,:,:);
F_dff_decisionPeriod_highCorr = F_dff_decisionPeriod(highCorrBoolIndex_AB_shared,:,:);

tempMappingCellIndex_highCorr = tempMappingCellIndex(highCorrBoolIndex_AB_shared,:);
tempMappingCellIndex_suite2p_highCorr = tempMappingCellIndex_suite2p(highCorrBoolIndex_AB_shared,:);

tempMappingCellIndex_lowCorr = tempMappingCellIndex(~highCorrBoolIndex_AB_shared,:);
tempMappingCellIndex_suite2p_lowCorr = tempMappingCellIndex_suite2p(~highCorrBoolIndex_AB_shared,:);


roiNum_AB_shared_highCorr = sum(highCorrBoolIndex_AB_shared);


%% test example roi correlation
% temptemp_cellIndex_suite2p_A = 41;%41
% temptemp_cellIndex_suite2p_B = 1148;%1148
% 
% temptemp_cellIndex_A = find(tempExactCellIndex_suite2p_currentFOV_A==temptemp_cellIndex_suite2p_A);
% temptemp_cellIndex_B = find(tempExactCellIndex_suite2p_currentFOV_B==temptemp_cellIndex_suite2p_B);
% 
% x = dff_seqMean_delay1_currentFOV_A(temptemp_cellIndex_A,:)';
% y = dff_seqMean_delay1_currentFOV_B(temptemp_cellIndex_B,:)';
% [temp_r,temp_p] = corr(x,y);

%% Sorting tempMappingCellIndex_suite2p by correlation
tempMappingCellIndex_suite2p_withCorr = [tempMappingCellIndex_suite2p,r_corr_AB];
tempMappingCellIndex_suite2p_highCorr_withCorr = [tempMappingCellIndex_suite2p_highCorr,r_corr_AB(highCorrBoolIndex_AB_shared)];
tempMappingCellIndex_suite2p_lowCorr_withCorr = [tempMappingCellIndex_suite2p_lowCorr,r_corr_AB(~highCorrBoolIndex_AB_shared)];

[M,I] = sort(tempMappingCellIndex_suite2p_withCorr(:,3),'descend'); %#ok<*ASGLU>

tempMappingCellIndex_suite2p_withCorr_sorted = tempMappingCellIndex_suite2p_withCorr(I,:);

[M,I] = sort(tempMappingCellIndex_suite2p_withCorr(:,2));
tempMappingCellIndex_suite2p_B_withCorr = tempMappingCellIndex_suite2p_withCorr(I,:);

%% Generate decodingDataSimplified_AB
decodingDataSimplified_AB = struct;
decodingDataSimplified_AB.F_dff_baselinePeriod = F_dff_baselinePeriod_highCorr;
decodingDataSimplified_AB.baselinePeriod_interval = decodingDataSimplified_A.baselinePeriod_interval;
decodingDataSimplified_AB.F_dff_length1_sample = F_dff_length1_sample_highCorr;
decodingDataSimplified_AB.length1_sample_interval = decodingDataSimplified_A.length1_sample_interval;
decodingDataSimplified_AB.F_dff_length2_sample = F_dff_length2_sample_highCorr;
decodingDataSimplified_AB.length2_sample_interval = decodingDataSimplified_A.length2_sample_interval;
decodingDataSimplified_AB.F_dff_length3_sample = F_dff_length3_sample_highCorr;
decodingDataSimplified_AB.length3_sample_interval = decodingDataSimplified_A.length3_sample_interval;

decodingDataSimplified_AB.F_dff_decisionPeriodB = F_dff_decisionPeriodB_highCorr;
decodingDataSimplified_AB.decisionPeriodB_interval = decodingDataSimplified_A.decisionPeriodB_interval;
decodingDataSimplified_AB.F_dff_decisionPeriodA = F_dff_decisionPeriodA_highCorr;
decodingDataSimplified_AB.decisionPeriodA_interval = decodingDataSimplified_A.decisionPeriodA_interval;
decodingDataSimplified_AB.F_dff_decisionPeriod = F_dff_decisionPeriod_highCorr;
decodingDataSimplified_AB.decisionPeriod_interval = decodingDataSimplified_A.decisionPeriod_interval;
decodingDataSimplified_AB.seqIndex = seqIndex;
decodingDataSimplified_AB.trial_para_currentSequence = trial_para_currentSequence;
decodingDataSimplified_AB.trialIndex_bool_memoryCorrect = trialIndex_bool_memoryCorrect;
decodingDataSimplified_AB.target_seqSet_inOne = decodingDataSimplified_A.target_seqSet_inOne;
decodingDataSimplified_AB.target_seqSet = decodingDataSimplified_A.target_seqSet;
decodingDataSimplified_AB.trial_para_ifSelectOffloading = trial_para_ifSelectOffloading;
decodingDataSimplified_AB.trial_para_isFilled = trial_para_isFilled;
decodingDataSimplified_AB.trial_para_choiceCondition_flag = trial_para_choiceCondition_flag;
decodingDataSimplified_AB.trial_para_isCorrect = trial_para_isCorrect;

% extra variable for merged sessions
extraForMerged = struct;
extraForMerged.roiNum_A = roiNum_A;
extraForMerged.roiNum_B = roiNum_B;
extraForMerged.roiNum_AB = roiNum_AB_shared_highCorr;
extraForMerged.trialNum_A = trialNum_A;
extraForMerged.trialNum_B = trialNum_B;
extraForMerged.trialNum_AB = trialNum_AB;
extraForMerged.tempMappingCellIndex = tempMappingCellIndex_highCorr;
extraForMerged.tempMappingCellIndex_suite2p = tempMappingCellIndex_suite2p_highCorr;
extraForMerged.tempMappingCellIndex_suite2p_withCorr_sorted = tempMappingCellIndex_suite2p_withCorr_sorted;

decodingDataSimplified_AB.extraForMerged = extraForMerged;

fprintf('roiNum_A=%d,roiNum_B=%d,roiNum_AB=%d.\n',roiNum_A,roiNum_B,roiNum_AB_shared_highCorr);

%% Save
if if_save == 1
    temp_fileName = [currentSession_A(14:22),'_and_',currentSession_B(14:22)];
    output_path = [output_shortPath '\twoSessions\' temp_fileName];
    save(output_path,'decodingDataSimplified_AB');
end


% if true
%     aa1=decodingDataSimplified_AB.extraForMerged.tempMappingCellIndex_suite2p;
%     aa2=tempMappingCellIndex_suite2p_highCorr;
%     
%     temptempBoolIndex = nan(1,size(aa1,1));
%     for tempi=1:size(aa1,1)
%         temp1 = find(aa1(tempi,1)==aa2(:,1));
%         
%         if aa1(tempi,2) == aa2(temp1,2)
%             temptempBoolIndex(tempi) = true;
%         else
%             temptempBoolIndex(tempi) = false;
%         end
%     end
%     temp2 = find(temptempBoolIndex==false);
%     sort(aa1(temp2,2)')
%     length(temp2)
% end

%% End