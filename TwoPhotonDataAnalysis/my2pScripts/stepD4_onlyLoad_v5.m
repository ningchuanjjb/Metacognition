%% Initialization
% clear
close all

% if_monkey_D0_Z1 = 1;

if_dff_singleSession1_twoSession2_allSession3 = 2;%2

if_dff_singleSession1_twoSession2 = if_dff_singleSession1_twoSession2_allSession3;

if_resample_seqBased0_locBased1 = 1;

if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2 = 0;%0

if_trainingSet_memoryCorrect0_allMemory1 = if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2;

if_memoryPrecision_accuracy0_sigma1 = 0;


if if_monkey_D0_Z1 == 0
    currentSessionIndex_AB = 8;%8-->5-->4-->3-->7
    currentSessionIndex_A = 19;%19-->11-->10-->5-->17
    currentSessionIndex_B = 20;
    
    %     currentSessionIndex_AB = 2;
    %     currentSessionIndex_A = 3;
    %     currentSessionIndex_B = 4;
    
%         currentSessionIndex_AB = 6;% for spatial organization
%         currentSessionIndex_A = 13;
%         currentSessionIndex_B = 14;
    
    %currentSessionIndex_AB = 7;%8-->5-->4-->3-->7
    %currentSessionIndex_A = 17;%19-->11-->10-->5-->17
    
%         currentSessionIndex_AB = 3;
%         currentSessionIndex_A = 5;
%         currentSessionIndex_B = 7;
    
    
%     currentSessionIndex = 20;% 0528A
%     currentSessionIndex = 16;% 0522A
    currentSessionIndex = 5;% 0504A
%     currentSessionIndex = 18;% 0526A
%     currentSessionIndex = 19;% 0527A
    %currentSessionIndex = 3;% 0502A
    %currentSessionIndex = 26;% 0605A
    
elseif if_monkey_D0_Z1 == 1
%     currentSessionIndex_AB = 9;
%     currentSessionIndex_A = 23;
%     currentSessionIndex_B = 24;

    currentSessionIndex_AB = 4;
    currentSessionIndex_A = 13;
    currentSessionIndex_B = 14;

%     currentSessionIndex_AB = 10;
%     currentSessionIndex_A = 25;
%     currentSessionIndex_B = 26;

%     currentSessionIndex_AB = 13;
%     currentSessionIndex_A = 31;
%     currentSessionIndex_B = 32;

    
%     currentSessionIndex_AB = 2;
%     currentSessionIndex_A = 9;
    
    %currentSessionIndex = 8; 
    currentSessionIndex = 13;
end


if exist('if_compute_summary','var') == 1
    currentSessionIndex_AB = currentSessionIndex_AB_summary;
    
%     if if_monkey_D0_Z1 == 0
%         multiFOV_matrix_summary = ...
%             [1,2;
%             3,4;
%             5,7;
%             8,10;
%             11,12;
%             13,14;
%             17,18;
%             19,20;
%             21,22;
%             23,24;
%             25,26;
%             28,29];        
%         
%     elseif if_monkey_D0_Z1 == 1
%         multiFOV_matrix_summary = ...
%             [7,8;
%             9,10;
%             11,12;
%             13,14;
%             15,16;
%             17,18;
%             19,20;
%             21,22;
%             23,24;
%             25,26;
%             27,28;
%             29,30;
%             31,32;
%             33,34;
%             35,36;
%             37,38];            
%     end    
    
    currentSessionIndex_A = multiFOV_matrix_summary(currentSessionIndex_AB,1);
    currentSessionIndex_B = multiFOV_matrix_summary(currentSessionIndex_AB,2);
    
end


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)
output_shortPath = 'D:\twoPhotonData_motionCorrected';

if if_dff_singleSession1_twoSession2_allSession3 == 1
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
    
    currentSession = currentSession_multi{currentSessionIndex};
    
    fprintf('currentSession = %s.\n',currentSession);
    
    temp_currentSession_path = [output_shortPath '\' currentSession];
    temp_if_max0_min1 = 0;
    output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
    
elseif if_dff_singleSession1_twoSession2_allSession3 == 2
    
    currentABSession_multi = string;
    
    if if_monkey_D0_Z1 == 0
        %currentABSession_multi = [currentABSession_multi; '20230510A_and_20230510B_and_20230511A'];%1
        %currentABSession_multi = [currentABSession_multi; '20230512A_and_20230513A'];%2
        %currentABSession_multi = [currentABSession_multi; '20230527A_and_20230528A'];%3
        %currentABSession_multi = [currentABSession_multi; '20230530A_and_20230531A'];%4
        
        currentABSession_multi = [currentABSession_multi; '20230426A_and_20230427A'];%1 few
        currentABSession_multi = [currentABSession_multi; '20230502A_and_20230503A'];%2
        currentABSession_multi = [currentABSession_multi; '20230504A_and_20230508A_and_20230509A'];%3
        
        currentABSession_multi = [currentABSession_multi; '20230510A_and_20230510B_and_20230511A'];%4
        currentABSession_multi = [currentABSession_multi; '20230512A_and_20230513A'];%5 few
        
        currentABSession_multi = [currentABSession_multi; '20230515A_and_20230516A'];%6
        currentABSession_multi = [currentABSession_multi; '20230525A_and_20230526A'];%7
        currentABSession_multi = [currentABSession_multi; '20230527A_and_20230528A'];%8
        currentABSession_multi = [currentABSession_multi; '20230530A_and_20230531A'];%9
        currentABSession_multi = [currentABSession_multi; '20230601A_and_20230602A'];%10
        
        currentABSession_multi = [currentABSession_multi; '20230604A_and_20230605A'];%11
        currentABSession_multi = [currentABSession_multi; '20230614A_and_20230615A'];%12
        
        
        %currentABSession_multi = [currentABSession_multi; '20230508A_and_20230509A'];%3 few
        
    elseif if_monkey_D0_Z1 == 1
        
        currentABSession_multi = [currentABSession_multi; '20240122A_and_20240123A'];%1
        currentABSession_multi = [currentABSession_multi; '20240124A_and_20240126A'];%2        
        currentABSession_multi = [currentABSession_multi; '20240129A_and_20240131A'];%3        
        currentABSession_multi = [currentABSession_multi; '20240202A_and_20240203A'];%4        
        currentABSession_multi = [currentABSession_multi; '20240207A_and_20240208A'];%5        
        currentABSession_multi = [currentABSession_multi; '20240210A_and_20240211A'];%6        
        currentABSession_multi = [currentABSession_multi; '20240216A_and_20240218A'];%7        
        currentABSession_multi = [currentABSession_multi; '20240220A_and_20240221A'];%8        
        currentABSession_multi = [currentABSession_multi; '20240226A_and_20240227A'];%9        
        currentABSession_multi = [currentABSession_multi; '20240229A_and_20240301A'];%10        
        currentABSession_multi = [currentABSession_multi; '20240304A_and_20240305A'];%11        
        currentABSession_multi = [currentABSession_multi; '20240307A_and_20240308A'];%12        
        currentABSession_multi = [currentABSession_multi; '20240312A_and_20240315A'];%13        
        currentABSession_multi = [currentABSession_multi; '20240319A_and_20240320A'];%14        
        currentABSession_multi = [currentABSession_multi; '20240322A_and_20240323A'];%15        
        currentABSession_multi = [currentABSession_multi; '20240329A_and_20240330A'];%16        
        currentABSession_multi = [currentABSession_multi; '20240402A_and_20240403A'];%17        
        currentABSession_multi = [currentABSession_multi; '20240410A_and_20240411A'];%18        
        
    end

    
    currentABSession_multi(1) = [];
    num_FOV_AB = length(currentABSession_multi);
    
    currentSession = currentABSession_multi{currentSessionIndex_AB};
    
    fprintf('currentSession = %s.\n',currentSession);
    
    output_path = 'D:\twoPhotonData_motionCorrected\twoSessions';
elseif if_dff_singleSession1_twoSession2_allSession3 == 3
    if if_resample_seqBased0_locBased1 == 0
        currentSession = 'allSessions_resampled_seqBased';
    elseif if_resample_seqBased0_locBased1 == 1
        currentSession = 'allSessions_resampled_locBased';
    end
    if if_monkey_D0_Z1 == 0
        currentSession = [currentSession '_monkeyD'];
    elseif if_monkey_D0_Z1 == 1
        currentSession = [currentSession '_monkeyZ'];        
    end
    fprintf('currentSession = %s.\n',currentSession);
    output_path = 'D:\twoPhotonData_motionCorrected\twoSessions';
end

%% Params setting
if_compute_train = 1;%1

if_plot = 1;

if_gAcc_behavior0_model1 = 0;
if_gAcc_behavior_singleSession0_allSession1 = 1;

if_behaviorAllSession_gAcc0_rProb1 = 0;

if_decoder_decodingAcc0_pProd1 = 0;

if_decodingAcc_threshold0_sort1 = 1;


if_gAcc_noChoice0_allMemory1 = 1;

if if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2 == 2
    if_gAcc_noChoice0_allMemory1 = 0;
end

if_seqLevel0_trialLevel1 = 1;
if_decoder_6location0_56seq1 = 0;%0

valid_length = 3;


KFold_num = 5;%10-->5

% For reproducibility, delay1: 1, 4, 5
% For reproducibility, delay2: 2, 6
% rng(4);%2, 16, 17
rng(5);%5

if_normal0_resample1_leaveOneSeqOut2 = 1;


resampleTrialCount = 7;%7-->5-->7-->10(AB session)
% resampleTrialCount_min = resampleTrialCount;
resampleIterCount = 32;%16-->32-->

% if if_normal0_resample1_leaveOneSeqOut2 == 0
%     resampleIterCount = 1;
% end

if if_dff_singleSession1_twoSession2_allSession3 < 3
    if if_normal0_resample1_leaveOneSeqOut2 == 0
        resampleIterCount = 1;
        fprintf('Normal mode (no cross validation).\n');
    elseif if_normal0_resample1_leaveOneSeqOut2 == 1
        if if_resample_seqBased0_locBased1 == 0
            fprintf('Resample mode (seqBased). resampleTrialCount=%d.\n',resampleTrialCount);
        elseif if_resample_seqBased0_locBased1 == 1
            fprintf('Resample mode (locBased).\n');
        end
    elseif if_normal0_resample1_leaveOneSeqOut2 == 2
        fprintf('Leave-one-seq-out mode.\n');
    end
elseif if_dff_singleSession1_twoSession2_allSession3 == 3
    if if_resample_seqBased0_locBased1 == 0
        fprintf('Resample mode (seqBased). resampleTrialCount=%d.\n',resampleTrialCount);
    elseif if_resample_seqBased0_locBased1 == 1
        fprintf('Resample mode (locBased).\n');
    end
end


t_decoder = templateSVM('Standardize',true,'KernelFunction','linear'); % To standardise input data
% t_decoder = templateLinear('Regularization','lasso', 'Lambda', 0.1);



%% Load data
t0 = tic;
numFrames = 6;
pointKindsNum = 4;
target_seqSet = get_target_seqSet_v1(numFrames,pointKindsNum);
numSeq = zeros(1,pointKindsNum);
for tempi=1:pointKindsNum
    numSeq(tempi) = length(target_seqSet{tempi});
end

if_locDecoder_meanPosterior0_meanAcc1 = 0;

% s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_v11_20230831A']);
% s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_v12_onlyDiagonal_20230901A']);
% s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_v13_loopN11N_20230904A']);
% s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_v13_loopN11N_onlyDiagonal_20230904B']);
% s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_v13_loopN11N_onlyDiagonal_allMemory_20230907A']);
% s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_lp0_v3_20230831A']);
% s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_lp0_v5_loopN11N_20230904A']);
% s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_lp0_v5_loopN11N_onlyDiagonal_20230904B']);
% s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_f_l_v6_20230831A']);

if if_monkey_D0_Z1 == 0
    s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_v13_loopN11N_onlyDiagonal_allMemory_20230907A']);
elseif if_monkey_D0_Z1 == 1
    s = load(['C:\ASDROOT\STUDY\temp\data\','FittingResults_CCM_noOrder_distri_kwp0_v14_loopN11N_onlyDiagonal_allMemory_Zelku_20240412A']);
end

FittingResults = s.FittingResults;
behaivorMatrix_model = FittingResults.intermediate.RetrievalPMF_total';
gAcc_noChoice_inOne_model = zeros(1,size(behaivorMatrix_model,1));
for tempi=1:size(behaivorMatrix_model,1)
    gAcc_noChoice_inOne_model(tempi) = behaivorMatrix_model(tempi,tempi);
end
gAcc_noChoice_model = cell(1,pointKindsNum);
for tempi=1:pointKindsNum
    temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    gAcc_noChoice_model{tempi} = gAcc_noChoice_inOne_model(temp_range)';
end

% trial_para = []; %#ok<*NASGU>

% if exist('decodingData','var') == 0 || if_compute_train == 1
%     load([output_path,'\','decodingData.mat'],'decodingData');
% end
if exist('decodingDataSimplified','var') == 0 || if_compute_train == 1
    %load([output_path,'\','decodingDataSimplified.mat'],'decodingDataSimplified');
    if if_dff_singleSession1_twoSession2_allSession3 == 1
        temp_load = load([output_path,'\','decodingDataSimplified.mat'],'decodingDataSimplified');
        decodingDataSimplified = temp_load.decodingDataSimplified;
    elseif if_dff_singleSession1_twoSession2_allSession3 == 2
        temp_load = load([output_path '\' currentSession]);
        decodingDataSimplified = temp_load.decodingDataSimplified_AB;
    elseif if_dff_singleSession1_twoSession2_allSession3 == 3
        temp_load = load([output_path,'\',currentSession]);
        decodingDataSimplified = temp_load.decodingDataSimplified_allSessions_resampled;
    end
end
% load([output_path,'\','trial_para.mat'],'trial_para');
% load([output_path,'\','markerParse_trialLevel.mat'],'markerParse_trialLevel');

% trial_para_isFilled = trial_para.isFilled;
% trial_para_ifSelectOffloading = trial_para.ifSelectOffloading;

% temp_fields = fieldnames(decodingData);
% for tempi=1:length(temp_fields)
%     eval([temp_fields{tempi},'=decodingData.(temp_fields{tempi});']);
% end
temp_fields = fieldnames(decodingDataSimplified);
for tempi=1:length(temp_fields)
    eval([temp_fields{tempi},'=decodingDataSimplified.(temp_fields{tempi});']);
end

if if_monkey_D0_Z1 == 0
    searchName_gAcc = 'from23-04-26to23-06-20_D_gAcc_1';
    searchName_rProb = 'from23-04-26to23-06-20_D_offloadingProb_1';    
    searchName_responseMatrix = 'D_responseMatrix_2p_B.mat';
elseif if_monkey_D0_Z1 == 1
    %searchName_gAcc = 'from01-08to01-23_Z_gAcc_1';
    %searchName_rProb = 'from01-08to01-23_Z_offloadingProb_1'; 
    %searchName_responseMatrix = 'Z_responseMatrix_early2p.mat';    
    searchName_gAcc = 'from01-08to03-30_Z_gAcc_1';
    searchName_rProb = 'from01-08to03-30_Z_offloadingProb_1'; 
    searchName_responseMatrix = 'Z_responseMatrix_almost2p.mat';        
end    
searchName_mappingParam = 'from23-04-26to23-06-20_D_mappingParam_1';

path_behavior = [output_shortPath '\behavior'];

% Load other processed results
s = load([path_behavior,'\',searchName_responseMatrix]);
responseMatrix_choice = s.responseMatrix_choice;
responseMatrix_noChoice = s.responseMatrix_noChoice;
responseMatrix_allMemory = responseMatrix_choice+responseMatrix_noChoice;

gAcc_allMemory_collapsed_inOne = zeros(1,size(responseMatrix_allMemory,1));
gAcc_allMemory_collapsed = cell(1,pointKindsNum);
for tempi=1:size(responseMatrix_allMemory,1)
    temp1 = responseMatrix_allMemory(tempi,tempi);
    gAcc_allMemory_collapsed_inOne(tempi) = temp1/sum(responseMatrix_allMemory(tempi,:));
end
for tempi=1:pointKindsNum
    temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
    gAcc_allMemory_collapsed{tempi} = gAcc_allMemory_collapsed_inOne(temp_range)';
end

load_gAcc = loadMat_single(searchName_gAcc, path_behavior);
gAcc_noChoice_collapsed_inOne = load_gAcc.gAcc_noChoice_collapsed_inOne;
gAcc_noChoice_collapsed = load_gAcc.gAcc_noChoice_collapsed;

if if_gAcc_noChoice0_allMemory1 == 0
    gAcc_target_collapsed_inOne = gAcc_noChoice_collapsed_inOne;
    gAcc_target_collapsed = gAcc_noChoice_collapsed;
elseif if_gAcc_noChoice0_allMemory1 == 1
    gAcc_target_collapsed_inOne = gAcc_allMemory_collapsed_inOne;
    gAcc_target_collapsed = gAcc_allMemory_collapsed;
end

% Load other processed results
load_mapping = loadMat_single(searchName_mappingParam, path_behavior);
seqSet_inOne = load_mapping.seqSet_inOne;


% Load other processed results
load_rProb = loadMat_single(searchName_rProb, path_behavior);
offloadingProb_collapsed = load_rProb.offloadingProb_all;

offloadingProb_inOne = [];
for tempi=1:pointKindsNum
    offloadingProb_inOne = [offloadingProb_inOne offloadingProb_collapsed{tempi}']; %#ok<*AGROW>
end

if if_behaviorAllSession_gAcc0_rProb1 == 1
    gAcc_target_collapsed_inOne = offloadingProb_inOne;
    gAcc_target_collapsed = offloadingProb_collapsed;
end

fprintf('if_gAcc_behavior0_model1=%d.\n',if_gAcc_behavior0_model1);
if if_gAcc_behavior0_model1 == 0
    fprintf('if_gAcc_behavior_singleSession0_allSession1=%d.\n',if_gAcc_behavior_singleSession0_allSession1);
end

temp_range = 1:sum(numSeq(1:valid_length));
[temp_r,temp_p]=corr(gAcc_noChoice_inOne_model(temp_range)',gAcc_target_collapsed_inOne(temp_range)');
% fprintf('r_length123_model_behavior  =%.4f, p_length123_model_behavior  =%.4f.\n',temp_r,temp_p);

cd(targetPATH)

%% Preparation
if if_dff_singleSession1_twoSession2_allSession3 < 3
    % F_dff_decisionBin1 = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(2):decisionPeriod_interval(3)),3);
    F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
    F_dff_decisionBin1 = double(F_dff_decisionBin1);
    F_dff_decisionBin1 = F_dff_decisionBin1 + eps;
    a = 1; %#ok<*NASGU>
    
    temp_F_dff_decisionBin = F_dff_decisionBin1;
    
    roiNum = size(temp_F_dff_decisionBin,1);
    
    if if_seqLevel0_trialLevel1 == 0
        % abandoned now
    elseif if_seqLevel0_trialLevel1 == 1
        
        %     temp_trialIndex_bool = trialIndex_lengthx_bool(2,:);
        %temp_trialIndex_bool = trialIndex_lengthx_bool_memoryCorrect(1,:);
        %     temp_trialIndex_bool = trialIndex_lengthx_bool_memoryCorrect(2,:);
        %temp_trialIndex_bool = trialIndex_lengthx_bool_memoryCorrect(1,:) | trialIndex_lengthx_bool_memoryCorrect(2,:) | trialIndex_lengthx_bool_memoryCorrect(3,:);
        %temp_trialIndex_bool = trialIndex_bool_memoryCorrect; % use memory correct trial
        %     temp_trialIndex_bool = true(1,size(F_dff_decisionBin1,2)); % use all trial
        %temp_trialIndex_bool = trialIndex_bool_choiceMemory; % use choice memory trial
        
        
        if if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2 == 0
            temp_trialIndex_bool = trialIndex_bool_memoryCorrect; % use memory correct trial
        elseif if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2 == 1
            temp_trialIndex_bool = trial_para_ifSelectOffloading==-1; % use all memory trial
        elseif if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2 == 2
            temp_trialIndex_bool = trialIndex_bool_memoryCorrect & (trial_para_choiceCondition_flag==0); % use memory correct trial
        end
        
        temp_F_dff_decisionBin = F_dff_decisionBin1(:,temp_trialIndex_bool);
        
        seqIndex_valid = seqIndex(temp_trialIndex_bool);
        
        boolIndex_location_trial = false(numFrames,length(seqIndex_valid));
        for tempi=1:length(seqIndex_valid)
            currentSequence = target_seqSet_inOne{seqIndex_valid(tempi)};
            boolIndex_location_trial(currentSequence,tempi) = true;
        end
        
        boolIndex_location_seq = false(numFrames,length(target_seqSet_inOne));
        for tempi=1:length(target_seqSet_inOne)
            currentSequence = target_seqSet_inOne{tempi};
            boolIndex_location_seq(currentSequence,tempi) = true;
        end
        
        boolIndex_location_allTrial = false(numFrames,length(seqIndex));
        for tempi=1:length(seqIndex)
            currentSequence = target_seqSet_inOne{seqIndex(tempi)};
            boolIndex_location_allTrial(currentSequence,tempi) = true;
        end
        
        
        boolIndex_location_allTrial_response = false(numFrames,length(seqIndex));
        for tempi=1:length(seqIndex)
            %boolIndex_location_allTrial_response(:,tempi) = ~trial_para.isFilled{tempi};
            boolIndex_location_allTrial_response(:,tempi) = ~trial_para_isFilled{tempi};
        end
        
        seqIndex;
        seqIndex_response = nan(size(seqIndex));
        for tempi=1:length(seqIndex)
            temp_seqBoolIndex = boolIndex_location_allTrial_response(:,tempi)';
            
            for tempj=1:size(boolIndex_location_seq,2)
                temp_seqBoolIndex2 = boolIndex_location_seq(:,tempj)';
                
                if sum(temp_seqBoolIndex==temp_seqBoolIndex2) == numFrames
                    seqIndex_response(tempi) = tempj;
                    break
                end
            end            
        end
        
        
        %%
        seqIndex_valid;
        seqCount_length = cell(3,1);
        boolIndex_location_seq_T = boolIndex_location_seq';
        for tempi=1:3
            temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
            seqCount_length{tempi} = zeros(1,length(temp_range));
            for tempj=1:length(temp_range)
                temptempIndex = temp_range(tempj);
                temptempSeqBoolIndex = boolIndex_location_seq_T(temptempIndex,:);
                
                seqCount_length{tempi}(tempj) = sum(seqIndex_valid==temp_range(tempj),'all');
            end
        end
        seqCount_length
        
        seqCount = [];
        for tempi=1:3
            seqCount = [seqCount seqCount_length{tempi}];
        end
    end
elseif if_dff_singleSession1_twoSession2_allSession3 == 3
    roiNum = decodingDataSimplified.roiNum_allSessions_sum;
    
    boolIndex_location_seq = false(numFrames,length(target_seqSet_inOne));
    for tempi=1:length(target_seqSet_inOne)
        currentSequence = target_seqSet_inOne{tempi};
        boolIndex_location_seq(currentSequence,tempi) = true;
    end
    
    seqCount_valid_allSessions_collapsed = decodingDataSimplified.seqCount_valid_allSessions_collapsed;
    seqCount_valid_allSessions_collapsed_min = min(seqCount_valid_allSessions_collapsed,[],2);
    
    seqCount_min_length = cell(3,1);
    boolIndex_location_seq_T = boolIndex_location_seq';
    for tempi=1:3
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
        %seqCount_min_length{tempi} = zeros(1,length(temp_range));
        
        seqCount_min_length{tempi} = seqCount_valid_allSessions_collapsed_min(temp_range)';
    end
    seqCount_min_length
    
    seqIndex_valid_resampled = seqIndex_valid_resample_allSessions_collapsed;
    seqIndex_valid = seqIndex_valid_resample_allSessions_collapsed;
    temp_F_dff_decisionBin_resampled = temp_F_dff_decisionBin_resample_allSessions_collaped;
end
a = 1;


%%

if_beta_n11n = 1;
if_memoryPrecisionTuning_decoder0_behavior1 = 1;

%% Load trial data
temp_trial_para_ifSelectOffloading = trial_para_ifSelectOffloading==1;

choiceBoolIndex = trial_para_choiceCondition_flag == 2;
allMemoryBoolIndex = ~temp_trial_para_ifSelectOffloading;
choiceMemoryBoolIndex = choiceBoolIndex & allMemoryBoolIndex;
choiceOffloadBoolIndex = choiceBoolIndex & temp_trial_para_ifSelectOffloading;

allMemoryCorrectBoolIndex = trialIndex_bool_memoryCorrect;
allMemoryErrorBoolIndex = allMemoryBoolIndex & ~trialIndex_bool_memoryCorrect;

choiceMemoryCorrectBoolIndex = choiceMemoryBoolIndex & trialIndex_bool_memoryCorrect;
choiceMemoryErrorBoolIndex = choiceMemoryBoolIndex & ~trialIndex_bool_memoryCorrect;

trial_num = length(choiceBoolIndex);

if if_dff_singleSession1_twoSession2 == 1
elseif if_dff_singleSession1_twoSession2 == 2
    FOVName_currentFOV = currentSession;
    FOVName_currentFOV2 = [FOVName_currentFOV(5:9),'+',FOVName_currentFOV(19:23)];
end

if if_dff_singleSession1_twoSession2 == 1
    
elseif if_dff_singleSession1_twoSession2 == 2
    cellIndex_suite2p = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,1);
    
    
    temp_range_B = FOVAllCellRange_multiFOV(currentSessionIndex_B,1):FOVAllCellRange_multiFOV(currentSessionIndex_B,2);
    cellIndex_suite2p_B_raw = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range_B);
    cellIndex_suite2p_B = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end);
    temp_range_AB = nan(1,length(cellIndex_suite2p_B));
    for tempi=1:length(cellIndex_suite2p_B)
        temp_range_AB(tempi) = find(cellIndex_suite2p_B_raw==cellIndex_suite2p_B(tempi));
    end
    temp_range_AB = temp_range_AB + temp_range_B(1) - 1;
    
    cellIndex_suite2p_B2 = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range_AB);
    
    r2_6loc = glm_r2_6locMean_delay1Bin_multiFOV(temp_range_AB);
    beta_6loc = glm_beta_6locMean_delay1Bin_multiFOV(temp_range_AB,:);
    std_beta_6loc = std(beta_6loc,1,2);
        
end

a = 1;


%% Tuning for seq-level, r2_gProb
F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
F_dff_decisionBin1_seqMerged = nan(size(F_dff_decisionBin1,1),sum(numSeq));
for tempi=1:sum(numSeq)
    temp_dff = F_dff_decisionBin1(:,seqIndex==tempi);
    F_dff_decisionBin1_seqMerged(:,tempi) = mean(temp_dff,2);
end


r2_gProb = nan(size(F_dff_decisionBin1,1),1);
beta_gProb = nan(size(F_dff_decisionBin1,1),1);
p_gProb = nan(size(F_dff_decisionBin1,1),1);
parfor tempi=1:size(F_dff_decisionBin1,1)
    x = F_dff_decisionBin1_seqMerged(tempi,:)';
    y = (1-offloadingProb_inOne)';
    temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
    r2_gProb(tempi) = temp_mdl.Rsquared.Adjusted;
    beta_gProb(tempi) = temp_mdl.Coefficients.Estimate(2);
    p_gProb(tempi) = temp_mdl.Coefficients.pValue(2);
end
beta_gProb_raw = beta_gProb;
if if_beta_n11n == 1
    beta_gProb = rescale(beta_gProb_raw,-1,1);
end

%% Tuning for seq-level, r2_seqPrecision
r2_seqPrecision = nan(size(F_dff_decisionBin1,1),1);
beta_seqPrecision = nan(size(F_dff_decisionBin1,1),1);
p_seqPrecision = nan(size(F_dff_decisionBin1,1),1);
parfor tempi=1:size(F_dff_decisionBin1,1)
    %y = seqPrecision_behavior;
    %y = seqPrecision_behavior_56;
    y = gAcc_target_collapsed_inOne';    
    x = F_dff_decisionBin1_seqMerged(tempi,1:length(y))'; %#ok<*PFBNS>
    temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
    r2_seqPrecision(tempi) = temp_mdl.Rsquared.Adjusted;
    beta_seqPrecision(tempi) = temp_mdl.Coefficients.Estimate(2);
    p_seqPrecision(tempi) = temp_mdl.Coefficients.pValue(2);
end
beta_seqPrecision_raw = beta_seqPrecision;
if if_beta_n11n == 1
    beta_seqPrecision = rescale(beta_seqPrecision_raw,-1,1);
end

[M,I] = max(r2_seqPrecision); %#ok<*ASGLU>


%% Tuning for seq-level, r2_gProb_baseline
temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
% temptemp_range = (baselinePeriod_interval(1)):baselinePeriod_interval(3);
F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,temptemp_range),3);
F_dff_baselineBin_seqMerged = nan(size(F_dff_baselineBin,1),sum(numSeq));
for tempi=1:sum(numSeq)
    temp_dff = F_dff_baselineBin(:,seqIndex==tempi);
    F_dff_baselineBin_seqMerged(:,tempi) = mean(temp_dff,2);
end


r2_gProb_baseline = nan(size(F_dff_baselineBin,1),1);
beta_gProb_baseline = nan(size(F_dff_decisionBin1,1),1);
p_gProb_baseline = nan(size(F_dff_baselineBin,1),1);
parfor tempi=1:size(F_dff_baselineBin,1)
    x = F_dff_baselineBin_seqMerged(tempi,:)';
    y = (1-offloadingProb_inOne)';
    temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
    r2_gProb_baseline(tempi) = temp_mdl.Rsquared.Adjusted;
    beta_gProb_baseline(tempi) = temp_mdl.Rsquared.Adjusted;
    p_gProb_baseline(tempi) = temp_mdl.Coefficients.pValue(2);
end
beta_gProb_baseline_raw = beta_gProb_baseline;
if if_beta_n11n == 1
    beta_gProb_baseline = rescale(beta_gProb_baseline_raw,-1,1);
end


%% Tuning for trial-level, r2_memoryPrecision
if if_memoryPrecisionTuning_decoder0_behavior1 == 0

elseif if_memoryPrecisionTuning_decoder0_behavior1 == 1
    temp_y = nan(trial_num,1);
    temp_y(allMemoryCorrectBoolIndex) = 1;
    temp_y(allMemoryErrorBoolIndex) = 0;
end

x = F_dff_decisionBin1';
%y = memoryPrecision_trialLevel_resampleMean;
y = temp_y;

r2 = zeros(roiNum,1);
p = zeros(roiNum,1);
beta = zeros(roiNum,1);
r = zeros(roiNum,1);
parfor tempi=1:roiNum
    % for tempi=1:roiNum
    warning('off');
    temp_mdl = fitglm(x(:,tempi),y);
    r2(tempi) = temp_mdl.Rsquared.Adjusted;
    p(tempi) = temp_mdl.Coefficients.pValue(2);
    beta(tempi) = temp_mdl.Coefficients.Estimate(2);
    r(tempi) = corr(x(~isnan(y),tempi),y(~isnan(y)));
    warning('on');
end
r2; %#ok<*VUNUS>
p;

[M,I] = sort(r2,'descend');

r2_memoryPrecision = r2;
p_memoryPrecision = p;
beta_memoryPrecision = beta;
r_memoryPrecision = r;

beta_memoryPrecision_raw = beta_memoryPrecision;
if if_beta_n11n == 1
    beta_memoryPrecision = rescale(beta_memoryPrecision_raw,-1,1);
end

%% Tuning for trial-level, r2_choiceMemory
x_raw = F_dff_decisionBin1';
x = x_raw(choiceBoolIndex,:);
y_raw = choiceMemoryBoolIndex';
y = y_raw(choiceBoolIndex);


r2 = zeros(roiNum,1);
p = zeros(roiNum,1);
beta = zeros(roiNum,1);
r = zeros(roiNum,1);
parfor tempi=1:roiNum
    % for tempi=1:roiNum
    warning('off');
    temp_mdl = fitglm(x(:,tempi),y);
    r2(tempi) = temp_mdl.Rsquared.Adjusted;
    p(tempi) = temp_mdl.Coefficients.pValue(2);
    beta(tempi) = temp_mdl.Coefficients.Estimate(2);
    r(tempi) = corr(x(:,tempi),y);
    warning('on');
end
r2; %#ok<*VUNUS>
p;

[M,I] = sort(r2,'descend');

r2_choiceMemory = r2;
p_choiceMemory = p;
beta_choiceMemory = beta;
r_choiceMemory = r;

beta_choiceMemory_raw = beta_choiceMemory;
if if_beta_n11n == 1
    beta_choiceMemory = rescale(beta_choiceMemory_raw,-1,1);
end


%% Tuning for trial-level, r2_choiceMemory_baseline
temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,temptemp_range),3);

x_raw = F_dff_baselineBin';
x = x_raw(choiceBoolIndex,:);
y_raw = choiceMemoryBoolIndex';
y = y_raw(choiceBoolIndex);


r2 = zeros(roiNum,1);
p = zeros(roiNum,1);
beta = zeros(roiNum,1);
r = zeros(roiNum,1);
parfor tempi=1:roiNum
    % for tempi=1:roiNum
    warning('off');
    temp_mdl = fitglm(x(:,tempi),y);
    r2(tempi) = temp_mdl.Rsquared.Adjusted;
    p(tempi) = temp_mdl.Coefficients.pValue(2);
    beta(tempi) = temp_mdl.Coefficients.Estimate(2);
    r(tempi) = corr(x(:,tempi),y);
    warning('on');
end
r2; %#ok<*VUNUS>
p;

[M,I] = sort(r2,'descend');

r2_choiceMemory_baseline = r2;
p_choiceMemory_baseline = p;
beta_choiceMemory_baseline = beta;
r_choiceMemory_baseline = r;

beta_choiceMemory_baseline_raw = beta_choiceMemory_baseline;
if if_beta_n11n == 1
    beta_choiceMemory_baseline = rescale(beta_choiceMemory_baseline_raw,-1,1);
end


%%
p_memoryPrecision;
p_choiceMemory;
p_choiceMemory_baseline;

p_seqPrecision;
p_gProb;
p_gProb_baseline;


%% End