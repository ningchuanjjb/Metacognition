%% Initialization
% clear
close all

% if_monkey_D0_Z1 = 1;

if_dff_singleSession1_twoSession2_allSession3 = 2;%2
if_resample_seqBased0_locBased1 = 1;

if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 = 0;%0

if_trainingSet_memoryCorrect0_allMemory1 = if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3;

if_memoryPrecision_accuracy0_sigma1 = 1;%0

if_usePureNeuronDecoding = 0;%0

if if_monkey_D0_Z1 == 0
%     currentSessionIndex_AB = 8;%S02B
    %currentSessionIndex_A = 19;
    %currentSessionIndex_B = 20;
        
%     currentSessionIndex_AB = 3;%S02
    %currentSessionIndex_A = 5;
    %currentSessionIndex_B = 7;
    
    currentSessionIndex_AB = 4;%S05, good for linear axis
    %currentSessionIndex_A = 8;
    %currentSessionIndex_B = 10;   
    
%     currentSessionIndex_AB = 7;%S09B
    %currentSessionIndex_A = 17;
    %currentSessionIndex_B = 18;       
    
    
    %currentSessionIndex_AB = 9;%S05C
    %currentSessionIndex_A = 21;
    %currentSessionIndex_B = 22;    
    
    %currentSessionIndex_AB = 11;%S07
    %currentSessionIndex_A = 25;
    %currentSessionIndex_B = 26;
    
%     currentSessionIndex_AB = 12;%S05D
    %currentSessionIndex_A = 28;
    %currentSessionIndex_B = 29;    
    
    %currentSessionIndex_AB = 2;
    
    
    %     currentSessionIndex_AB = 2;
    %     currentSessionIndex_A = 3;
    %     currentSessionIndex_B = 4;
    
%         currentSessionIndex_AB = 6;% for spatial organization
%         currentSessionIndex_A = 13;
%         currentSessionIndex_B = 14;
    
    %currentSessionIndex_AB = 7;%8-->5-->4-->3-->7
    %currentSessionIndex_A = 17;%19-->11-->10-->5-->17
    
    
    
    currentSessionIndex = 20;% 0528A
%     currentSessionIndex = 16;% 0522A
%     currentSessionIndex = 5;% 0504A
%     currentSessionIndex = 18;% 0526A
%     currentSessionIndex = 19;% 0527A
    %currentSessionIndex = 3;% 0502A
    %currentSessionIndex = 26;% 0605A
    
elseif if_monkey_D0_Z1 == 1
    
%     currentSessionIndex_AB = 2;%S06B
    %currentSessionIndex_AB = 4;%S06XA
    %currentSessionIndex_AB = 6;%S10A
%     currentSessionIndex_AB = 9;%S10B, good (Fig. 2 & 3)
    
%     currentSessionIndex_AB = 10;%S06C
    currentSessionIndex_AB = 13;%S06RA, good (Fig. 5)
    %currentSessionIndex_AB = 14;%S09E
    %currentSessionIndex_AB = 16;%S05B

    
    %     currentSessionIndex_AB = 9;
    %     currentSessionIndex_A = 23;
    %     currentSessionIndex_B = 24;
    
    %currentSessionIndex_AB = 4;
    %currentSessionIndex_A = 13;
    %currentSessionIndex_B = 14;
    
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



if if_monkey_D0_Z1 == 0
    temp_multiFOV_matrix_summary = ...
        [1,2;
        3,4;
        5,7;
        8,10;
        11,12;
        13,14;
        17,18;
        19,20;
        21,22;
        23,24;
        25,26;
        28,29];
    
elseif if_monkey_D0_Z1 == 1
    temp_multiFOV_matrix_summary = ...
        [7,8;
        9,10;
        11,12;
        13,14;
        15,16;
        17,18;
        19,20;
        21,22;
        23,24;
        25,26;
        27,28;
        29,30;
        31,32;
        33,34;
        35,36;
        37,38];
end

currentSessionIndex_A = temp_multiFOV_matrix_summary(currentSessionIndex_AB,1);
currentSessionIndex_B = temp_multiFOV_matrix_summary(currentSessionIndex_AB,2);


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
        
        currentABSession_multi = [currentABSession_multi; '20230426A_and_20230427A'];%1 few, S16
        currentABSession_multi = [currentABSession_multi; '20230502A_and_20230503A'];%2, S13
        currentABSession_multi = [currentABSession_multi; '20230504A_and_20230508A_and_20230509A'];%3, S02
        
        currentABSession_multi = [currentABSession_multi; '20230510A_and_20230510B_and_20230511A'];%4, S05
        currentABSession_multi = [currentABSession_multi; '20230512A_and_20230513A'];%5 few, S09
        
        currentABSession_multi = [currentABSession_multi; '20230515A_and_20230516A'];%6, S24
        currentABSession_multi = [currentABSession_multi; '20230525A_and_20230526A'];%7, S09B
        currentABSession_multi = [currentABSession_multi; '20230527A_and_20230528A'];%8, S02B
        currentABSession_multi = [currentABSession_multi; '20230530A_and_20230531A'];%9, S05C
        currentABSession_multi = [currentABSession_multi; '20230601A_and_20230602A'];%10, S13B
        
        currentABSession_multi = [currentABSession_multi; '20230604A_and_20230605A'];%11, S07
        currentABSession_multi = [currentABSession_multi; '20230614A_and_20230615A'];%12, S05D
        
        
        %currentABSession_multi = [currentABSession_multi; '20230508A_and_20230509A'];%3 few
        
    elseif if_monkey_D0_Z1 == 1
        
        currentABSession_multi = [currentABSession_multi; '20240122A_and_20240123A'];%1, S09B
        currentABSession_multi = [currentABSession_multi; '20240124A_and_20240126A'];%2, S06B        
        currentABSession_multi = [currentABSession_multi; '20240129A_and_20240131A'];%3, S07A        
        currentABSession_multi = [currentABSession_multi; '20240202A_and_20240203A'];%4, S06XA        
        currentABSession_multi = [currentABSession_multi; '20240207A_and_20240208A'];%5, S05A        
        currentABSession_multi = [currentABSession_multi; '20240210A_and_20240211A'];%6, S10A        
        currentABSession_multi = [currentABSession_multi; '20240216A_and_20240218A'];%7, S09C        
        currentABSession_multi = [currentABSession_multi; '20240220A_and_20240221A'];%8, S06XB        
        currentABSession_multi = [currentABSession_multi; '20240226A_and_20240227A'];%9, S10B        
        currentABSession_multi = [currentABSession_multi; '20240229A_and_20240301A'];%10, S06C        
        currentABSession_multi = [currentABSession_multi; '20240304A_and_20240305A'];%11, S09D        
        currentABSession_multi = [currentABSession_multi; '20240307A_and_20240308A'];%12, S10C        
        currentABSession_multi = [currentABSession_multi; '20240312A_and_20240315A'];%13, S06RA        
        currentABSession_multi = [currentABSession_multi; '20240319A_and_20240320A'];%14, S09E        
        currentABSession_multi = [currentABSession_multi; '20240322A_and_20240323A'];%15, S07B        
        currentABSession_multi = [currentABSession_multi; '20240329A_and_20240330A'];%16, S05B  
        
        currentABSession_multi = [currentABSession_multi; '20240402A_and_20240403A'];%17, bad      
        currentABSession_multi = [currentABSession_multi; '20240410A_and_20240411A'];%18, bad   
        
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
if_gAcc_behavior_singleSession0_allSession1 = 1;%1

if_behaviorAllSession_gAcc0_rProb1 = 0;

if_decoder_decodingAcc0_pProd1 = 0;

if_decodingAcc_threshold0_sort1 = 1;


if_gAcc_noChoice0_allMemory1_choiceMemory2 = 1;

if if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 2
    if_gAcc_noChoice0_allMemory1_choiceMemory2 = 0;
end
if if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 3
    if_gAcc_noChoice0_allMemory1_choiceMemory2 = 2;
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

gAcc_choice_collapsed = load_gAcc.gAcc_choice_collapsed;
gAcc_choice_collapsed_inOne = [];
for tempi=1:length(gAcc_choice_collapsed)
    gAcc_choice_collapsed_inOne = [gAcc_choice_collapsed_inOne gAcc_choice_collapsed{tempi}'];
end

if if_gAcc_noChoice0_allMemory1_choiceMemory2 == 0
    gAcc_target_collapsed_inOne = gAcc_noChoice_collapsed_inOne;
    gAcc_target_collapsed = gAcc_noChoice_collapsed;
elseif if_gAcc_noChoice0_allMemory1_choiceMemory2 == 1
    gAcc_target_collapsed_inOne = gAcc_allMemory_collapsed_inOne;
    gAcc_target_collapsed = gAcc_allMemory_collapsed;
elseif if_gAcc_noChoice0_allMemory1_choiceMemory2 == 2
    gAcc_target_collapsed_inOne = gAcc_choice_collapsed_inOne;
    gAcc_target_collapsed = gAcc_choice_collapsed;        
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

%% if_usePureNeuronDecoding
neuronBoolIndex = true(size(F_dff_decisionPeriodA,1),1); %#ok<*PREALL>

if if_usePureNeuronDecoding == 1
    neuronBoolIndex = neuronBoolIndex2;
end


%% Preparation
if if_dff_singleSession1_twoSession2_allSession3 < 3
    % F_dff_decisionBin1 = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(2):decisionPeriod_interval(3)),3);
    F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
    F_dff_decisionBin1 = double(F_dff_decisionBin1);
    F_dff_decisionBin1 = F_dff_decisionBin1 + eps;
        
    F_dff_decisionBin1 = F_dff_decisionBin1(neuronBoolIndex,:);
    
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
        
        
        if if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 0
            temp_trialIndex_bool = trialIndex_bool_memoryCorrect; % use memory correct trial
        elseif if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 1
            temp_trialIndex_bool = trial_para_ifSelectOffloading==-1; % use all memory trial
        elseif if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 2
            temp_trialIndex_bool = trialIndex_bool_memoryCorrect & (trial_para_choiceCondition_flag==0); % use noChoice correct trial
        elseif if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 3
            temp_trialIndex_bool = trialIndex_bool_memoryCorrect & (trial_para_choiceCondition_flag==2); % use choiceMemory correct trial            
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

%% Train decoder
if if_compute_train == 1
    t_fun_seqProbSVM_6location = tic;
    
    if if_dff_singleSession1_twoSession2_allSession3 < 3
        if if_seqLevel0_trialLevel1 == 0
            % abandoned now
        elseif if_seqLevel0_trialLevel1 == 1
            if if_decoder_6location0_56seq1 == 0
                
                fun_wrapped_locSVM_Name_v = autoGetFunName_myScripts('fun_wrapped_locSVM', [targetPATH '\functions']);
                fun_wrapped_locSVM = str2func(fun_wrapped_locSVM_Name_v);
                
                options_wrapped_locSVM = struct;
                
                options_wrapped_locSVM.if_normal0_resample1_leaveOneSeqOut2 = if_normal0_resample1_leaveOneSeqOut2;
                options_wrapped_locSVM.resampleTrialCount = resampleTrialCount;
                options_wrapped_locSVM.resampleIterCount = resampleIterCount;
                options_wrapped_locSVM.if_resample_seqBased0_locBased1 = if_resample_seqBased0_locBased1;
                
                options_wrapped_locSVM.numSeq = numSeq;
                options_wrapped_locSVM.numFrames = numFrames;
                options_wrapped_locSVM.targetPATH = targetPATH;
                
                options_wrapped_locSVM.boolIndex_location_seq = boolIndex_location_seq;
                options_wrapped_locSVM.boolIndex_location_seq_T = boolIndex_location_seq_T;
                options_wrapped_locSVM.boolIndex_location_trial = boolIndex_location_trial;
                options_wrapped_locSVM.seqCount = seqCount;
                options_wrapped_locSVM.seqIndex_valid = seqIndex_valid;
                options_wrapped_locSVM.temp_F_dff_decisionBin = temp_F_dff_decisionBin;
                
                options_wrapped_locSVM.gAcc_target_collapsed_inOne = gAcc_target_collapsed_inOne;
                options_wrapped_locSVM.KFold_num = KFold_num;
                options_wrapped_locSVM.if_locDecoder_meanPosterior0_meanAcc1 = if_locDecoder_meanPosterior0_meanAcc1;
                options_wrapped_locSVM.t_decoder = t_decoder;
                options_wrapped_locSVM.if_decodingAcc_threshold0_sort1 = if_decodingAcc_threshold0_sort1;
                
                output_wrapped_locSVM = fun_wrapped_locSVM(options_wrapped_locSVM);
                
                svm_train_length1_outputs = output_wrapped_locSVM.svm_train_length1_outputs;
                svm_train_length2_outputs = output_wrapped_locSVM.svm_train_length2_outputs;
                svm_train_length3_outputs = output_wrapped_locSVM.svm_train_length3_outputs;
                
                temp_trialIndex_valid_resample = output_wrapped_locSVM.temp_trialIndex_valid_resample;
                
                a = 1;
                
            elseif if_decoder_6location0_56seq1 == 1
                % abandoned now
            end
        end
        %fprintf('t_fun_seqProbSVM_6location = %.1f secs.\n',toc(t_fun_seqProbSVM_6location));
        
    elseif if_dff_singleSession1_twoSession2_allSession3 == 3
        %% Get resampled SVM
        fun_resampled_svm_train_lengthx_outputs_Name_v = autoGetFunName_myScripts('fun_resampled_svm_train_lengthx_outputs', [targetPATH '\functions']);
        fun_resampled_svm_train_lengthx_outputs = str2func(fun_resampled_svm_train_lengthx_outputs_Name_v);
        
        options_resampledSVM = struct;
        options_resampledSVM.numSeq = numSeq;
        options_resampledSVM.gAcc_target_collapsed_inOne = gAcc_target_collapsed_inOne;
        options_resampledSVM.seqIndex_valid_resampled = seqIndex_valid_resampled;
        options_resampledSVM.temp_F_dff_decisionBin_resampled = temp_F_dff_decisionBin_resampled;
        options_resampledSVM.boolIndex_location_seq = boolIndex_location_seq;
        options_resampledSVM.KFold_num = KFold_num;
        options_resampledSVM.resampleIterCount = resampleIterCount;
        options_resampledSVM.numFrames = numFrames;
        options_resampledSVM.if_locDecoder_meanPosterior0_meanAcc1 = if_locDecoder_meanPosterior0_meanAcc1;
        options_resampledSVM.targetPATH = targetPATH;
        options_resampledSVM.t_decoder = t_decoder;
        options_resampledSVM.if_normal0_resample1_leaveOneSeqOut2 = if_normal0_resample1_leaveOneSeqOut2;
        options_resampledSVM.if_decodingAcc_threshold0_sort1 = if_decodingAcc_threshold0_sort1;
        
        output_resampledSVM = fun_resampled_svm_train_lengthx_outputs(options_resampledSVM);
        
        svm_train_length1_outputs = output_resampledSVM.svm_train_length1_outputs;
        svm_train_length2_outputs = output_resampledSVM.svm_train_length2_outputs;
        svm_train_length3_outputs = output_resampledSVM.svm_train_length3_outputs;
        
        a = 1;
    end
elseif if_compute_train == 0
    
end
a = 1;
% rand(1)
a = 1;

% %% Compute seq-level mean location distribution
% for target_length=1:3
%     if target_length == 2
%         temp_svm = svm_train_length2_outputs;
%     elseif target_length == 3
%         temp_svm = svm_train_length3_outputs;
%     elseif target_length == 1
%         temp_svm = svm_train_length1_outputs;
%     end
%
%     temp_trialNum = size(temp_svm.temp_svm_Y_valid_T,1);
%     temp_seqIndex_trial = zeros(temp_trialNum,1);
%     for tempi=1:temp_trialNum
%         temp_seqBoolIndex_trial = temp_svm.temp_svm_Y_valid_T(tempi,:);
%         for tempj=1:numSeq(target_length)
%             temp_seqIndex = temp_svm.seq_range(tempj);
%             temp_seqBoolIndex = temp_svm.boolIndex_location_seq(:,temp_seqIndex)';
%
%             if sum(temp_seqBoolIndex==temp_seqBoolIndex_trial) == numFrames
%                 temp_seqIndex_trial(tempi) = temp_svm.seq_range(tempj);
%                 break
%             end
%         end
%     end
%     temp_svm.seqIndex_trial = temp_seqIndex_trial;
%
%
%     Posterior_2d_n11n_lengthx_mean = zeros(numSeq(target_length),numFrames);
%     Posterior_2d_w_lengthx_mean = zeros(numSeq(target_length),numFrames);
%     for tempi=1:numSeq(target_length)
%         temp_seqIndex = temp_svm.seq_range(tempi);
%         tempSeqBoolIndex = temp_svm.boolIndex_location_seq(:,temp_seqIndex)';
%
%         tempBoolIndex = temp_seqIndex_trial == temp_seqIndex;
%         temp_Posterior_2d_n11n = temp_svm.Posterior_2d_n11n(tempBoolIndex,:);
%         temp_Posterior_2d_w = temp_svm.Posterior_2d_w(tempBoolIndex,:);
%
%         if if_locDecoder_meanPosterior0_meanAcc1 == 0
%             Posterior_2d_n11n_lengthx_mean(tempi,:) = mean(temp_Posterior_2d_n11n,1);
%             Posterior_2d_w_lengthx_mean(tempi,:) = mean(temp_Posterior_2d_w,1);
%         elseif if_locDecoder_meanPosterior0_meanAcc1 == 1
%             tempBoolIndex2 = temp_Posterior_2d_n11n>0.5;
%             for tempj=1:numFrames
%                 tempBoolIndex3 = tempBoolIndex2(:,tempj) == tempSeqBoolIndex(tempj);
%                 Posterior_2d_n11n_lengthx_mean(tempi,tempj) = sum(tempBoolIndex3)/length(tempBoolIndex3);
%             end
%             tempBoolIndex2 = temp_Posterior_2d_w>0.5;
%             for tempj=1:numFrames
%                 tempBoolIndex3 = tempBoolIndex2(:,tempj) == tempSeqBoolIndex(tempj);
%                 Posterior_2d_w_lengthx_mean(tempi,tempj) = sum(tempBoolIndex3)/length(tempBoolIndex3);
%             end
%         end
%     end
%     temp_svm.Posterior_2d_n11n_lengthx_mean = Posterior_2d_n11n_lengthx_mean;
%     temp_svm.Posterior_2d_w_lengthx_mean = Posterior_2d_w_lengthx_mean;
%
%     if target_length == 2
%         svm_train_length2_outputs = temp_svm;
%     elseif target_length == 3
%         svm_train_length3_outputs = temp_svm;
%     elseif target_length == 1
%         svm_train_length1_outputs = temp_svm;
%     end
% end

%% Compute seq-level probability production
coeff_w_power_length1 = svm_train_length1_outputs.coeff_w_power;
coeff_w_power_length2 = svm_train_length2_outputs.coeff_w_power;
coeff_w_power_length3 = svm_train_length3_outputs.coeff_w_power;

for target_length=1:valid_length
    if target_length == 1
        svm_train_lengthx_outputs = svm_train_length1_outputs;
    elseif target_length == 2
        svm_train_lengthx_outputs = svm_train_length2_outputs;
    elseif target_length == 3
        svm_train_lengthx_outputs = svm_train_length3_outputs;
    end
    Posterior_2d_n11n_lengthx_mean = svm_train_lengthx_outputs.Posterior_2d_n11n_lengthx_mean;
    seq_range = svm_train_lengthx_outputs.seq_range;
    temp_boolIndex_location_seq = boolIndex_location_seq(:,seq_range)';
    temp_p = Posterior_2d_n11n_lengthx_mean;
    temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
    temp_p_seq = prod(temp_p,2);
    svm_train_lengthx_outputs.p_seq_prod = temp_p_seq';
    
    if target_length == 1
        svm_train_length1_outputs = svm_train_lengthx_outputs;
    elseif target_length == 2
        svm_train_length2_outputs = svm_train_lengthx_outputs;
    elseif target_length == 3
        svm_train_length3_outputs = svm_train_lengthx_outputs;
    end
    
end

if if_decoder_decodingAcc0_pProd1 == 0
    svm_posterior_length1 = svm_train_length1_outputs.svm_posterior_lengthx;
    svm_posterior_length2 = svm_train_length2_outputs.svm_posterior_lengthx;
    svm_posterior_length3 = svm_train_length3_outputs.svm_posterior_lengthx;
elseif if_decoder_decodingAcc0_pProd1 == 1
    svm_posterior_length1 = svm_train_length1_outputs.p_seq_prod;
    svm_posterior_length2 = svm_train_length2_outputs.p_seq_prod;
    svm_posterior_length3 = svm_train_length3_outputs.p_seq_prod;
end
fprintf('if_decoder_decodingAcc0_pProd1=%d.\n',if_decoder_decodingAcc0_pProd1);

svm_posterior_length1_mean = mean(svm_posterior_length1,'omitnan');
svm_posterior_length2_mean = mean(svm_posterior_length2,'omitnan');
svm_posterior_length3_mean = mean(svm_posterior_length3,'omitnan');

%%
if if_dff_singleSession1_twoSession2_allSession3 < 3
    trial_num = length(trialIndex_bool_memoryCorrect);
    % trialIndex_bool_allMemory = ~(trial_para.ifSelectOffloading(1:trial_para.trial_count)==1);
    trialIndex_bool_allMemory = ~(trial_para_ifSelectOffloading(1:trial_num)==1);
    trialIndex_bool_memoryError = trialIndex_bool_allMemory & (~trialIndex_bool_memoryCorrect);
    
    F_dff_decisionBin1;
    
    %gAcc_allLength = cell(1,3);
    seqCount_allLength_memoryCorrect = cell(1,4);
    seqCount_allLength_allMemory = cell(1,4);
    seqCount_allLength_gAcc = cell(1,4);
    seqCount_allLength_memoryError = cell(1,4);
    for target_length=1:4
        
        temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
        
        temp_seqCount_memoryCorrect = zeros(1,numSeq(target_length));
        temp_seqCount_allMemory = zeros(1,numSeq(target_length));
        temp_seqCount_memoryError = zeros(1,numSeq(target_length));
        for tempi=1:numSeq(target_length)
            tempj = temp_range(tempi);
            temp_seqCount_memoryCorrect(tempi) = sum(seqIndex(temp_trialIndex_bool)==tempj,'all');
            if if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 0
                temp_seqCount_allMemory(tempi) = sum(seqIndex(trialIndex_bool_allMemory)==tempj,'all');
            elseif if_trainingSet_memoryCorrect0_allMemory1_noChoiceCorrect2_CMC3 == 2
                temp_seqCount_allMemory(tempi) = sum(seqIndex(trial_para_choiceCondition_flag==0)==tempj,'all');
            end
            temp_seqCount_memoryError(tempi) = sum(seqIndex(trialIndex_bool_memoryError)==tempj,'all');
        end
        
        seqCount_allLength_memoryCorrect{target_length} = temp_seqCount_memoryCorrect;
        seqCount_allLength_allMemory{target_length} = temp_seqCount_allMemory;
        seqCount_allLength_gAcc{target_length} = temp_seqCount_memoryCorrect./temp_seqCount_allMemory;
        %gAcc_allLength{target_length} = gAcc;
        
        seqCount_allLength_memoryError{target_length} = temp_seqCount_memoryError;
    end
end

%% Compute correlation with behavior

if if_gAcc_behavior0_model1 == 0
    if if_gAcc_behavior_singleSession0_allSession1 == 1
        gAcc_length1 = gAcc_target_collapsed{1}';
        gAcc_length2 = gAcc_target_collapsed{2}';
        gAcc_length3 = gAcc_target_collapsed{3}';
    elseif if_gAcc_behavior_singleSession0_allSession1 == 0
        gAcc_length1 = seqCount_allLength_gAcc{1};
        gAcc_length2 = seqCount_allLength_gAcc{2};
        gAcc_length3 = seqCount_allLength_gAcc{3};
        
        gAcc_target_collapsed = seqCount_allLength_gAcc;
        gAcc_target_collapsed_inOne = ...
            [seqCount_allLength_gAcc{1},...
            seqCount_allLength_gAcc{2},...
            seqCount_allLength_gAcc{3},...
            seqCount_allLength_gAcc{4}];
        
        
    end
elseif if_gAcc_behavior0_model1 == 1
    gAcc_length1 = gAcc_noChoice_model{1}';
    gAcc_length2 = gAcc_noChoice_model{2}';
    gAcc_length3 = gAcc_noChoice_model{3}';
end


a = 1;
% rand(1)
a = 1;
% svm_posterior_length1
x = [svm_posterior_length1';svm_posterior_length2(~isnan(svm_posterior_length2))';svm_posterior_length3(~isnan(svm_posterior_length3))'];
y = [gAcc_length1';gAcc_length2(~isnan(svm_posterior_length2))';gAcc_length3(~isnan(svm_posterior_length3))'];
[r_123,p_123] = corr(x,y);
if if_gAcc_behavior0_model1 == 0
    fprintf('r_length123_decoder_behavior=%.4f, p_length123_decoder_behavior=%.4f.\n',r_123,p_123);
elseif if_gAcc_behavior0_model1 == 1
    fprintf('r_length123_decoder_model   =%.4f, p_length123_decoder_model   =%.4f.\n',r_123,p_123);
end

[r_1,p_1] = corr(svm_posterior_length1',gAcc_length1');
[r_2,p_2] = corr(svm_posterior_length2(~isnan(svm_posterior_length2))',gAcc_length2(~isnan(svm_posterior_length2))');
[r_3,p_3] = corr(svm_posterior_length3(~isnan(svm_posterior_length3))',gAcc_length3(~isnan(svm_posterior_length3))');
fprintf('r_1=%.4f, p_1=%.4f, mean=%.4f.\n',r_1,p_1,mean(svm_posterior_length1')); %#ok<*UDIM>
fprintf('r_2=%.4f, p_2=%.4f, mean=%.4f.\n',r_2,p_2,mean(svm_posterior_length2(~isnan(svm_posterior_length2))'));
fprintf('r_3=%.4f, p_3=%.4f, mean=%.4f.\n',r_3,p_3,mean(svm_posterior_length3(~isnan(svm_posterior_length3))'));



Posterior_2d_n11n_mean = [svm_train_length1_outputs.Posterior_2d_n11n_lengthx_mean;...
    svm_train_length2_outputs.Posterior_2d_n11n_lengthx_mean;...
    svm_train_length3_outputs.Posterior_2d_n11n_lengthx_mean];

Posterior_2d_w_mean = [svm_train_length1_outputs.Posterior_2d_w_lengthx_mean;...
    svm_train_length2_outputs.Posterior_2d_w_lengthx_mean;...
    svm_train_length3_outputs.Posterior_2d_w_lengthx_mean];


temp_range = 1:sum(numSeq(1:valid_length));
Posterior_2d_model = FittingResults.intermediate.locationDistri_seq(temp_range,:);

Posterior_2d_model_raw = Posterior_2d_model;

r_n11n = zeros(1,length(temp_range));
p_n11n = zeros(1,length(temp_range));
for tempi=1:length(temp_range)
    [r_n11n(tempi),p_n11n(tempi)] = corr(Posterior_2d_n11n_mean(tempi,:)',Posterior_2d_model(tempi,:)');
end


r_w = zeros(1,length(temp_range));
p_w = zeros(1,length(temp_range));
for tempi=1:length(temp_range)
    [r_w(tempi),p_w(tempi)] = corr(Posterior_2d_w_mean(tempi,:)',Posterior_2d_model(tempi,:)');
end

% fprintf('num(p_n11n<0.05)=%d, sum(p_w<0.05)=%d.\n',sum(p_n11n<0.05),sum(p_w<0.05));
fprintf('num(p_posterior_seq_n11n<0.05)=%d/%d.\n',sum(p_n11n<0.05),sum(~isnan(p_n11n)));


r_loc_n11n = zeros(1,numFrames);
p_loc_n11n = zeros(1,numFrames);
for tempi=1:numFrames
    tempBoolIndex = ~isnan(Posterior_2d_n11n_mean(:,1));
    [r_loc_n11n(tempi),p_loc_n11n(tempi)] = ...
        corr(Posterior_2d_n11n_mean(tempBoolIndex,tempi),Posterior_2d_model(tempBoolIndex,tempi));
end

r_loc_w = zeros(1,numFrames);
p_loc_w = zeros(1,numFrames);
for tempi=1:numFrames
    tempBoolIndex = ~isnan(Posterior_2d_w_mean(:,1));
    [r_loc_w(tempi),p_loc_w(tempi)] = ...
        corr(Posterior_2d_w_mean(tempBoolIndex,tempi),Posterior_2d_model(tempBoolIndex,tempi));
end

% fprintf('r_loc_n11n= %s.\n',num2str(round(r_loc_n11n*10000)/10000));
% fprintf('r_loc_w=    %s.\n',num2str(round(r_loc_w*10000)/10000));
% fprintf('r_posterior_loc_n11n= %s.\n',num2str(round(r_loc_n11n*10000)/10000));

locMSE_n11n = mean((Posterior_2d_n11n_mean-Posterior_2d_model).^2,1,'omitnan');
locMSE_w = mean((Posterior_2d_w_mean-Posterior_2d_model).^2,1,'omitnan');
% fprintf('locMSE_n11n=%s.\n',num2str(round(locMSE_n11n*10000)/10000));
% fprintf('locMSE_w=   %s.\n',num2str(round(locMSE_w*10000)/10000));

a = 1;

%% get beta_loc_lengthx
% if if_normal0_resample1_leaveOneSeqOut2 == 0
%     KFold_num = length(svm_train_length1_outputs.temp_Mdl_CV_binary{tempi}.Trained);
%     beta_loc_lengthx = zeros(roiNum,numFrames);
%     for tempLength=1:valid_length
%         if tempLength == 1
%             temp_svm = svm_train_length1_outputs;
%         elseif tempLength == 2
%             temp_svm = svm_train_length2_outputs;
%         elseif tempLength == 3
%             temp_svm = svm_train_length3_outputs;
%         end
%         for tempi=1:numFrames
%             temp_beta = zeros(roiNum,KFold_num);
%             for tempk=1:KFold_num
%                 temp_beta(:,tempk) = temp_svm.temp_Mdl_CV_binary{tempi}.Trained{tempk}.BinaryLearners{1}.Beta;
%             end
%             temp_beta_median = median(temp_beta,2);
%             beta_loc_lengthx(:,tempi) = temp_beta_median;
%         end
%         if tempLength == 1
%             beta_loc_length1 = beta_loc_lengthx;
%         elseif tempLength == 2
%             beta_loc_length2 = beta_loc_lengthx;
%         elseif tempLength == 3
%             beta_loc_length3 = beta_loc_lengthx;
%         end
%     end
%     beta_loc_length1;
%     beta_loc_length2;
%     beta_loc_length3;
%
%     topX = 5;
%
%     for tempLength=1:valid_length
%         if tempLength == 1
%             beta_loc_lengthx = beta_loc_length1;
%         elseif tempLength == 2
%             beta_loc_lengthx = beta_loc_length2;
%         elseif tempLength == 3
%             beta_loc_lengthx = beta_loc_length3;
%         end
%
%         beta_loc_lengthx_topXIndex = zeros(topX,numFrames);
%         for tempi=1:numFrames
%             [B,I] = sort(beta_loc_lengthx(:,tempi),'descend');
%             beta_loc_lengthx_topXIndex(:,tempi) = I(1:topX);
%         end
%
%         if tempLength == 1
%             beta_loc_length1_topXIndex = beta_loc_lengthx_topXIndex;
%         elseif tempLength == 2
%             beta_loc_length2_topXIndex = beta_loc_lengthx_topXIndex;
%         elseif tempLength == 3
%             beta_loc_length3_topXIndex = beta_loc_lengthx_topXIndex;
%         end
%     end
%
%     svm_beta = struct;
%     svm_beta.beta_loc_length1 = beta_loc_length1;
%     svm_beta.beta_loc_length2 = beta_loc_length2;
%     svm_beta.beta_loc_length3 = beta_loc_length3;
%     svm_beta.beta_loc_length1_topXIndex = beta_loc_length1_topXIndex;
%     svm_beta.beta_loc_length2_topXIndex = beta_loc_length2_topXIndex;
%     svm_beta.beta_loc_length3_topXIndex = beta_loc_length3_topXIndex;
%     svm_beta.topX = topX;
%     svm_beta.roiNum = roiNum;
%
% end


fprintf('t_script = %.1f secs.\n',toc(t0));

a = 1;
% rand(1)
a = 1;


if if_plot == 1
    test_plot_r123_r_loc_Name_v = autoGetFunName_myScripts('test_plot_r123_r_loc', [targetPATH]);
    script_test_plot_r123_r_loc = str2func(test_plot_r123_r_loc_Name_v);
    script_test_plot_r123_r_loc();
end

%% Precision, sigma-based
test_sigmaBased_precision_Name_v = autoGetFunName_myScripts('test_temp_sigmaBased_precision', [targetPATH]);
script_test_sigmaBased_precision = str2func(test_sigmaBased_precision_Name_v);
script_test_sigmaBased_precision();




%% End