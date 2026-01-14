% Chuan's 19th script (20251214)
% This script: To conduct trial parcellation of two-photon data, and then save them.
%% Initialization
% clear
close all
% home;% To scroll down in command window

if_compute = 1;

if if_compute == 1
    %clear
    if_compute = 1;
    if_load = 1;
else
    if_load = 0;
end

% if_trialStartDetrend = 1;

% if_load = 0;

if_profile = 0;

if_save_decoding_and_PTB_data = 0;
if_save_glmData = 0;

if_plotSingleCell = 1;
ifSave_fig = 0;

if_plotPCA = 0;
order_glm = 0;% 0 means whole seq

if_compute_selectiveCellNum_decisionBin_lengthxSeq_MON = 0;

if_compute_fun_rProb_glm = 1;
if_memoryPrecision_accuracy0_sigma1 = 0;

if_plot_svm_beta = 0;
if_plot_memoryMetaMismatch = 0;
if_plot_rProb_glm = 0;

plot_lengthFlag = 1;
if_plot_location0_seq1 = 0;

if_substractBaseLine = 0;%1-->0
if_delay1Selective_old0_new1 = 1;% old is T + delay1(0-->1100), new is pure delay1(ChoiceCue-1100-->ChoiceCue)
if_selective_seq0_loc1 = 1;%1

if_plotSelectiveCell = 1;
significantThreshold = 0.05;%0.01-->0.05-->0.15-->0.001-->0.03-->0.05
significantThreshold_lowBound = -1;
%selectiveCellIndex_suite2p_lengthx

if_sortRasterPlot_delay1 = 0;
% if_sortRasterPlotSeq_decision = 1;
% if_plot_MNOtrial0_MNOseq1_afterMtrial2 = 2;
if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 = 1;%2
if_plotSelectiveCell_0T_1delay1_2delay2_3MNO = 1;


plotRoiNum = 0;
temp_id_suite2p = 0;

temp_id = temp_id_suite2p + 1;

windowSize = 5;%3-->5
% compensateImageDelay = 1;%1-->8-->3, to compensate imaging delay

compensateImageDelay = 0;
% compensateImageDelay = compensateImageDelay + 1;


currentSession_multi = string;

% currentSession_multi = [currentSession_multi; '113Recording_20240111A_Zelku_Site09A'];
% currentSession_multi = [currentSession_multi; '113Recording_20240112A_Zelku_Site06A'];
% currentSession_multi = [currentSession_multi; '113Recording_20240115A_Zelku_Site06A'];
% currentSession_multi = [currentSession_multi; '113Recording_20240117A_Zelku_Site14A'];
% currentSession_multi = [currentSession_multi; '113Recording_20240118A_Zelku_Site18A'];
% currentSession_multi = [currentSession_multi; '113Recording_20240119A_Zelku_Site17A'];
% currentSession_multi = [currentSession_multi; '113Recording_20240122A_Zelku_Site09B'];
% currentSession_multi = [currentSession_multi; '113Recording_20240123A_Zelku_Site09B_sameFOV0122'];
% currentSession_multi = [currentSession_multi; '113Recording_20240124A_Zelku_Site06B'];
% currentSession_multi = [currentSession_multi; '113Recording_20240126A_Zelku_Site06B_sameFOV0124'];
% currentSession_multi = [currentSession_multi; '113Recording_20240129A_Zelku_Site07A'];
% currentSession_multi = [currentSession_multi; '113Recording_20240131A_Zelku_Site07A_sameFOV0129'];
% currentSession_multi = [currentSession_multi; '113Recording_20240202A_Zelku_Site06XA'];
% currentSession_multi = [currentSession_multi; '113Recording_20240203A_Zelku_Site06XA_sameFOV0202'];
% currentSession_multi = [currentSession_multi; '113Recording_20240207A_Zelku_Site05A'];
% currentSession_multi = [currentSession_multi; '113Recording_20240208A_Zelku_Site05A_sameFOV0207'];
% currentSession_multi = [currentSession_multi; '113Recording_20240210A_Zelku_Site10A'];
% currentSession_multi = [currentSession_multi; '113Recording_20240211A_Zelku_Site10A_sameFOV0210'];
% currentSession_multi = [currentSession_multi; '113Recording_20240216A_Zelku_Site09C'];
% currentSession_multi = [currentSession_multi; '113Recording_20240218A_Zelku_Site09C_sameFOV0216'];
% currentSession_multi = [currentSession_multi; '113Recording_20240220A_Zelku_Site06XB'];
% currentSession_multi = [currentSession_multi; '113Recording_20240221A_Zelku_Site06XB_sameFOV0220'];
% currentSession_multi = [currentSession_multi; '113Recording_20240226A_Zelku_Site10B'];
% currentSession_multi = [currentSession_multi; '113Recording_20240227A_Zelku_Site10B_sameFOV0226'];
% currentSession_multi = [currentSession_multi; '113Recording_20240229A_Zelku_Site06C'];
% currentSession_multi = [currentSession_multi; '113Recording_20240301A_Zelku_Site06C_sameFOV0229'];
% currentSession_multi = [currentSession_multi; '113Recording_20240304A_Zelku_Site09D'];
% currentSession_multi = [currentSession_multi; '113Recording_20240305A_Zelku_Site09D_sameFOV0304'];
% currentSession_multi = [currentSession_multi; '113Recording_20240307A_Zelku_Site10C'];
% currentSession_multi = [currentSession_multi; '113Recording_20240308A_Zelku_Site10C_sameFOV0307'];
% currentSession_multi = [currentSession_multi; '113Recording_20240312A_Zelku_Site06RA'];
currentSession_multi = [currentSession_multi; '113Recording_20240315A_Zelku_Site06RA_sameFOV0312'];
% currentSession_multi = [currentSession_multi; '113Recording_20240319A_Zelku_Site09E'];
% currentSession_multi = [currentSession_multi; '113Recording_20240320A_Zelku_Site09E_sameFOV0319'];
% currentSession_multi = [currentSession_multi; '113Recording_20240322A_Zelku_Site07B'];
% currentSession_multi = [currentSession_multi; '113Recording_20240323A_Zelku_Site07B_sameFOV0322'];
% currentSession_multi = [currentSession_multi; '113Recording_20240329A_Zelku_Site05B'];
% currentSession_multi = [currentSession_multi; '113Recording_20240330A_Zelku_Site05B_sameFOV0329'];
% currentSession_multi = [currentSession_multi; '113Recording_20240402A_Zelku_Site14B'];
% currentSession_multi = [currentSession_multi; '113Recording_20240403A_Zelku_Site14B_sameFOV0402'];
% currentSession_multi = [currentSession_multi; '113Recording_20240410A_Zelku_Site17B'];
% currentSession_multi = [currentSession_multi; '113Recording_20240411A_Zelku_Site17B_sameFOV0410'];


% currentSession_multi = [currentSession_multi; '113Recording_20230426A_Ding_Site16'];
% currentSession_multi = [currentSession_multi; '113Recording_20230427A_Ding_Site16_sameFOV0426'];
% currentSession_multi = [currentSession_multi; '113Recording_20230502A_Ding_Site13'];
% currentSession_multi = [currentSession_multi; '113Recording_20230503A_Ding_Site13_sameFOV0502'];
% currentSession_multi = [currentSession_multi; '113Recording_20230504A_Ding_Site02'];
% currentSession_multi = [currentSession_multi; '113Recording_20230508A_Ding_Site02_sameFOV0509'];
% currentSession_multi = [currentSession_multi; '113Recording_20230509A_Ding_Site02'];% 660000 frames, easy to crash
% % % % 
% currentSession_multi = [currentSession_multi; '113Recording_20230510A_Ding_Site05_sameFOV0511'];
% currentSession_multi = [currentSession_multi; '113Recording_20230510B_Ding_Site05_sameFOV0511'];
% currentSession_multi = [currentSession_multi; '113Recording_20230511A_Ding_Site05'];
% currentSession_multi = [currentSession_multi; '113Recording_20230512A_Ding_Site09'];
% currentSession_multi = [currentSession_multi; '113Recording_20230513A_Ding_Site09_sameFOV0512'];
% 
% currentSession_multi = [currentSession_multi; '113Recording_20230515A_Ding_Site24_sameFOV0516'];
% currentSession_multi = [currentSession_multi; '113Recording_20230516A_Ding_Site24'];
% currentSession_multi = [currentSession_multi; '113Recording_20230517A_Ding_Site16B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230522A_Ding_Site05B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230525A_Ding_Site09B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230526A_Ding_Site09B_sameFOV0525'];
% % 
% currentSession_multi = [currentSession_multi; '113Recording_20230527A_Ding_Site02B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230528A_Ding_Site02B_sameFOV0527'];
% currentSession_multi = [currentSession_multi; '113Recording_20230530A_Ding_Site05C'];
% currentSession_multi = [currentSession_multi; '113Recording_20230531A_Ding_Site05C_sameFOV0530'];
% currentSession_multi = [currentSession_multi; '113Recording_20230601A_Ding_Site13B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230602A_Ding_Site13B_sameFOV0601'];
% % 
% currentSession_multi = [currentSession_multi; '113Recording_20230604A_Ding_Site07'];
% currentSession_multi = [currentSession_multi; '113Recording_20230605A_Ding_Site07_sameFOV0604'];
% currentSession_multi = [currentSession_multi; '113Recording_20230612A_Ding_Site14'];
% currentSession_multi = [currentSession_multi; '113Recording_20230614A_Ding_Site05D'];
% currentSession_multi = [currentSession_multi; '113Recording_20230615A_Ding_Site05D_sameFOV0614'];
% currentSession_multi = [currentSession_multi; '113Recording_20230619A_Ding_Site02C'];
% currentSession_multi = [currentSession_multi; '113Recording_20230620A_Ding_Site05E'];

currentSession_multi(1) = [];
num_FOV = length(currentSession_multi);

time0 = tic;
for tempSessionIndex=1:num_FOV
    currentSession = currentSession_multi{tempSessionIndex};
    fprintf('currentSession = %s.\n',currentSession);
    
    
    if if_profile == 1
        profile on
    end
    t0 = tic;
    
    
    if contains(currentSession,'Ding')
        searchName_gAcc = 'from23-04-26to23-06-20_D_gAcc_1';
        searchName_rProb = 'from23-04-26to23-06-20_D_offloadingProb_1';
    else
        %searchName_gAcc = 'from01-08to01-23_Z_gAcc_1';
        %searchName_rProb = 'from01-08to01-23_Z_offloadingProb_1';    
        searchName_gAcc = 'from01-08to01-31_Z_gAcc_1';
        searchName_rProb = 'from01-08to01-31_Z_offloadingProb_1';            
    end
    
    
    
    targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
    cd(targetPATH)
    
    % output_shortPath = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected';
    output_shortPath = 'D:\twoPhotonData_motionCorrected';
    temp_currentSession_path = [output_shortPath '\' currentSession];
    temp_if_max0_min1 = 0;
    output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
    path_plane = [output_path,'\plane0'];
    
    if if_load == 1
        load([output_path,'\','markerParse_trialLevel.mat'],'markerParse_trialLevel');
    end
    
    output_path_singleCell = [output_path,'\singleCell'];
    if ~exist(output_path_singleCell, 'dir')
        mkdir(output_path_singleCell);
    end
    
    
    % rawData_path = 'D:\twoPhotonRawData\ToBeProcessed';
    rawData_path = 'R:\twoPhotonRawData\ToBeProcessed';
    currentSession_path = [rawData_path '\' currentSession];
    PTB_name = '.mat';
    PTB_fullName = autoGetMarkerFileName(['2*',PTB_name], currentSession_path);
    
    if if_load == 1
        basic_para = [];
        trial_para = [];
        load(PTB_fullName,'basic_para','trial_para');
        a = 1;
        
        %path_behavior = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\behavior';
        path_behavior = [output_shortPath '\behavior'];
        
        % Load other processed results
        load_gAcc = loadMat_single(searchName_gAcc, path_behavior);
        gAcc_noChoice_collapsed_inOne = load_gAcc.gAcc_noChoice_collapsed_inOne;
        seqPrecision_behavior = load_gAcc.seqPrecision_behavior;
        
        
        if if_memoryPrecision_accuracy0_sigma1 == 0
            precision_inOne = gAcc_noChoice_collapsed_inOne;
        elseif if_memoryPrecision_accuracy0_sigma1 == 1
            precision_inOne = seqPrecision_behavior';
        end
        
        % Load other processed results
        load_rProb = loadMat_single(searchName_rProb, path_behavior);
        offloadingProb_collapsed = load_rProb.offloadingProb_all;
        
        offloadingProb_inOne = [];
        pointKindsNum = 4;
        for tempi=1:pointKindsNum
            offloadingProb_inOne = [offloadingProb_inOne offloadingProb_collapsed{tempi}'];
        end
        
        
        cd(targetPATH)
        
        fileName_Fall = 'Fall.mat';
        fileName_iscell = 'iscell.npy';
        fullFileName_Fall = [path_plane,'\',fileName_Fall];
        fullFileName_iscell = [path_plane,'\',fileName_iscell];
        
        % load(fullFileName_Fall,'F_dff');
        load(fullFileName_Fall,'F_dff_raw');
        load(fullFileName_Fall,'F0_raw');
        iscell = readNPY(fullFileName_iscell);
        
        % iscell(:,1) = 1;
        
        F_dff = F_dff_raw(iscell(:,1)==1, :);
        F0 = F0_raw(iscell(:,1)==1, :);
        clear F_dff_raw
        clear F0_raw
        
        s = load(fullFileName_Fall,'stat');
        roi_stats_raw = s.stat;
        temp_cellIndex = find(iscell(:,1)==1);
        roi_stats = roi_stats_raw(temp_cellIndex);
        
        markerParse_trialLevel;
        iscell;
        F_dff;
        % F_dff_raw;
        
        % windowSize = 5;%3-->5
        x = F_dff;
        F_dff = smoothdata(x,2,'gaussian',windowSize);
        % F_dff = smoothdata(x,2,'sgolay',windowSize,'Degree',4);
        % F_dff = smoothdata(x,2,'rloess',windowSize);
        % F_dff = halfGaussianFilter_v1(F_dff,windowSize);
        
        % F_dff(:,1:end-1) = F_dff(:,2:end); % to compensate imaging delay
        
        % compensateImageDelay = 3;%1-->8-->3, to compensate imaging delay
        F_dff(:,1:end-compensateImageDelay) = F_dff(:,1+compensateImageDelay:end); % to compensate imaging delay
        
        
        F0 = smoothdata(F0,2,'gaussian',windowSize);
        F0(:,1:end-compensateImageDelay) = F0(:,1+compensateImageDelay:end); % to compensate imaging delay
        
        
        temp_if_max0_min1 = 0;
        %template_path = autoGetFileName_general('template*.tif', output_path,temp_if_max0_min1);
        template_path = autoGetFileName_general('maxProjection*.tif', output_path,temp_if_max0_min1);        
        template = double(loadtiff(template_path));
    end
    
    % cellBoolIndex = iscell(:,2)>1;
    % cellIndex = find(iscell(:,2)>1);
    cellIndex = find(iscell(:,1)==1);
    cellIndex_suite2p = cellIndex - 1;
    
    if isempty(markerParse_trialLevel{end}.currentTrialTotalMarkers_frameIndex) == 1
        markerParse_trialLevel(end) = [];
    end
    
    trialIndex_length1 = [];
    trialIndex_length2 = [];
    trialIndex_length3 = [];
    trialIndex_length4 = [];
    
    for tempi=1:length(markerParse_trialLevel)
        temp_seq_length = length(markerParse_trialLevel{tempi}.currentSequence);
        %trial_para.isCorrect()
        if temp_seq_length == 1
            trialIndex_length1 = [trialIndex_length1 tempi]; %#ok<*AGROW>
        elseif temp_seq_length == 2
            trialIndex_length2 = [trialIndex_length2 tempi];
        elseif temp_seq_length == 3
            trialIndex_length3 = [trialIndex_length3 tempi];
        elseif temp_seq_length == 4
            trialIndex_length4 = [trialIndex_length4 tempi];
        end
    end
    
    sequence_length1 = zeros(1,length(trialIndex_length1));
    for tempi=1:length(trialIndex_length1)
        temp_trialIndex = trialIndex_length1(tempi);
        sequence_length1(tempi) = markerParse_trialLevel{temp_trialIndex}.currentSequence;
    end
    
    sequence_length2 = zeros(2,length(trialIndex_length2));
    for tempi=1:length(trialIndex_length2)
        temp_trialIndex = trialIndex_length2(tempi);
        sequence_length2(:,tempi) = markerParse_trialLevel{temp_trialIndex}.currentSequence;
    end
    
    sequence_length3 = zeros(3,length(trialIndex_length3));
    for tempi=1:length(trialIndex_length3)
        temp_trialIndex = trialIndex_length3(tempi);
        sequence_length3(:,tempi) = markerParse_trialLevel{temp_trialIndex}.currentSequence;
    end
    
    sequence_length4 = zeros(4,length(trialIndex_length4));
    for tempi=1:length(trialIndex_length4)
        temp_trialIndex = trialIndex_length4(tempi);
        sequence_length4(:,tempi) = markerParse_trialLevel{temp_trialIndex}.currentSequence;
    end
    
    a = 1;
    numFrames = 6;
    
    responseSeq_length1_onehot = false(length(trialIndex_length1) ,numFrames);
    for tempi=1:length(trialIndex_length1)
        temp_trialIndex = trialIndex_length1(tempi);
        responseSeq_length1_onehot(tempi,:) = ~trial_para.isFilled{temp_trialIndex};
    end
    
    responseSeq_length2_onehot = false(length(trialIndex_length2) ,numFrames);
    for tempi=1:length(trialIndex_length2)
        temp_trialIndex = trialIndex_length2(tempi);
        responseSeq_length2_onehot(tempi,:) = ~trial_para.isFilled{temp_trialIndex};
    end
    
    responseSeq_length3_onehot = false(length(trialIndex_length3) ,numFrames);
    for tempi=1:length(trialIndex_length3)
        temp_trialIndex = trialIndex_length3(tempi);
        responseSeq_length3_onehot(tempi,:) = ~trial_para.isFilled{temp_trialIndex};
    end
    
    responseSeq_length4_onehot = false(length(trialIndex_length4) ,numFrames);
    for tempi=1:length(trialIndex_length4)
        temp_trialIndex = trialIndex_length4(tempi);
        responseSeq_length4_onehot(tempi,:) = ~trial_para.isFilled{temp_trialIndex};
    end
    
    % numFrames = 6;
    pointKindsNum = 4;
    target_seqSet = get_target_seqSet_v1(numFrames,pointKindsNum);
    numSeq = zeros(1,pointKindsNum);
    for tempi=1:pointKindsNum
        numSeq(tempi) = length(target_seqSet{tempi});
    end
    
    seqIndex = zeros(1,trial_para.trial_count);
    target_seqSet;
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
    target_seqSet_inOne = [];
    for tempi=1:size(target_seqSet,2)
        target_seqSet_inOne = [target_seqSet_inOne target_seqSet{tempi}'];
    end
    
    
    %% trialIndex_lengthx_seq
    trialIndex_length1_seq = cell(numFrames,1);
    for tempi=1:numFrames
        temp_boolIndex = (sequence_length1 == tempi);
        temp_index = trialIndex_length1(temp_boolIndex);
        trialIndex_length1_seq{tempi} = temp_index;
    end
    a = 1;
    
    trialIndex_length2_seq = cell(numSeq(2),1);
    for tempi=1:numSeq(2)
        temp_boolIndex = false(1,size(sequence_length2,2));
        for tempj=1:size(sequence_length2,2)
            if sum(ismember(sequence_length2(:,tempj)',target_seqSet{2}{tempi})) == 2
                temp_boolIndex(tempj) = true;
            end
        end
        temp_index = trialIndex_length2(temp_boolIndex);
        trialIndex_length2_seq{tempi} = temp_index;
    end
    
    trialIndex_length3_seq = cell(numSeq(3),1);
    for tempi=1:numSeq(3)
        temp_boolIndex = false(1,size(sequence_length3,2));
        for tempj=1:size(sequence_length3,2)
            if sum(ismember(sequence_length3(:,tempj)',target_seqSet{3}{tempi})) == 3
                temp_boolIndex(tempj) = true;
            end
        end
        temp_index = trialIndex_length3(temp_boolIndex);
        trialIndex_length3_seq{tempi} = temp_index;
    end
    
    trialIndex_length4_seq = cell(numSeq(4),1);
    for tempi=1:numSeq(4)
        temp_boolIndex = false(1,size(sequence_length4,2));
        for tempj=1:size(sequence_length4,2)
            if sum(ismember(sequence_length4(:,tempj)',target_seqSet{4}{tempi})) == 4
                temp_boolIndex(tempj) = true;
            end
        end
        temp_index = trialIndex_length4(temp_boolIndex);
        trialIndex_length4_seq{tempi} = temp_index;
    end
    
    %% trialIndex_lengthx_location
    trialIndex_length1_location = trialIndex_length1_seq;
    
    trialIndex_length2_location = cell(numFrames,1);
    temp_length = 2;
    for tempi=1:numFrames
        tempIndex = [];
        for tempj=1:numSeq(temp_length)
            if ismember(tempi,target_seqSet{temp_length}{tempj})
                tempIndex = [tempIndex tempj];                     %#ok<*AGROW>
            end
        end
        temp_trialIndex = trialIndex_length2_seq(tempIndex);
        for tempk=1:length(tempIndex)
            trialIndex_length2_location{tempi} = [trialIndex_length2_location{tempi} temp_trialIndex{tempk}]; %#ok<*SAGROW>
        end
    end
    
    trialIndex_length3_location = cell(numFrames,1);
    temp_length = 3;
    for tempi=1:numFrames
        tempIndex = [];
        for tempj=1:numSeq(temp_length)
            if ismember(tempi,target_seqSet{temp_length}{tempj})
                tempIndex = [tempIndex tempj];                     %#ok<*AGROW>
            end
        end
        temp_trialIndex = trialIndex_length3_seq(tempIndex);
        for tempk=1:length(tempIndex)
            trialIndex_length3_location{tempi} = [trialIndex_length3_location{tempi} temp_trialIndex{tempk}]; %#ok<*SAGROW>
        end
    end
    
    trialIndex_length4_location = cell(numFrames,1);
    temp_length = 4;
    for tempi=1:numFrames
        tempIndex = [];
        for tempj=1:numSeq(temp_length)
            if ismember(tempi,target_seqSet{temp_length}{tempj})
                tempIndex = [tempIndex tempj];                     %#ok<*AGROW>
            end
        end
        temp_trialIndex = trialIndex_length4_seq(tempIndex);
        for tempk=1:length(tempIndex)
            trialIndex_length4_location{tempi} = [trialIndex_length4_location{tempi} temp_trialIndex{tempk}]; %#ok<*SAGROW>
        end
    end
    
    %% others
    choice_BoolFlag = trial_para.choiceCondition_flag==2;
    ifSelectOffloading = trial_para.ifSelectOffloading(1:trial_para.trial_count);
    ifSelectOffloading_bool = ifSelectOffloading==1;
    isCorrect_bool = trial_para.isCorrect==1;
    
    trialIndex_lengthx_bool_choice = false(4,trial_para.trial_count);
    
    trialIndex_length1_bool = false(1,trial_para.trial_count);
    trialIndex_length1_bool(trialIndex_length1) = true;
    trialIndex_lengthx_bool_choice(1,:) = choice_BoolFlag & trialIndex_length1_bool;
    
    trialIndex_length2_bool = false(1,trial_para.trial_count);
    trialIndex_length2_bool(trialIndex_length2) = true;
    trialIndex_lengthx_bool_choice(2,:) = choice_BoolFlag & trialIndex_length2_bool;
    
    trialIndex_length3_bool = false(1,trial_para.trial_count);
    trialIndex_length3_bool(trialIndex_length3) = true;
    trialIndex_lengthx_bool_choice(3,:) = choice_BoolFlag & trialIndex_length3_bool;
    
    trialIndex_length4_bool = false(1,trial_para.trial_count);
    trialIndex_length4_bool(trialIndex_length4) = true;
    trialIndex_lengthx_bool_choice(4,:) = choice_BoolFlag & trialIndex_length4_bool;
    
    trialIndex_lengthx_bool = false(4,trial_para.trial_count);
    trialIndex_lengthx_bool(1,:) = trialIndex_length1_bool;
    trialIndex_lengthx_bool(2,:) = trialIndex_length2_bool;
    trialIndex_lengthx_bool(3,:) = trialIndex_length3_bool;
    trialIndex_lengthx_bool(4,:) = trialIndex_length4_bool;
    
    trialIndex_lengthx_bool_memory = false(4,trial_para.trial_count);
    trialIndex_lengthx_bool_memory(1,:) = (~ifSelectOffloading_bool) & trialIndex_length1_bool;
    trialIndex_lengthx_bool_memory(2,:) = (~ifSelectOffloading_bool) & trialIndex_length2_bool;
    trialIndex_lengthx_bool_memory(3,:) = (~ifSelectOffloading_bool) & trialIndex_length3_bool;
    trialIndex_lengthx_bool_memory(4,:) = (~ifSelectOffloading_bool) & trialIndex_length4_bool;
    
    trialIndex_lengthx_bool_memoryCorrect = false(4,trial_para.trial_count);
    trialIndex_lengthx_bool_memoryCorrect(1,:) = (~ifSelectOffloading_bool) & trialIndex_length1_bool & isCorrect_bool;
    trialIndex_lengthx_bool_memoryCorrect(2,:) = (~ifSelectOffloading_bool) & trialIndex_length2_bool & isCorrect_bool;
    trialIndex_lengthx_bool_memoryCorrect(3,:) = (~ifSelectOffloading_bool) & trialIndex_length3_bool & isCorrect_bool;
    trialIndex_lengthx_bool_memoryCorrect(4,:) = (~ifSelectOffloading_bool) & trialIndex_length4_bool & isCorrect_bool;
    
    trialIndex_bool_memoryCorrect = (~ifSelectOffloading_bool) & isCorrect_bool;
    trialIndex_bool_memoryError = (~ifSelectOffloading_bool) & (~isCorrect_bool);
    
    trialIndex_bool_afterMemoryCorrect = false(1,trial_para.trial_count);
    trialIndex_bool_afterMemoryCorrect(2:end) = trialIndex_bool_memoryCorrect(1:end-1);
    trialIndex_bool_afterMemoryError = false(1,trial_para.trial_count);
    trialIndex_bool_afterMemoryError(2:end) = trialIndex_bool_memoryError(1:end-1);
    
    trialIndex_length1_memoryCorrect = find(trialIndex_lengthx_bool_memoryCorrect(1,:) == true);
    trialIndex_length2_memoryCorrect = find(trialIndex_lengthx_bool_memoryCorrect(2,:) == true);
    trialIndex_length3_memoryCorrect = find(trialIndex_lengthx_bool_memoryCorrect(3,:) == true);
    trialIndex_length4_memoryCorrect = find(trialIndex_lengthx_bool_memoryCorrect(4,:) == true);
    
    a = 1;
    
    validIndex_trialIndex_length1 = ismember(trialIndex_length1,trialIndex_length1_memoryCorrect);
    validIndex_trialIndex_length2 = ismember(trialIndex_length2,trialIndex_length2_memoryCorrect);
    validIndex_trialIndex_length3 = ismember(trialIndex_length3,trialIndex_length3_memoryCorrect);
    validIndex_trialIndex_length4 = ismember(trialIndex_length4,trialIndex_length4_memoryCorrect);
    
    
    
    trialIndex_lengthx_bool_choiceMemory = false(4,trial_para.trial_count);
    trialIndex_lengthx_bool_choiceMemory(1,:) = (~ifSelectOffloading_bool) & trialIndex_length1_bool & choice_BoolFlag;
    trialIndex_lengthx_bool_choiceMemory(2,:) = (~ifSelectOffloading_bool) & trialIndex_length2_bool & choice_BoolFlag;
    trialIndex_lengthx_bool_choiceMemory(3,:) = (~ifSelectOffloading_bool) & trialIndex_length3_bool & choice_BoolFlag;
    trialIndex_lengthx_bool_choiceMemory(4,:) = (~ifSelectOffloading_bool) & trialIndex_length4_bool & choice_BoolFlag;
    
    trialIndex_bool_choiceMemory = (~ifSelectOffloading_bool) & choice_BoolFlag;
    
    trialIndex_length1_choiceMemory = find(trialIndex_lengthx_bool_choiceMemory(1,:) == true);
    trialIndex_length2_choiceMemory = find(trialIndex_lengthx_bool_choiceMemory(2,:) == true);
    trialIndex_length3_choiceMemory = find(trialIndex_lengthx_bool_choiceMemory(3,:) == true);
    trialIndex_length4_choiceMemory = find(trialIndex_lengthx_bool_choiceMemory(4,:) == true);
    
    choiceMemoryIndex_trialIndex_length1 = ismember(trialIndex_length1,trialIndex_length1_choiceMemory);
    choiceMemoryIndex_trialIndex_length2 = ismember(trialIndex_length2,trialIndex_length2_choiceMemory);
    choiceMemoryIndex_trialIndex_length3 = ismember(trialIndex_length3,trialIndex_length3_choiceMemory);
    choiceMemoryIndex_trialIndex_length4 = ismember(trialIndex_length4,trialIndex_length4_choiceMemory);
    
    
    
    
    trialIndex_lengthx_bool_choiceMemoryCorrect = false(4,trial_para.trial_count);
    trialIndex_lengthx_bool_choiceMemoryCorrect(1,:) = (~ifSelectOffloading_bool) & trialIndex_length1_bool & choice_BoolFlag & isCorrect_bool;
    trialIndex_lengthx_bool_choiceMemoryCorrect(2,:) = (~ifSelectOffloading_bool) & trialIndex_length2_bool & choice_BoolFlag & isCorrect_bool;
    trialIndex_lengthx_bool_choiceMemoryCorrect(3,:) = (~ifSelectOffloading_bool) & trialIndex_length3_bool & choice_BoolFlag & isCorrect_bool;
    trialIndex_lengthx_bool_choiceMemoryCorrect(4,:) = (~ifSelectOffloading_bool) & trialIndex_length4_bool & choice_BoolFlag & isCorrect_bool;
    
    trialIndex_length1_choiceMemoryCorrect = find(trialIndex_lengthx_bool_choiceMemoryCorrect(1,:) == true);
    trialIndex_length2_choiceMemoryCorrect = find(trialIndex_lengthx_bool_choiceMemoryCorrect(2,:) == true);
    trialIndex_length3_choiceMemoryCorrect = find(trialIndex_lengthx_bool_choiceMemoryCorrect(3,:) == true);
    trialIndex_length4_choiceMemoryCorrect = find(trialIndex_lengthx_bool_choiceMemoryCorrect(4,:) == true);
    
    choiceMemoryCorrectIndex_trialIndex_length1 = ismember(trialIndex_length1,trialIndex_length1_choiceMemoryCorrect);
    choiceMemoryCorrectIndex_trialIndex_length2 = ismember(trialIndex_length2,trialIndex_length2_choiceMemoryCorrect);
    choiceMemoryCorrectIndex_trialIndex_length3 = ismember(trialIndex_length3,trialIndex_length3_choiceMemoryCorrect);
    choiceMemoryCorrectIndex_trialIndex_length4 = ismember(trialIndex_length4,trialIndex_length4_choiceMemoryCorrect);
    
    a = 1;
    
    trialIndex_lengthx_bool_noChoice = false(4,trial_para.trial_count);
    trialIndex_lengthx_bool_noChoice(1,:) = ~choice_BoolFlag & trialIndex_length1_bool;
    trialIndex_lengthx_bool_noChoice(2,:) = ~choice_BoolFlag & trialIndex_length2_bool;
    trialIndex_lengthx_bool_noChoice(3,:) = ~choice_BoolFlag & trialIndex_length3_bool;
    trialIndex_lengthx_bool_noChoice(4,:) = ~choice_BoolFlag & trialIndex_length4_bool;
    
    trialIndex_bool_noChoice = ~choice_BoolFlag;
    
    trialIndex_length1_noChoice = find(trialIndex_lengthx_bool_noChoice(1,:) == true);
    trialIndex_length2_noChoice = find(trialIndex_lengthx_bool_noChoice(2,:) == true);
    trialIndex_length3_noChoice = find(trialIndex_lengthx_bool_noChoice(3,:) == true);
    trialIndex_length4_noChoice = find(trialIndex_lengthx_bool_noChoice(4,:) == true);
    
    
    trialIndex_lengthx_bool_noChoiceCorrect = false(4,trial_para.trial_count);
    trialIndex_lengthx_bool_noChoiceCorrect(1,:) = ~choice_BoolFlag & trialIndex_length1_bool & isCorrect_bool;
    trialIndex_lengthx_bool_noChoiceCorrect(2,:) = ~choice_BoolFlag & trialIndex_length2_bool & isCorrect_bool;
    trialIndex_lengthx_bool_noChoiceCorrect(3,:) = ~choice_BoolFlag & trialIndex_length3_bool & isCorrect_bool;
    trialIndex_lengthx_bool_noChoiceCorrect(4,:) = ~choice_BoolFlag & trialIndex_length4_bool & isCorrect_bool;
    
    trialIndex_length1_noChoiceCorrect = find(trialIndex_lengthx_bool_noChoiceCorrect(1,:) == true);
    trialIndex_length2_noChoiceCorrect = find(trialIndex_lengthx_bool_noChoiceCorrect(2,:) == true);
    trialIndex_length3_noChoiceCorrect = find(trialIndex_lengthx_bool_noChoiceCorrect(3,:) == true);
    trialIndex_length4_noChoiceCorrect = find(trialIndex_lengthx_bool_noChoiceCorrect(4,:) == true);
    
    trialIndex_lengthx_bool_offload = false(4,trial_para.trial_count);
    trialIndex_lengthx_bool_offload(1,:) = ifSelectOffloading_bool & trialIndex_length1_bool;
    trialIndex_lengthx_bool_offload(2,:) = ifSelectOffloading_bool & trialIndex_length2_bool;
    trialIndex_lengthx_bool_offload(3,:) = ifSelectOffloading_bool & trialIndex_length3_bool;
    trialIndex_lengthx_bool_offload(4,:) = ifSelectOffloading_bool & trialIndex_length4_bool;
    
    trialIndex_lengthx_bool_choiceOffload = trialIndex_lengthx_bool_offload;
    
    trialIndex_bool_choiceOffload = ifSelectOffloading_bool;
    
    trialIndex_length1_choiceOffload = find(trialIndex_lengthx_bool_offload(1,:) == true);
    trialIndex_length2_choiceOffload = find(trialIndex_lengthx_bool_offload(2,:) == true);
    trialIndex_length3_choiceOffload = find(trialIndex_lengthx_bool_offload(3,:) == true);
    trialIndex_length4_choiceOffload = find(trialIndex_lengthx_bool_offload(4,:) == true);
    
    offloadIndex_trialIndex_length1 = ismember(trialIndex_length1,trialIndex_length1_choiceOffload);
    offloadIndex_trialIndex_length2 = ismember(trialIndex_length2,trialIndex_length2_choiceOffload);
    offloadIndex_trialIndex_length3 = ismember(trialIndex_length3,trialIndex_length3_choiceOffload);
    offloadIndex_trialIndex_length4 = ismember(trialIndex_length4,trialIndex_length4_choiceOffload);
    
    a = 1;
    
    %% Get F_dff_baselinePeriod
    baselinePeriod_frameIndex = zeros(sum(trialIndex_bool_choiceMemory),3);
    for tempi=1:trial_para.trial_count
        temp_trialIndex = tempi;
        temp_markers = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers;
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(1);
        baselinePeriod_frameIndex(tempi,:) = [temp_frameIndex-21 temp_frameIndex temp_frameIndex+12];
        %baselinePeriod_frameIndex(tempi,:) = [temp_frameIndex-12 temp_frameIndex temp_frameIndex+12];
    end
    baselinePeriod_interval = baselinePeriod_frameIndex(1,:)-baselinePeriod_frameIndex(1,1)+1;
    
    [baselinePeriod_frameRangeMin,~] = bounds(baselinePeriod_frameIndex,2);
    temp_frameRange = repmat(baselinePeriod_frameRangeMin,1,baselinePeriod_interval(end)-baselinePeriod_interval(1)+1);
    baselinePeriod_frameRange = temp_frameRange + ((baselinePeriod_interval(1):baselinePeriod_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,baselinePeriod_frameRange'),[size(F_dff,1),size(baselinePeriod_frameRange')]);
    F_dff_baselinePeriod = permute(temp_dff,[1 3 2]);
    
    F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,baselinePeriod_interval(1):baselinePeriod_interval(3)),3);
    
    a = 1;
    F_dff_baselineInit = mean(F_dff_baselinePeriod(:,:,1),2);
    if if_substractBaseLine == 1
        a = 1;
        
        F0_baseline = F_dff_baselineInit.*F0;
        F_dff;
        F0;
        
        F = (F_dff.*F0)+F0;
        F_dff = (F-F0-F0_baseline)./(F0+F0_baseline);
        
        %F_dff = F_dff - F0_baseline;
        
        temp_dff = reshape(F_dff(:,baselinePeriod_frameRange'),[size(F_dff,1),size(baselinePeriod_frameRange')]);
        F_dff_baselinePeriod = permute(temp_dff,[1 3 2]);
        F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,baselinePeriod_interval(1):baselinePeriod_interval(3)),3);
        
        %F_dff_baselinePeriod = F_dff_baselinePeriod - F_dff_baselineInit;
        %F_dff_baselineBin = F_dff_baselineBin - F_dff_baselineInit;
        %F_dff = F_dff - F_dff_baselineInit;
    end
    
    %% Extract length1
    %DELAY 1, DELAY 2 = 1500~1700
    % Period 1: prefixation --> TRIALID --> TARGET 1 ITEM x ON --> DELAY 1 ON --> layout disappear(true delay1)--> DELAY 1 ON + 1500 ms
    % Period 2A: DELAY 1 OFF - 500 ms--> DELAY 1 OFF --> DELAY 1 OFF + 500 ms    
    % Period 2: SELECTING AND DELAY 2 ON - 500 ms--> SELECTING AND DELAY 2 ON --> SELECTING AND DELAY 2 ON + 1500 ms
    % Period 3: SUBMIT - 600 ms --> SUBMIT --> SUBMIT + 600 ms
    %F_dff
    %F_dff_length1_period1
    %F_dff_length1_period2
    
    length1_period1_frameIndex = zeros(length(trialIndex_length1),5);%4,6
    length1_period2A_frameIndex = zeros(length(trialIndex_length1),3);
    length1_period2_frameIndex = zeros(length(trialIndex_length1),3);
    length1_period3_frameIndex = zeros(length(trialIndex_length1),3);
    for tempi=1:length(trialIndex_length1)
        temp_trialIndex = trialIndex_length1(tempi);
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(1:3);
        
        %temp_frameIndex(2) = temp_frameIndex(2) + 1;% base on the observation that T1 marker would be ahead of actual stimuli
        
        %temp_frameIndex = [temp_frameIndex(1:2) temp_frameIndex(2)+6 temp_frameIndex(end)];
        %length1_period1_frameIndex(tempi,:) = [temp_frameIndex temp_frameIndex(end)+45];
        %length1_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex temp_frameIndex(end)+45];
        %length1_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex temp_frameIndex(end)+12 temp_frameIndex(end)+45];
        length1_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex(1:2) temp_frameIndex(end)+12 temp_frameIndex(end)+45];
        
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(4);
        length1_period2A_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+15];
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(5);
        length1_period2_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+45];
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(end);
        length1_period3_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+15];%18
    end
    % length1_period1_interval = length1_period1_frameIndex(1,:)-length1_period1_frameIndex(1,1)+1;
    % length1_period2_interval = length1_period2_frameIndex(1,:)-length1_period2_frameIndex(1,1)+1;
    % length1_period3_interval = length1_period3_frameIndex(1,:)-length1_period3_frameIndex(1,1)+1;
    
    length1_period1_interval = round(median(length1_period1_frameIndex-length1_period1_frameIndex(:,1)+1,1));
    length1_period2A_interval = round(median(length1_period2A_frameIndex-length1_period2A_frameIndex(:,1)+1,1));    
    length1_period2_interval = round(median(length1_period2_frameIndex-length1_period2_frameIndex(:,1)+1,1));
    length1_period3_interval = round(median(length1_period3_frameIndex-length1_period3_frameIndex(:,1)+1,1));
    
    %length1_period1_interval = max(length1_period1_frameIndex-length1_period1_frameIndex(:,1)+1,[],1);
    %length1_period2_interval = max(length1_period2_frameIndex-length1_period2_frameIndex(:,1)+1,[],1);
    %length1_period3_interval = max(length1_period3_frameIndex-length1_period3_frameIndex(:,1)+1,[],1);
    
    length1_period1_frameIndex(:,3) = length1_period1_frameIndex(:,3) + 1;% base on the observation that T1 marker would be ahead of actual stimuli
    
    % align to fixation
    % [length1_period1_frameRangeMin,~] = bounds(length1_period1_frameIndex,2);
    % temp_frameRange = repmat(length1_period1_frameRangeMin,1,delength1_period1_interval(end)-length1_period1_interval(1)+1);
    % length1_period1_frameRange = temp_frameRange + ((length1_period1_interval(1):length1_period1_interval(end))-1);
    
    % align to T1
    length1_period1_frameRangeMin = length1_period1_frameIndex(:,3);
    temp_frameRange = repmat(length1_period1_frameRangeMin,1,length1_period1_interval(end)-length1_period1_interval(1)+1);
    length1_period1_frameRange_alignT1 = temp_frameRange + ((length1_period1_interval(1):length1_period1_interval(end))-length1_period1_interval(3));
    
    a = 1;
    % align to Delay
    length1_period1_frameRangeMin = length1_period1_frameIndex(:,4);
    temp_frameRange = repmat(length1_period1_frameRangeMin,1,length1_period1_interval(end)-length1_period1_interval(1)+1);
    length1_period1_frameRange_alignDelay = temp_frameRange + ((length1_period1_interval(1):length1_period1_interval(end))-length1_period1_interval(4));
    
    temp1 = length1_period1_frameRange_alignT1(:,1:(length1_period1_interval(4)-1));
    temp2 = length1_period1_frameRange_alignDelay(:,length1_period1_interval(4):length1_period1_interval(end));
    length1_period1_frameRange = [temp1 temp2];
    
    
    % b1 = reshape(F_dff(2,length1_period1_frameRange(3,:)'),size(length1_period1_frameRange(3,:)'));
    temp_dff = reshape(F_dff(:,length1_period1_frameRange'),[size(F_dff,1),size(length1_period1_frameRange')]);
    F_dff_length1_period1 = permute(temp_dff,[1 3 2]);
    
    a = 1;
    
    [length1_period2A_frameRangeMin,~] = bounds(length1_period2A_frameIndex,2);
    temp_frameRange = repmat(length1_period2A_frameRangeMin,1,length1_period2A_interval(end)-length1_period2A_interval(1)+1);
    length1_period2A_frameRange = temp_frameRange + ((length1_period2A_interval(1):length1_period2A_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length1_period2A_frameRange'),[size(F_dff,1),size(length1_period2A_frameRange')]);
    F_dff_length1_period2A = permute(temp_dff,[1 3 2]);
    
    [length1_period2_frameRangeMin,~] = bounds(length1_period2_frameIndex,2);
    temp_frameRange = repmat(length1_period2_frameRangeMin,1,length1_period2_interval(end)-length1_period2_interval(1)+1);
    length1_period2_frameRange = temp_frameRange + ((length1_period2_interval(1):length1_period2_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length1_period2_frameRange'),[size(F_dff,1),size(length1_period2_frameRange')]);
    F_dff_length1_period2 = permute(temp_dff,[1 3 2]);
    
    [length1_period3_frameRangeMin,~] = bounds(length1_period3_frameIndex,2);
    temp_frameRange = repmat(length1_period3_frameRangeMin,1,length1_period3_interval(end)-length1_period3_interval(1)+1);
    length1_period3_frameRange = temp_frameRange + ((length1_period3_interval(1):length1_period3_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length1_period3_frameRange'),[size(F_dff,1),size(length1_period3_frameRange')]);
    F_dff_length1_period3 = permute(temp_dff,[1 3 2]);
    
    a = 1;
    %% Extract length2
    %DELAY 1, DELAY 2 = 1500~1700
    % Period 1: prefixation --> TRIALID --> TARGET 1 ITEM x ON --> TARGET 2 ITEM x ON --> DELAY 1 ON --> layout disappear--> DELAY 1 ON + 1500 ms
    % Period 2A: DELAY 1 OFF - 500 ms--> DELAY 1 OFF --> DELAY 1 OFF + 500 ms        
    % Period 2: SELECTING AND DELAY 2 ON - 500 ms--> SELECTING AND DELAY 2 ON --> SELECTING AND DELAY 2 ON + 1500 ms
    % Period 3: SUBMIT - 600 ms --> SUBMIT --> SUBMIT + 600 ms
    %F_dff
    %F_dff_length2_period1
    %F_dff_length2_period2
    length2_period1_frameIndex = zeros(length(trialIndex_length2),6);%5,7
    length2_period2A_frameIndex = zeros(length(trialIndex_length2),3);    
    length2_period2_frameIndex = zeros(length(trialIndex_length2),3);
    length2_period3_frameIndex = zeros(length(trialIndex_length2),3);
    for tempi=1:length(trialIndex_length2)
        temp_trialIndex = trialIndex_length2(tempi);
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(1:4);
        
        %temp_frameIndex(2) = temp_frameIndex(2) + 1;% base on the observation that T1 marker would be ahead of actual stimuli
        
        %length2_period1_frameIndex(tempi,:) = [temp_frameIndex temp_frameIndex(end)+45];
        %length2_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex temp_frameIndex(end)+45];
        %length2_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex temp_frameIndex(end)+12 temp_frameIndex(end)+45];
        length2_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex(1:3) temp_frameIndex(end)+12 temp_frameIndex(end)+45];        
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(5);
        length2_period2A_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+15];
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(6);
        length2_period2_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+45];
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(end);
        length2_period3_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+15];%18
    end
    % length2_period1_interval = length2_period1_frameIndex(1,:)-length2_period1_frameIndex(1,1)+1;
    % length2_period2_interval = length2_period2_frameIndex(1,:)-length2_period2_frameIndex(1,1)+1;
    % length2_period3_interval = length2_period3_frameIndex(1,:)-length2_period3_frameIndex(1,1)+1;
    length2_period1_interval = round(median(length2_period1_frameIndex-length2_period1_frameIndex(:,1)+1,1));
    length2_period2A_interval = round(median(length2_period2A_frameIndex-length2_period2A_frameIndex(:,1)+1,1));    
    length2_period2_interval = round(median(length2_period2_frameIndex-length2_period2_frameIndex(:,1)+1,1));
    length2_period3_interval = round(median(length2_period3_frameIndex-length2_period3_frameIndex(:,1)+1,1));
    %length2_period1_interval = max(length2_period1_frameIndex-length2_period1_frameIndex(:,1)+1,[],1);
    %length2_period2_interval = max(length2_period2_frameIndex-length2_period2_frameIndex(:,1)+1,[],1);
    %length2_period3_interval = max(length2_period3_frameIndex-length2_period3_frameIndex(:,1)+1,[],1);
    
    length2_period1_frameIndex(:,3) = length2_period1_frameIndex(:,3) + 1;% base on the observation that T1 marker would be ahead of actual stimuli
    
    % align to fixation
    % [length2_period1_frameRangeMin,~] = bounds(length2_period1_frameIndex,2);
    % temp_frameRange = repmat(length2_period1_frameRangeMin,1,length2_period1_interval(end)-length2_period1_interval(1)+1);
    % length2_period1_frameRange = temp_frameRange + ((length2_period1_interval(1):length2_period1_interval(end))-1);
    
    % align to T1
    length2_period1_frameRangeMin = length2_period1_frameIndex(:,3);
    temp_frameRange = repmat(length2_period1_frameRangeMin,1,length2_period1_interval(end)-length2_period1_interval(1)+1);
    length2_period1_frameRange_alignT1 = temp_frameRange + ((length2_period1_interval(1):length2_period1_interval(end))-length2_period1_interval(3));
    
    % align to T2
    length2_period1_frameRangeMin = length2_period1_frameIndex(:,4);
    temp_frameRange = repmat(length2_period1_frameRangeMin,1,length2_period1_interval(end)-length2_period1_interval(1)+1);
    length2_period1_frameRange_alignT2 = temp_frameRange + ((length2_period1_interval(1):length2_period1_interval(end))-length2_period1_interval(4));
    
    % align to Delay
    length2_period1_frameRangeMin = length2_period1_frameIndex(:,5);
    temp_frameRange = repmat(length2_period1_frameRangeMin,1,length2_period1_interval(end)-length2_period1_interval(1)+1);
    length2_period1_frameRange_alignDelay = temp_frameRange + ((length2_period1_interval(1):length2_period1_interval(end))-length2_period1_interval(5));
    
    temp1 = length2_period1_frameRange_alignT1(:,1:(length2_period1_interval(4)-1));
    temp2 = length2_period1_frameRange_alignT2(:,length2_period1_interval(4):(length2_period1_interval(5)-1));
    temp3 = length2_period1_frameRange_alignDelay(:,length2_period1_interval(5):length2_period1_interval(end));
    length2_period1_frameRange = [temp1 temp2 temp3];
    
    
    temp_dff = reshape(F_dff(:,length2_period1_frameRange'),[size(F_dff,1),size(length2_period1_frameRange')]);
    F_dff_length2_period1 = permute(temp_dff,[1 3 2]);
    
    
    [length2_period2A_frameRangeMin,~] = bounds(length2_period2A_frameIndex,2);
    temp_frameRange = repmat(length2_period2A_frameRangeMin,1,length2_period2A_interval(end)-length2_period2A_interval(1)+1);
    length2_period2A_frameRange = temp_frameRange + ((length2_period2A_interval(1):length2_period2A_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length2_period2A_frameRange'),[size(F_dff,1),size(length2_period2A_frameRange')]);
    F_dff_length2_period2A = permute(temp_dff,[1 3 2]);

    
    [length2_period2_frameRangeMin,~] = bounds(length2_period2_frameIndex,2);
    temp_frameRange = repmat(length2_period2_frameRangeMin,1,length2_period2_interval(end)-length2_period2_interval(1)+1);
    length2_period2_frameRange = temp_frameRange + ((length2_period2_interval(1):length2_period2_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length2_period2_frameRange'),[size(F_dff,1),size(length2_period2_frameRange')]);
    F_dff_length2_period2 = permute(temp_dff,[1 3 2]);
    
    [length2_period3_frameRangeMin,~] = bounds(length2_period3_frameIndex,2);
    temp_frameRange = repmat(length2_period3_frameRangeMin,1,length2_period3_interval(end)-length2_period3_interval(1)+1);
    length2_period3_frameRange = temp_frameRange + ((length2_period3_interval(1):length2_period3_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length2_period3_frameRange'),[size(F_dff,1),size(length2_period3_frameRange')]);
    F_dff_length2_period3 = permute(temp_dff,[1 3 2]);
    
    %% Extract length3
    %DELAY 1, DELAY 2 = 1500~1700
    % Period 1: prefixation --> TRIALID --> TARGET 1 ITEM x ON --> TARGET 2 ITEM x ON --> TARGET 3 ITEM x ON --> DELAY 1 ON --> layout disappear--> DELAY 1 ON + 1500 ms
    % Period 2A: DELAY 1 OFF - 500 ms--> DELAY 1 OFF --> DELAY 1 OFF + 500 ms            
    % Period 2: SELECTING AND DELAY 2 ON - 500 ms--> SELECTING AND DELAY 2 ON --> SELECTING AND DELAY 2 ON + 1500 ms
    % Period 3: SUBMIT - 600 ms --> SUBMIT --> SUBMIT + 600 ms
    
    %F_dff
    %F_dff_length3_period1
    %F_dff_length3_period2
    length3_period1_frameIndex = zeros(length(trialIndex_length3),7);%6,8
    length3_period2A_frameIndex = zeros(length(trialIndex_length3),3);    
    length3_period2_frameIndex = zeros(length(trialIndex_length3),3);
    length3_period3_frameIndex = zeros(length(trialIndex_length3),3);
    for tempi=1:length(trialIndex_length3)
        temp_trialIndex = trialIndex_length3(tempi);
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(1:5);
        
        %temp_frameIndex(2) = temp_frameIndex(2) + 1;% base on the observation that T1 marker would be ahead of actual stimuli
        
        %length3_period1_frameIndex(tempi,:) = [temp_frameIndex temp_frameIndex(end)+45];
        %length3_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex temp_frameIndex(end)+12 temp_frameIndex(end)+45];
        length3_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex(1:4) temp_frameIndex(end)+12 temp_frameIndex(end)+45];
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(6);
        length3_period2A_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+15];
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(7);
        length3_period2_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+45];
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(end);
        length3_period3_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+15];%18
    end
    % length3_period1_interval = length3_period1_frameIndex(1,:)-length3_period1_frameIndex(1,1)+1;
    % length3_period2_interval = length3_period2_frameIndex(1,:)-length3_period2_frameIndex(1,1)+1;
    % length3_period3_interval = length3_period3_frameIndex(1,:)-length3_period3_frameIndex(1,1)+1;
    length3_period1_interval = round(median(length3_period1_frameIndex-length3_period1_frameIndex(:,1)+1,1));
    length3_period2A_interval = round(median(length3_period2A_frameIndex-length3_period2A_frameIndex(:,1)+1,1));    
    length3_period2_interval = round(median(length3_period2_frameIndex-length3_period2_frameIndex(:,1)+1,1));
    length3_period3_interval = round(median(length3_period3_frameIndex-length3_period3_frameIndex(:,1)+1,1));
    %length3_period1_interval = max(length3_period1_frameIndex-length3_period1_frameIndex(:,1)+1,[],1);
    %length3_period2_interval = max(length3_period2_frameIndex-length3_period2_frameIndex(:,1)+1,[],1);
    %length3_period3_interval = max(length3_period3_frameIndex-length3_period3_frameIndex(:,1)+1,[],1);
    
    length3_period1_frameIndex(:,3) = length3_period1_frameIndex(:,3) + 1;% base on the observation that T1 marker would be ahead of actual stimuli
    
    % % align to fixation
    % [length3_period1_frameRangeMin,~] = bounds(length3_period1_frameIndex,2);
    % temp_frameRange = repmat(length3_period1_frameRangeMin,1,length3_period1_interval(end)-length3_period1_interval(1)+1);
    % length3_period1_frameRange = temp_frameRange + ((length3_period1_interval(1):length3_period1_interval(end))-1);
    
    % align to T1
    length3_period1_frameRangeMin = length3_period1_frameIndex(:,3);
    temp_frameRange = repmat(length3_period1_frameRangeMin,1,length3_period1_interval(end)-length3_period1_interval(1)+1);
    length3_period1_frameRange_alignT1 = temp_frameRange + ((length3_period1_interval(1):length3_period1_interval(end))-length3_period1_interval(3));
    
    % align to T2
    length3_period1_frameRangeMin = length3_period1_frameIndex(:,4);
    temp_frameRange = repmat(length3_period1_frameRangeMin,1,length3_period1_interval(end)-length3_period1_interval(1)+1);
    length3_period1_frameRange_alignT2 = temp_frameRange + ((length3_period1_interval(1):length3_period1_interval(end))-length3_period1_interval(4));
    
    % align to T3
    length3_period1_frameRangeMin = length3_period1_frameIndex(:,5);
    temp_frameRange = repmat(length3_period1_frameRangeMin,1,length3_period1_interval(end)-length3_period1_interval(1)+1);
    length3_period1_frameRange_alignT3 = temp_frameRange + ((length3_period1_interval(1):length3_period1_interval(end))-length3_period1_interval(5));
    
    % align to Delay
    length3_period1_frameRangeMin = length3_period1_frameIndex(:,6);
    temp_frameRange = repmat(length3_period1_frameRangeMin,1,length3_period1_interval(end)-length3_period1_interval(1)+1);
    length3_period1_frameRange_alignDelay = temp_frameRange + ((length3_period1_interval(1):length3_period1_interval(end))-length3_period1_interval(6));
    
    temp1 = length3_period1_frameRange_alignT1(:,1:(length3_period1_interval(4)-1));
    temp2 = length3_period1_frameRange_alignT2(:,length3_period1_interval(4):(length3_period1_interval(5)-1));
    temp3 = length3_period1_frameRange_alignT3(:,length3_period1_interval(5):(length3_period1_interval(6)-1));
    temp4 = length3_period1_frameRange_alignDelay(:,length3_period1_interval(6):length3_period1_interval(end));
    length3_period1_frameRange = [temp1 temp2 temp3 temp4];
    
    
    temp_dff = reshape(F_dff(:,length3_period1_frameRange'),[size(F_dff,1),size(length3_period1_frameRange')]);
    F_dff_length3_period1 = permute(temp_dff,[1 3 2]);
    
    [length3_period2_frameRangeMin,~] = bounds(length3_period2_frameIndex,2);
    temp_frameRange = repmat(length3_period2_frameRangeMin,1,length3_period2_interval(end)-length3_period2_interval(1)+1);
    length3_period2_frameRange = temp_frameRange + ((length3_period2_interval(1):length3_period2_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length3_period2_frameRange'),[size(F_dff,1),size(length3_period2_frameRange')]);
    F_dff_length3_period2 = permute(temp_dff,[1 3 2]);
    
    
    [length3_period2A_frameRangeMin,~] = bounds(length3_period2A_frameIndex,2);
    temp_frameRange = repmat(length3_period2A_frameRangeMin,1,length3_period2A_interval(end)-length3_period2A_interval(1)+1);
    length3_period2A_frameRange = temp_frameRange + ((length3_period2A_interval(1):length3_period2A_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length3_period2A_frameRange'),[size(F_dff,1),size(length3_period2A_frameRange')]);
    F_dff_length3_period2A = permute(temp_dff,[1 3 2]);

    
    [length3_period3_frameRangeMin,~] = bounds(length3_period3_frameIndex,2);
    temp_frameRange = repmat(length3_period3_frameRangeMin,1,length3_period3_interval(end)-length3_period3_interval(1)+1);
    length3_period3_frameRange = temp_frameRange + ((length3_period3_interval(1):length3_period3_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length3_period3_frameRange'),[size(F_dff,1),size(length3_period3_frameRange')]);
    F_dff_length3_period3 = permute(temp_dff,[1 3 2]);
    
    
    %% Extract length4
    %DELAY 1, DELAY 2 = 1500~1700
    % Period 1: prefixation --> TRIALID --> TARGET 1 ITEM x ON --> TARGET 2 ITEM x ON --> TARGET 3 ITEM x ON --> TARGET 4 ITEM x ON --> DELAY 1 ON --> layout disappear--> DELAY 1 ON + 1500 ms
    % Period 2A: DELAY 1 OFF - 500 ms--> DELAY 1 OFF --> DELAY 1 OFF + 500 ms            
    % Period 2: SELECTING AND DELAY 2 ON - 500 ms--> SELECTING AND DELAY 2 ON --> SELECTING AND DELAY 2 ON + 1500 ms
    % Period 3: SUBMIT - 600 ms --> SUBMIT --> SUBMIT + 600 ms
    
    %F_dff
    %F_dff_length4_period1
    %F_dff_length4_period2
    length4_period1_frameIndex = zeros(length(trialIndex_length4),8);
    length4_period2A_frameIndex = zeros(length(trialIndex_length4),3);    
    length4_period2_frameIndex = zeros(length(trialIndex_length4),3);
    length4_period3_frameIndex = zeros(length(trialIndex_length4),3);
    for tempi=1:length(trialIndex_length4)
        temp_trialIndex = trialIndex_length4(tempi);
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(1:6);
        
        length4_period1_frameIndex(tempi,:) = [temp_frameIndex(1)-21 temp_frameIndex(1:5) temp_frameIndex(end)+12 temp_frameIndex(end)+45];
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(7);
        length4_period2A_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+15];
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(8);
        length4_period2_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+45];
        
        temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(end);
        length4_period3_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+15];%18
    end
    length4_period1_interval = round(median(length4_period1_frameIndex-length4_period1_frameIndex(:,1)+1,1));
    length4_period2A_interval = round(median(length4_period2A_frameIndex-length4_period2A_frameIndex(:,1)+1,1));    
    length4_period2_interval = round(median(length4_period2_frameIndex-length4_period2_frameIndex(:,1)+1,1));
    length4_period3_interval = round(median(length4_period3_frameIndex-length4_period3_frameIndex(:,1)+1,1));
    
    length4_period1_frameIndex(:,3) = length4_period1_frameIndex(:,3) + 1;% base on the observation that T1 marker would be ahead of actual stimuli
    
    % % align to fixation
    % [length4_period1_frameRangeMin,~] = bounds(length4_period1_frameIndex,2);
    % temp_frameRange = repmat(length4_period1_frameRangeMin,1,length4_period1_interval(end)-length4_period1_interval(1)+1);
    % length4_period1_frameRange = temp_frameRange + ((length4_period1_interval(1):length4_period1_interval(end))-1);
    
    % align to T1
    length4_period1_frameRangeMin = length4_period1_frameIndex(:,3);
    temp_frameRange = repmat(length4_period1_frameRangeMin,1,length4_period1_interval(end)-length4_period1_interval(1)+1);
    length4_period1_frameRange_alignT1 = temp_frameRange + ((length4_period1_interval(1):length4_period1_interval(end))-length4_period1_interval(3));
    
    % align to T2
    length4_period1_frameRangeMin = length4_period1_frameIndex(:,4);
    temp_frameRange = repmat(length4_period1_frameRangeMin,1,length4_period1_interval(end)-length4_period1_interval(1)+1);
    length4_period1_frameRange_alignT2 = temp_frameRange + ((length4_period1_interval(1):length4_period1_interval(end))-length4_period1_interval(4));
    
    % align to T3
    length4_period1_frameRangeMin = length4_period1_frameIndex(:,5);
    temp_frameRange = repmat(length4_period1_frameRangeMin,1,length4_period1_interval(end)-length4_period1_interval(1)+1);
    length4_period1_frameRange_alignT3 = temp_frameRange + ((length4_period1_interval(1):length4_period1_interval(end))-length4_period1_interval(5));
    
    % align to T4
    length4_period1_frameRangeMin = length4_period1_frameIndex(:,6);
    temp_frameRange = repmat(length4_period1_frameRangeMin,1,length4_period1_interval(end)-length4_period1_interval(1)+1);
    length4_period1_frameRange_alignT4 = temp_frameRange + ((length4_period1_interval(1):length4_period1_interval(end))-length4_period1_interval(6));
    
    % align to Delay
    length4_period1_frameRangeMin = length4_period1_frameIndex(:,7);
    temp_frameRange = repmat(length4_period1_frameRangeMin,1,length4_period1_interval(end)-length4_period1_interval(1)+1);
    length4_period1_frameRange_alignDelay = temp_frameRange + ((length4_period1_interval(1):length4_period1_interval(end))-length4_period1_interval(7));
    
    temp1 = length4_period1_frameRange_alignT1(:,1:(length4_period1_interval(4)-1));
    temp2 = length4_period1_frameRange_alignT2(:,length4_period1_interval(4):(length4_period1_interval(5)-1));
    temp3 = length4_period1_frameRange_alignT3(:,length4_period1_interval(5):(length4_period1_interval(6)-1));
    temp4 = length4_period1_frameRange_alignT4(:,length4_period1_interval(6):(length4_period1_interval(7)-1));    
    temp5 = length4_period1_frameRange_alignDelay(:,length4_period1_interval(7):length4_period1_interval(end));
    length4_period1_frameRange = [temp1 temp2 temp3 temp4 temp5];
    
    
    temp_dff = reshape(F_dff(:,length4_period1_frameRange'),[size(F_dff,1),size(length4_period1_frameRange')]);
    F_dff_length4_period1 = permute(temp_dff,[1 3 2]);
    
    [length4_period2_frameRangeMin,~] = bounds(length4_period2_frameIndex,2);
    temp_frameRange = repmat(length4_period2_frameRangeMin,1,length4_period2_interval(end)-length4_period2_interval(1)+1);
    length4_period2_frameRange = temp_frameRange + ((length4_period2_interval(1):length4_period2_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length4_period2_frameRange'),[size(F_dff,1),size(length4_period2_frameRange')]);
    F_dff_length4_period2 = permute(temp_dff,[1 3 2]);
    
    
    [length4_period2A_frameRangeMin,~] = bounds(length4_period2A_frameIndex,2);
    temp_frameRange = repmat(length4_period2A_frameRangeMin,1,length4_period2A_interval(end)-length4_period2A_interval(1)+1);
    length4_period2A_frameRange = temp_frameRange + ((length4_period2A_interval(1):length4_period2A_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length4_period2A_frameRange'),[size(F_dff,1),size(length4_period2A_frameRange')]);
    F_dff_length4_period2A = permute(temp_dff,[1 3 2]);

    
    [length4_period3_frameRangeMin,~] = bounds(length4_period3_frameIndex,2);
    temp_frameRange = repmat(length4_period3_frameRangeMin,1,length4_period3_interval(end)-length4_period3_interval(1)+1);
    length4_period3_frameRange = temp_frameRange + ((length4_period3_interval(1):length4_period3_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,length4_period3_frameRange'),[size(F_dff,1),size(length4_period3_frameRange')]);
    F_dff_length4_period3 = permute(temp_dff,[1 3 2]);    
    
    
    decodingData = struct;
    
    decodingData.F_dff_length1_period1 = F_dff_length1_period1;
    decodingData.F_dff_length1_period2A = F_dff_length1_period2A;    
    decodingData.F_dff_length1_period2 = F_dff_length1_period2;
    decodingData.F_dff_length2_period1 = F_dff_length2_period1;
    decodingData.F_dff_length2_period2A = F_dff_length2_period2A;    
    decodingData.F_dff_length2_period2 = F_dff_length2_period2;
    decodingData.F_dff_length3_period1 = F_dff_length3_period1;
    decodingData.F_dff_length3_period2A = F_dff_length3_period2A;    
    decodingData.F_dff_length3_period2 = F_dff_length3_period2;
    decodingData.F_dff_length4_period1 = F_dff_length4_period1;
    decodingData.F_dff_length4_period2A = F_dff_length4_period2A;    
    decodingData.F_dff_length4_period2 = F_dff_length4_period2;    
    
    decodingData.sequence_length1 = sequence_length1;
    decodingData.sequence_length2 = sequence_length2;
    decodingData.sequence_length3 = sequence_length3;
    decodingData.sequence_length4 = sequence_length4;
    
    decodingData.length1_period1_interval = length1_period1_interval;
    decodingData.length2_period1_interval = length2_period1_interval;
    decodingData.length3_period1_interval = length3_period1_interval;
    decodingData.length4_period1_interval = length4_period1_interval;    
    decodingData.lengthx_period2_interval = length1_period2_interval;
    
    decodingData.trialIndex_length1 = trialIndex_length1;
    decodingData.trialIndex_length2 = trialIndex_length2;
    decodingData.trialIndex_length3 = trialIndex_length3;
    decodingData.trialIndex_length4 = trialIndex_length4;
    
    
    a = 1;
    
    % F_dff_length1_period2;
    % F_dff_length2_period1;
    % F_dff_length2_period2;
    % F_dff_length3_period1;
    % F_dff_length3_period2;
    %
    % sequence_length1;
    % sequence_length2;
    % sequence_length3;
    %
    % length1_period1_interval;
    % length2_period1_interval;
    % length3_period1_interval;
    % lengthx_period2_interval = length1_period2_interval;
    %
    % trialIndex_length1;
    % trialIndex_length2;
    % trialIndex_length3;
    
    a = 1;
    
    %% Extract decision
    trial_para;
    trialIndex_length1;
    
    
    %trialIndex_lengthx_bool_choiceMemory = false(3,trial_para.trial_count);
    %trialIndex_lengthx_bool_choiceOffload = false(3,trial_para.trial_count);
    %for tempi=1:3
    %    trialIndex_lengthx_bool_choiceMemory(tempi,:) = trialIndex_lengthx_bool_choice(tempi,:) & (ifSelectOffloading==-1);
    %    trialIndex_lengthx_bool_choiceOffload(tempi,:) = trialIndex_lengthx_bool_choice(tempi,:) & (ifSelectOffloading==1);
    %end
    
    % trialIndex_bool_choiceMemory = sum(trialIndex_lengthx_bool_choiceMemory,1)>0;
    % trialIndex_bool_choiceOffload = sum(trialIndex_lengthx_bool_choiceOffload,1)>0;
    %
    % trialIndex_lengthx_bool_noChoice;
    % trialIndex_bool_noChoice = sum(trialIndex_lengthx_bool_noChoice,1)>0;
    a = 1;
    
    
    % Decision period : SELECTING AND DELAY 2 ON - 1500 ms--> SELECTING AND DELAY 2 ON --> SELECTING AND DELAY 2 ON + 1500 ms
    %F_dff
    %F_dff_decisionPeriod_choiceMemory
    %F_dff_decisionPeriod_choiceOffload
    
    F_dff;
    
    % Get F_dff_decisionPeriod_choiceMemory
    decisionPeriodA_choiceMemory_frameIndex = zeros(sum(trialIndex_bool_choiceMemory),3);
    decisionPeriod_choiceMemory_frameIndex = zeros(sum(trialIndex_bool_choiceMemory),3);
    for tempi=1:sum(trialIndex_bool_choiceMemory)
        temp_trialIndex = find(trialIndex_bool_choiceMemory == 1,tempi);
        temp_trialIndex = temp_trialIndex(end);
                       
        temp_markers = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers;
        
        for tempj=1:length(temp_markers)
            if strcmp(temp_markers{tempj}, 'DELAY 1 OFF') == 1
                temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(tempj);
                break
            end
        end
        decisionPeriodA_choiceMemory_frameIndex(tempi,:) = [temp_frameIndex-33 temp_frameIndex temp_frameIndex+15];        
        
        for tempj=1:length(temp_markers)
            if strcmp(temp_markers{tempj}, 'SELECTING AND DELAY 2 ON') == 1
                temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(tempj);
                break
            end
        end
        decisionPeriod_choiceMemory_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+45];

    end
    decisionPeriodA_interval = decisionPeriodA_choiceMemory_frameIndex(1,:)-decisionPeriodA_choiceMemory_frameIndex(1,1)+1;    
    decisionPeriod_interval = decisionPeriod_choiceMemory_frameIndex(1,:)-decisionPeriod_choiceMemory_frameIndex(1,1)+1;
    
    [decisionPeriod_choiceMemory_frameRangeMin,~] = bounds(decisionPeriod_choiceMemory_frameIndex,2);
    temp_frameRange = repmat(decisionPeriod_choiceMemory_frameRangeMin,1,decisionPeriod_interval(end)-decisionPeriod_interval(1)+1);
    decisionPeriod_choiceMemory_frameRange = temp_frameRange + ((decisionPeriod_interval(1):decisionPeriod_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,decisionPeriod_choiceMemory_frameRange'),[size(F_dff,1),size(decisionPeriod_choiceMemory_frameRange')]);
    F_dff_decisionPeriod_choiceMemory = permute(temp_dff,[1 3 2]);
    
    
    [decisionPeriodA_choiceMemory_frameRangeMin,~] = bounds(decisionPeriodA_choiceMemory_frameIndex,2);
    temp_frameRange = repmat(decisionPeriodA_choiceMemory_frameRangeMin,1,decisionPeriodA_interval(end)-decisionPeriodA_interval(1)+1);
    decisionPeriodA_choiceMemory_frameRange = temp_frameRange + ((decisionPeriodA_interval(1):decisionPeriodA_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,decisionPeriodA_choiceMemory_frameRange'),[size(F_dff,1),size(decisionPeriodA_choiceMemory_frameRange')]);
    F_dff_decisionPeriodA_choiceMemory = permute(temp_dff,[1 3 2]);
    
    a = 1;
    
    % Get F_dff_decisionPeriod_choiceOffload
    decisionPeriodA_choiceOffload_frameIndex = zeros(sum(trialIndex_bool_choiceOffload),3);
    decisionPeriod_choiceOffload_frameIndex = zeros(sum(trialIndex_bool_choiceOffload),3);    
    for tempi=1:sum(trialIndex_bool_choiceOffload)
        temp_trialIndex = find(trialIndex_bool_choiceOffload == 1,tempi);
        temp_trialIndex = temp_trialIndex(end);
        
        temp_markers = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers;
        
        for tempj=1:length(temp_markers)
            if strcmp(temp_markers{tempj}, 'DELAY 1 OFF') == 1
                temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(tempj);
                break
            end
        end
        decisionPeriodA_choiceOffload_frameIndex(tempi,:) = [temp_frameIndex-33 temp_frameIndex temp_frameIndex+15];        
        
        for tempj=1:length(temp_markers)
            if strcmp(temp_markers{tempj}, 'SELECTING AND DELAY 2 ON') == 1
                temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(tempj);
                break
            end
        end
        decisionPeriod_choiceOffload_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+45];
    end
    decisionPeriodA_interval = decisionPeriodA_choiceOffload_frameIndex(1,:)-decisionPeriodA_choiceOffload_frameIndex(1,1)+1;    
    decisionPeriod_interval = decisionPeriod_choiceOffload_frameIndex(1,:)-decisionPeriod_choiceOffload_frameIndex(1,1)+1;
    
    [decisionPeriod_choiceOffload_frameRangeMin,~] = bounds(decisionPeriod_choiceOffload_frameIndex,2);
    temp_frameRange = repmat(decisionPeriod_choiceOffload_frameRangeMin,1,decisionPeriod_interval(end)-decisionPeriod_interval(1)+1);
    decisionPeriod_choiceOffload_frameRange = temp_frameRange + ((decisionPeriod_interval(1):decisionPeriod_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,decisionPeriod_choiceOffload_frameRange'),[size(F_dff,1),size(decisionPeriod_choiceOffload_frameRange')]);
    F_dff_decisionPeriod_choiceOffload = permute(temp_dff,[1 3 2]);
    
    
    [decisionPeriodA_choiceOffload_frameRangeMin,~] = bounds(decisionPeriodA_choiceOffload_frameIndex,2);
    temp_frameRange = repmat(decisionPeriodA_choiceOffload_frameRangeMin,1,decisionPeriodA_interval(end)-decisionPeriodA_interval(1)+1);
    decisionPeriodA_choiceOffload_frameRange = temp_frameRange + ((decisionPeriodA_interval(1):decisionPeriodA_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,decisionPeriodA_choiceOffload_frameRange'),[size(F_dff,1),size(decisionPeriodA_choiceOffload_frameRange')]);
    F_dff_decisionPeriodA_choiceOffload = permute(temp_dff,[1 3 2]);
    
    
    
    % Get F_dff_decisionPeriod_noChoice
    decisionPeriodA_noChoice_frameIndex = zeros(sum(trialIndex_bool_noChoice),3);    
    decisionPeriod_noChoice_frameIndex = zeros(sum(trialIndex_bool_noChoice),3);
    for tempi=1:sum(trialIndex_bool_noChoice)
        temp_trialIndex = find(trialIndex_bool_noChoice == 1,tempi);
        temp_trialIndex = temp_trialIndex(end);
        
        temp_markers = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers;
        
        for tempj=1:length(temp_markers)
            if strcmp(temp_markers{tempj}, 'DELAY 1 OFF') == 1
                temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(tempj);
                break
            end
        end
        decisionPeriodA_noChoice_frameIndex(tempi,:) = [temp_frameIndex-33 temp_frameIndex temp_frameIndex+15];        
        
        for tempj=1:length(temp_markers)
            if strcmp(temp_markers{tempj}, 'SELECTING AND DELAY 2 ON') == 1
                temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(tempj);
                break
            end
        end
        decisionPeriod_noChoice_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+45];
    end
    decisionPeriodA_interval = decisionPeriodA_noChoice_frameIndex(1,:)-decisionPeriodA_noChoice_frameIndex(1,1)+1;    
    decisionPeriod_interval = decisionPeriod_noChoice_frameIndex(1,:)-decisionPeriod_noChoice_frameIndex(1,1)+1;
    
    [decisionPeriod_noChoice_frameRangeMin,~] = bounds(decisionPeriod_noChoice_frameIndex,2);
    temp_frameRange = repmat(decisionPeriod_noChoice_frameRangeMin,1,decisionPeriod_interval(end)-decisionPeriod_interval(1)+1);
    decisionPeriod_noChoice_frameRange = temp_frameRange + ((decisionPeriod_interval(1):decisionPeriod_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,decisionPeriod_noChoice_frameRange'),[size(F_dff,1),size(decisionPeriod_noChoice_frameRange')]);
    F_dff_decisionPeriod_noChoice = permute(temp_dff,[1 3 2]);
    
    [decisionPeriodA_noChoice_frameRangeMin,~] = bounds(decisionPeriodA_noChoice_frameIndex,2);
    temp_frameRange = repmat(decisionPeriodA_noChoice_frameRangeMin,1,decisionPeriodA_interval(end)-decisionPeriodA_interval(1)+1);
    decisionPeriodA_noChoice_frameRange = temp_frameRange + ((decisionPeriodA_interval(1):decisionPeriodA_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,decisionPeriodA_noChoice_frameRange'),[size(F_dff,1),size(decisionPeriodA_noChoice_frameRange')]);
    F_dff_decisionPeriodA_noChoice = permute(temp_dff,[1 3 2]);
    
    a = 1;
    F_dff_decisionPeriod_choiceMemory;
    F_dff_decisionPeriod_choiceOffload;
    F_dff_decisionPeriod_noChoice;
    
    %% do some filter in F_dff
    
    
    F_dff_length1_period1_raw = F_dff_length1_period1;
    F_dff_length1_period2A_raw = F_dff_length1_period2A;
    F_dff_length1_period2_raw = F_dff_length1_period2;
    F_dff_length1_period3_raw = F_dff_length1_period3;    
    F_dff_length2_period1_raw = F_dff_length2_period1;
    F_dff_length2_period2A_raw = F_dff_length2_period2A;    
    F_dff_length2_period2_raw = F_dff_length2_period2;
    F_dff_length2_period3_raw = F_dff_length2_period3;    
    F_dff_length3_period1_raw = F_dff_length3_period1;
    F_dff_length3_period2A_raw = F_dff_length3_period2A;    
    F_dff_length3_period2_raw = F_dff_length3_period2;
    F_dff_length3_period3_raw = F_dff_length3_period3;
    F_dff_length4_period1_raw = F_dff_length4_period1;
    F_dff_length4_period2A_raw = F_dff_length4_period2A;    
    F_dff_length4_period2_raw = F_dff_length4_period2;
    F_dff_length4_period3_raw = F_dff_length4_period3;
    

    F_dff_length1_period1 = F_dff_length1_period1_raw(:,validIndex_trialIndex_length1,:);
    F_dff_length1_period2A = F_dff_length1_period2A_raw(:,validIndex_trialIndex_length1,:);        
    F_dff_length1_period2 = F_dff_length1_period2_raw(:,validIndex_trialIndex_length1,:);  
    F_dff_length1_period3 = F_dff_length1_period3_raw(:,validIndex_trialIndex_length1,:);
                            
    F_dff_length2_period1 = F_dff_length2_period1_raw(:,validIndex_trialIndex_length2,:);    
    F_dff_length2_period2A = F_dff_length2_period2A_raw(:,validIndex_trialIndex_length2,:);            
    F_dff_length2_period2 = F_dff_length2_period2_raw(:,validIndex_trialIndex_length2,:);    
    F_dff_length2_period3 = F_dff_length2_period3_raw(:,validIndex_trialIndex_length2,:);
    
    F_dff_length3_period1 = F_dff_length3_period1_raw(:,validIndex_trialIndex_length3,:);
    F_dff_length3_period2A = F_dff_length3_period2A_raw(:,validIndex_trialIndex_length3,:);            
    F_dff_length3_period2 = F_dff_length3_period2_raw(:,validIndex_trialIndex_length3,:);    
    F_dff_length3_period3 = F_dff_length3_period3_raw(:,validIndex_trialIndex_length3,:);
    
    F_dff_length4_period1 = F_dff_length4_period1_raw(:,validIndex_trialIndex_length4,:);
    F_dff_length4_period2A = F_dff_length4_period2A_raw(:,validIndex_trialIndex_length4,:);            
    F_dff_length4_period2 = F_dff_length4_period2_raw(:,validIndex_trialIndex_length4,:);    
    F_dff_length4_period3 = F_dff_length4_period3_raw(:,validIndex_trialIndex_length4,:);
    
        
    %% perform one-way anova to test cell selectivity
    % sequence_length1_memoryCorrect
    sequence_length1_memoryCorrect = sequence_length1(validIndex_trialIndex_length1);
    sequence_length1_memoryCorrect_onehot = false(size(sequence_length1_memoryCorrect,2),numFrames);
    for tempi=1:size(sequence_length1_memoryCorrect,2)
        sequence_length1_memoryCorrect_onehot(tempi,sequence_length1_memoryCorrect(tempi)) = true;
    end
    
    % sequence_length2_memoryCorrect
    sequence_length2_memoryCorrect = sequence_length2(:,validIndex_trialIndex_length2);
    sequence_length2_memoryCorrect_onehot = false(size(sequence_length2_memoryCorrect,2),numFrames);
    for tempi=1:size(sequence_length2_memoryCorrect,2)
        sequence_length2_memoryCorrect_onehot(tempi,sequence_length2_memoryCorrect(:,tempi)') = true;
    end
    seqSetIndex_length2_memoryCorrect = zeros(1,sum(validIndex_trialIndex_length2));
    tempj = 0;
    for tempi=1:size(sequence_length2,2)
        if validIndex_trialIndex_length2(tempi) == true
            tempj = tempj + 1;
            temp_seq = sequence_length2(:,tempi)';
            temp_index = 0;
            for tempk=1:numSeq(2)
                if sum(ismember(temp_seq,target_seqSet{2}{tempk})) == 2
                    temp_index = tempk;
                    break
                end
            end
            temp_index;
            seqSetIndex_length2_memoryCorrect(tempj) = temp_index;
        end
    end
    
    % sequence_length3_memoryCorrect
    sequence_length3_memoryCorrect = sequence_length3(:,validIndex_trialIndex_length3);
    sequence_length3_memoryCorrect_onehot = false(size(sequence_length3_memoryCorrect,2),numFrames);
    for tempi=1:size(sequence_length3_memoryCorrect,2)
        sequence_length3_memoryCorrect_onehot(tempi,sequence_length3_memoryCorrect(:,tempi)') = true;
    end
    seqSetIndex_length3_memoryCorrect = zeros(1,sum(validIndex_trialIndex_length3));
    tempj = 0;
    for tempi=1:size(sequence_length3,2)
        if validIndex_trialIndex_length3(tempi) == true
            tempj = tempj + 1;
            temp_seq = sequence_length3(:,tempi)';
            temp_index = 0;
            for tempk=1:numSeq(3)
                if sum(ismember(temp_seq,target_seqSet{3}{tempk})) == 3
                    temp_index = tempk;
                    break
                end
            end
            temp_index;
            seqSetIndex_length3_memoryCorrect(tempj) = temp_index;
        end
    end
    
    % sequence_length4_memoryCorrect
    sequence_length4_memoryCorrect = sequence_length4(:,validIndex_trialIndex_length4);
    sequence_length4_memoryCorrect_onehot = false(size(sequence_length4_memoryCorrect,2),numFrames);
    for tempi=1:size(sequence_length4_memoryCorrect,2)
        sequence_length4_memoryCorrect_onehot(tempi,sequence_length4_memoryCorrect(:,tempi)') = true;
    end
    seqSetIndex_length4_memoryCorrect = zeros(1,sum(validIndex_trialIndex_length4));
    tempj = 0;
    for tempi=1:size(sequence_length4,2)
        if validIndex_trialIndex_length4(tempi) == true
            tempj = tempj + 1;
            temp_seq = sequence_length4(:,tempi)';
            temp_index = 0;
            for tempk=1:numSeq(4)
                if sum(ismember(temp_seq,target_seqSet{4}{tempk})) == 4
                    temp_index = tempk;
                    break
                end
            end
            temp_index;
            seqSetIndex_length4_memoryCorrect(tempj) = temp_index;
        end
    end
    
    
    
    %% Get F_dff_decisionPeriod
    decisionPeriodA_frameIndex = zeros(sum(trialIndex_bool_choiceMemory),3);    
    decisionPeriod_frameIndex = zeros(sum(trialIndex_bool_choiceMemory),3);
    for tempi=1:trial_para.trial_count
        temp_trialIndex = tempi;
        temp_markers = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers;
        
        for tempj=1:length(temp_markers)
            if strcmp(temp_markers{tempj}, 'DELAY 1 OFF') == 1
                temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(tempj);
                break
            end
        end
        decisionPeriodA_frameIndex(tempi,:) = [temp_frameIndex-33 temp_frameIndex temp_frameIndex+15];
        
        
        for tempj=1:length(temp_markers)
            if strcmp(temp_markers{tempj}, 'SELECTING AND DELAY 2 ON') == 1
                temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(tempj);
                break
            end
        end
        decisionPeriod_frameIndex(tempi,:) = [temp_frameIndex-15 temp_frameIndex temp_frameIndex+45];
    end
    decisionPeriodA_interval = decisionPeriodA_frameIndex(1,:)-decisionPeriodA_frameIndex(1,1)+1;    
    decisionPeriod_interval = decisionPeriod_frameIndex(1,:)-decisionPeriod_frameIndex(1,1)+1;
    
    [decisionPeriod_frameRangeMin,~] = bounds(decisionPeriod_frameIndex,2);
    temp_frameRange = repmat(decisionPeriod_frameRangeMin,1,decisionPeriod_interval(end)-decisionPeriod_interval(1)+1);
    decisionPeriod_frameRange = temp_frameRange + ((decisionPeriod_interval(1):decisionPeriod_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,decisionPeriod_frameRange'),[size(F_dff,1),size(decisionPeriod_frameRange')]);
    F_dff_decisionPeriod = permute(temp_dff,[1 3 2]);
        
    
    [decisionPeriodA_frameRangeMin,~] = bounds(decisionPeriodA_frameIndex,2);
    temp_frameRange = repmat(decisionPeriodA_frameRangeMin,1,decisionPeriodA_interval(end)-decisionPeriodA_interval(1)+1);
    decisionPeriodA_frameRange = temp_frameRange + ((decisionPeriodA_interval(1):decisionPeriodA_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,decisionPeriodA_frameRange'),[size(F_dff,1),size(decisionPeriodA_frameRange')]);
    F_dff_decisionPeriodA = permute(temp_dff,[1 3 2]);
    
    %F_dff_decisionBin = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(1):decisionPeriod_interval(2)),3);
    %F_dff_decision = F_dff_decisionPeriod(:,:,decisionPeriod_interval(1):decisionPeriod_interval(2));    
    %F_dff_delay1Bin = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(1):decisionPeriod_interval(2)),3);
    %F_dff_delay2Bin = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(2):decisionPeriod_interval(3)),3);
    
    F_dff_decisionBin = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
    F_dff_decision = F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));    
    F_dff_delay1Bin = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
    F_dff_delay2Bin = mean(F_dff_decisionPeriod(:,:,decisionPeriod_interval(2):decisionPeriod_interval(3)),3);
    
    a = 1;
    F_dff_baseline = F_dff_baselinePeriod(:,:,baselinePeriod_interval(1):baselinePeriod_interval(3));
    F_dff_delay1 = F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
    F_dff_delay2 = F_dff_decisionPeriod(:,:,decisionPeriod_interval(2):decisionPeriod_interval(3));
    
    %     %% Get F_dff_baselinePeriod
    %     baselinePeriod_frameIndex = zeros(sum(trialIndex_bool_choiceMemory),3);
    %     for tempi=1:trial_para.trial_count
    %         temp_trialIndex = tempi;
    %         temp_markers = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers;
    %         temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(1);
    %         baselinePeriod_frameIndex(tempi,:) = [temp_frameIndex-21 temp_frameIndex temp_frameIndex+12];
    %         %baselinePeriod_frameIndex(tempi,:) = [temp_frameIndex-12 temp_frameIndex temp_frameIndex+12];
    %     end
    %     baselinePeriod_interval = baselinePeriod_frameIndex(1,:)-baselinePeriod_frameIndex(1,1)+1;
    %
    %     [baselinePeriod_frameRangeMin,~] = bounds(baselinePeriod_frameIndex,2);
    %     temp_frameRange = repmat(baselinePeriod_frameRangeMin,1,baselinePeriod_interval(end)-baselinePeriod_interval(1)+1);
    %     baselinePeriod_frameRange = temp_frameRange + ((baselinePeriod_interval(1):baselinePeriod_interval(end))-1);
    %
    %     temp_dff = reshape(F_dff(:,baselinePeriod_frameRange'),[size(F_dff,1),size(baselinePeriod_frameRange')]);
    %     F_dff_baselinePeriod = permute(temp_dff,[1 3 2]);
    %
    %     F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,baselinePeriod_interval(1):baselinePeriod_interval(3)),3);
    
    if true
        %% Get F_dff_decisionPeriodB (forward delay1)
        decisionPeriodB_frameIndex = zeros(sum(trialIndex_bool_choiceMemory),2);
        for tempi=1:trial_para.trial_count
            temp_trialIndex = tempi;
            temp_markers = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers;
            
            for tempj=1:length(temp_markers)
                if strcmp(temp_markers{tempj}, 'DELAY 1 ON') == 1
                    temp_frameIndex = markerParse_trialLevel{temp_trialIndex}.currentTrialTotalMarkers_frameIndex(tempj);
                    break
                end
            end
            decisionPeriodB_frameIndex(tempi,:) = [temp_frameIndex+12 temp_frameIndex+45];
                        
        end
        decisionPeriodB_interval = decisionPeriodB_frameIndex(1,:)-decisionPeriodB_frameIndex(1,1)+1;
        
        
        [decisionPeriodB_frameRangeMin,~] = bounds(decisionPeriodB_frameIndex,2);
        temp_frameRange = repmat(decisionPeriodB_frameRangeMin,1,decisionPeriodB_interval(end)-decisionPeriodB_interval(1)+1);
        decisionPeriodB_frameRange = temp_frameRange + ((decisionPeriodB_interval(1):decisionPeriodB_interval(end))-1);
        
        temp_dff = reshape(F_dff(:,decisionPeriodB_frameRange'),[size(F_dff,1),size(decisionPeriodB_frameRange')]);
        F_dff_decisionPeriodB = permute(temp_dff,[1 3 2]);
        
    end
    
    %% Anova
    if if_compute == 1
        
        frameAnova1Name_v = autoGetFunName_myScripts('frameAnova1', [targetPATH '\functions']);
        fun_frameAnova1 = str2func(frameAnova1Name_v);
        
        a = 1;
        
        % selectiveCellNum_length1_TBin
        F_dff_length1_T1 = F_dff_length1_period1(:,:,length1_period1_interval(3):(length1_period1_interval(3)+17));
        trialLabel = sequence_length1_memoryCorrect;
        [selectiveCellBoolIndex_length1_TBin,selectiveCellNum_length1_TBin] = ...
            fun_frameAnova1(F_dff_length1_T1,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{1},if_selective_seq0_loc1);
        
        a = 1;
        
        % selectiveCellNum_length2_TBin
        trialLabel = seqSetIndex_length2_memoryCorrect;
        F_dff_length2_T1 = F_dff_length2_period1(:,:,length2_period1_interval(3):(length2_period1_interval(3)+17));
        [selectiveCellBoolIndex_length2_T1Bin,selectiveCellNum_length2_T1Bin] = ...
            fun_frameAnova1(F_dff_length2_T1,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{2},if_selective_seq0_loc1);
        F_dff_length2_T2 = F_dff_length2_period1(:,:,length2_period1_interval(4):(length2_period1_interval(4)+17));
        [selectiveCellBoolIndex_length2_T2Bin,selectiveCellNum_length2_T2Bin] = ...
            fun_frameAnova1(F_dff_length2_T2,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{2},if_selective_seq0_loc1);
        
        selectiveCellBoolIndex_length2_TBin = selectiveCellBoolIndex_length2_T1Bin | selectiveCellBoolIndex_length2_T2Bin;
        % selectiveCellBoolIndex_length2_TBin = selectiveCellBoolIndex_length2_T1Bin & selectiveCellBoolIndex_length2_T2Bin;
        selectiveCellNum_length2_TBin = sum(selectiveCellBoolIndex_length2_TBin);
        
        
        % selectiveCellNum_length3_TBin
        trialLabel = seqSetIndex_length3_memoryCorrect;
        F_dff_length3_T1 = F_dff_length3_period1(:,:,length3_period1_interval(3):(length3_period1_interval(3)+17));
        [selectiveCellBoolIndex_length3_T1Bin,selectiveCellNum_length3_T1Bin] = ...
            fun_frameAnova1(F_dff_length3_T1,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{3},if_selective_seq0_loc1);
        F_dff_length3_T2 = F_dff_length3_period1(:,:,length3_period1_interval(4):(length3_period1_interval(4)+17));
        [selectiveCellBoolIndex_length3_T2Bin,selectiveCellNum_length3_T2Bin] = ...
            fun_frameAnova1(F_dff_length3_T2,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{3},if_selective_seq0_loc1);
        F_dff_length3_T3 = F_dff_length3_period1(:,:,length3_period1_interval(5):(length3_period1_interval(5)+17));
        [selectiveCellBoolIndex_length3_T3Bin,selectiveCellNum_length3_T3Bin] = ...
            fun_frameAnova1(F_dff_length3_T3,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{3},if_selective_seq0_loc1);
        
        % selectiveCellBoolIndex_length3_TBin = selectiveCellBoolIndex_length3_T1Bin;
        selectiveCellBoolIndex_length3_TBin = selectiveCellBoolIndex_length3_T1Bin | selectiveCellBoolIndex_length3_T2Bin | selectiveCellBoolIndex_length3_T3Bin;
        % selectiveCellBoolIndex_length3_TBin = selectiveCellBoolIndex_length3_T1Bin & selectiveCellBoolIndex_length3_T2Bin & selectiveCellBoolIndex_length3_T3Bin;
        selectiveCellNum_length3_TBin = sum(selectiveCellBoolIndex_length3_TBin);
        
        
        % selectiveCellNum_length4_TBin
        trialLabel = seqSetIndex_length4_memoryCorrect;
        F_dff_length4_T1 = F_dff_length4_period1(:,:,length4_period1_interval(3):(length4_period1_interval(3)+17));
        [selectiveCellBoolIndex_length4_T1Bin,selectiveCellNum_length4_T1Bin] = ...
            fun_frameAnova1(F_dff_length4_T1,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{4},if_selective_seq0_loc1);
        F_dff_length4_T2 = F_dff_length4_period1(:,:,length4_period1_interval(4):(length4_period1_interval(4)+17));
        [selectiveCellBoolIndex_length4_T2Bin,selectiveCellNum_length4_T2Bin] = ...
            fun_frameAnova1(F_dff_length4_T2,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{4},if_selective_seq0_loc1);
        F_dff_length4_T3 = F_dff_length4_period1(:,:,length4_period1_interval(5):(length4_period1_interval(5)+17));
        [selectiveCellBoolIndex_length4_T3Bin,selectiveCellNum_length4_T3Bin] = ...
            fun_frameAnova1(F_dff_length4_T3,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{4},if_selective_seq0_loc1);
        F_dff_length4_T4 = F_dff_length4_period1(:,:,length4_period1_interval(6):(length4_period1_interval(6)+17));
        [selectiveCellBoolIndex_length4_T4Bin,selectiveCellNum_length4_T4Bin] = ...
            fun_frameAnova1(F_dff_length4_T4,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{4},if_selective_seq0_loc1);
        
        selectiveCellBoolIndex_length4_TBin = selectiveCellBoolIndex_length4_T1Bin | selectiveCellBoolIndex_length4_T2Bin | selectiveCellBoolIndex_length4_T3Bin | selectiveCellBoolIndex_length4_T4Bin;
        selectiveCellNum_length4_TBin = sum(selectiveCellBoolIndex_length4_TBin);

        
        a = 1;
        
        % selectiveCellNum_length1_delay1Bin
        if if_delay1Selective_old0_new1 == 0
            F_dff_length1_delay1 = F_dff_length1_period1(:,:,length1_period1_interval(3):length1_period1_interval(end));
        elseif if_delay1Selective_old0_new1 == 1
            F_dff_length1_delay1 = F_dff_decision(:,trialIndex_lengthx_bool_memoryCorrect(1,:),:);      
        end
        trialLabel = sequence_length1_memoryCorrect;
        [selectiveCellBoolIndex_length1_delay1Bin,selectiveCellNum_length1_delay1Bin] = ...
            fun_frameAnova1(F_dff_length1_delay1,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{1},if_selective_seq0_loc1);
        
        % selectiveCellNum_length1_delay2Bin
        F_dff_length1_delay2 = F_dff_length1_period2(:,:,length1_period2_interval(2):length1_period2_interval(3));        
        trialLabel = sequence_length1_memoryCorrect;
        [selectiveCellBoolIndex_length1_delay2Bin,selectiveCellNum_length1_delay2Bin] = ...
            fun_frameAnova1(F_dff_length1_delay2,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{1},if_selective_seq0_loc1);
        
        % selectiveCellNum_length2_delay1Bin
        if if_delay1Selective_old0_new1 == 0
            F_dff_length2_delay1 = F_dff_length2_period1(:,:,length2_period1_interval(3):length2_period1_interval(end));
        elseif if_delay1Selective_old0_new1 == 1
            F_dff_length2_delay1 = F_dff_decision(:,trialIndex_lengthx_bool_memoryCorrect(2,:),:);      
        end        
        trialLabel = seqSetIndex_length2_memoryCorrect;
        [selectiveCellBoolIndex_length2_delay1Bin,selectiveCellNum_length2_delay1Bin] = ...
            fun_frameAnova1(F_dff_length2_delay1,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{2},if_selective_seq0_loc1);
        
        % selectiveCellNum_length2_delay2Bin
        F_dff_length2_delay2 = F_dff_length2_period2(:,:,length2_period2_interval(2):length2_period2_interval(3));
        trialLabel = seqSetIndex_length2_memoryCorrect;
        [selectiveCellBoolIndex_length2_delay2Bin,selectiveCellNum_length2_delay2Bin] = ...
            fun_frameAnova1(F_dff_length2_delay2,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{2},if_selective_seq0_loc1);
        
        % selectiveCellNum_length3_delay1Bin
        if if_delay1Selective_old0_new1 == 0
            F_dff_length3_delay1 = F_dff_length3_period1(:,:,length3_period1_interval(3):length3_period1_interval(end));
        elseif if_delay1Selective_old0_new1 == 1
            F_dff_length3_delay1 = F_dff_decision(:,trialIndex_lengthx_bool_memoryCorrect(3,:),:);      
        end         
        trialLabel = seqSetIndex_length3_memoryCorrect;
        [selectiveCellBoolIndex_length3_delay1Bin,selectiveCellNum_length3_delay1Bin] = ...
            fun_frameAnova1(F_dff_length3_delay1,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{3},if_selective_seq0_loc1);
        
        % selectiveCellNum_length3_delay2Bin
        F_dff_length3_delay2 = F_dff_length3_period2(:,:,length3_period2_interval(2):length3_period2_interval(3));
        trialLabel = seqSetIndex_length3_memoryCorrect;
        [selectiveCellBoolIndex_length3_delay2Bin,selectiveCellNum_length3_delay2Bin] = ...
            fun_frameAnova1(F_dff_length3_delay2,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{3},if_selective_seq0_loc1);
        
        % selectiveCellNum_length4_delay1Bin
        if if_delay1Selective_old0_new1 == 0
            F_dff_length4_delay1 = F_dff_length4_period1(:,:,length4_period1_interval(3):length4_period1_interval(end));
        elseif if_delay1Selective_old0_new1 == 1
            F_dff_length4_delay1 = F_dff_decision(:,trialIndex_lengthx_bool_memoryCorrect(4,:),:);      
        end                 
        trialLabel = seqSetIndex_length4_memoryCorrect;
        [selectiveCellBoolIndex_length4_delay1Bin,selectiveCellNum_length4_delay1Bin] = ...
            fun_frameAnova1(F_dff_length4_delay1,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{4},if_selective_seq0_loc1);
        
        % selectiveCellNum_length4_delay2Bin
        F_dff_length4_delay2 = F_dff_length4_period2(:,:,length4_period2_interval(2):length4_period2_interval(3));
        trialLabel = seqSetIndex_length4_memoryCorrect;
        [selectiveCellBoolIndex_length4_delay2Bin,selectiveCellNum_length4_delay2Bin] = ...
            fun_frameAnova1(F_dff_length4_delay2,trialLabel,significantThreshold,significantThreshold_lowBound,target_seqSet{4},if_selective_seq0_loc1);
        
        
        % selectiveCellNum_decisionBin
        F_dff_decisionBin_choiceMemory = F_dff_decisionPeriodA_choiceMemory(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        F_dff_decisionBin_choiceOffload = F_dff_decisionPeriodA_choiceOffload(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        F_dff_decisionBin_noChoice = F_dff_decisionPeriodA_noChoice(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
        
        % selectiveCellBoolIndex_decisionBin_MO
        F_dff_decisionBin_MO = [F_dff_decisionBin_choiceMemory, F_dff_decisionBin_choiceOffload];
        temp_context = zeros(1,size(F_dff_decisionBin_MO,2));
        temp_context(1:size(F_dff_decisionBin_choiceMemory,2)) = 1;
        temp_context(size(F_dff_decisionBin_choiceMemory,2)+...
            (1:size(F_dff_decisionBin_choiceOffload,2))) = 2;
        trialLabel = temp_context;
        [selectiveCellBoolIndex_decisionBin_MO,selectiveCellNum_decisionBin_MO] = ...
            fun_frameAnova1(F_dff_decisionBin_MO,trialLabel,significantThreshold,significantThreshold_lowBound);
        selectiveCellIndex_suite2p_decisionBin_MO = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_MO);
        
        % selectiveCellBoolIndex_decisionBin_MN
        F_dff_decisionBin_MN = [F_dff_decisionBin_choiceMemory, F_dff_decisionBin_noChoice];
        temp_context = zeros(1,size(F_dff_decisionBin_MN,2));
        temp_context(1:size(F_dff_decisionBin_choiceMemory,2)) = 1;
        temp_context(size(F_dff_decisionBin_choiceMemory,2)+...
            (1:size(F_dff_decisionBin_noChoice,2))) = 3;
        trialLabel = temp_context;
        [selectiveCellBoolIndex_decisionBin_MN,selectiveCellNum_decisionBin_MN] = ...
            fun_frameAnova1(F_dff_decisionBin_MN,trialLabel,significantThreshold,significantThreshold_lowBound);
        
        % selectiveCellBoolIndex_decisionBin_ON
        F_dff_decisionBin_ON = [F_dff_decisionBin_choiceOffload, F_dff_decisionBin_noChoice];
        temp_context = zeros(1,size(F_dff_decisionBin_ON,2));
        temp_context(1:size(F_dff_decisionBin_choiceOffload,2)) = 2;
        temp_context(size(F_dff_decisionBin_choiceOffload,2)+...
            (1:size(F_dff_decisionBin_noChoice,2))) = 3;
        trialLabel = temp_context;
        [selectiveCellBoolIndex_decisionBin_ON,selectiveCellNum_decisionBin_ON] = ...
            fun_frameAnova1(F_dff_decisionBin_ON,trialLabel,significantThreshold,significantThreshold_lowBound);
        
        % selectiveCellBoolIndex_decisionBin_andMON
        selectiveCellBoolIndex_decisionBin_andMON = selectiveCellBoolIndex_decisionBin_MO & selectiveCellBoolIndex_decisionBin_MN & selectiveCellBoolIndex_decisionBin_ON;
        selectiveCellNum_decisionBin_andMON = sum(selectiveCellBoolIndex_decisionBin_andMON);
        selectiveCellIndex_suite2p_decisionBin_andMON = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_andMON);
        % fprintf('selectiveCellNum_decisionBin_andMON (choiceMemory & choiceOffload & noChoice) = %d.\n',selectiveCellNum_decisionBin_andMON);
        % selectiveCellBoolIndex_decisionBin_orMON
        selectiveCellBoolIndex_decisionBin_orMON = selectiveCellBoolIndex_decisionBin_MO | selectiveCellBoolIndex_decisionBin_MN | selectiveCellBoolIndex_decisionBin_ON;
        selectiveCellNum_decisionBin_orMON = sum(selectiveCellBoolIndex_decisionBin_orMON);
        selectiveCellIndex_suite2p_decisionBin_orMON = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_orMON);
        % fprintf('selectiveCellNum_decisionBin_orMON (choiceMemory | choiceOffload | noChoice) = %d.\n',selectiveCellNum_decisionBin_orMON);
        fprintf('selectiveCellNum_decisionBin_MON (and,or) = [%d,%d].\n',selectiveCellNum_decisionBin_andMON,selectiveCellNum_decisionBin_orMON);
        
        a = 1;
        
        
        
        % selectiveCellBoolIndex_baselineBin_MO
        F_dff_baselinePeriod_choiceMemory = F_dff_baselinePeriod(:,trialIndex_bool_choiceMemory,:);
        F_dff_baselinePeriod_choiceOffload = F_dff_baselinePeriod(:,trialIndex_bool_choiceOffload,:);
        F_dff_baselineBin_MO = [F_dff_baselinePeriod_choiceMemory, F_dff_baselinePeriod_choiceOffload];
        temp_context = zeros(1,size(F_dff_baselineBin_MO,2));
        temp_context(1:size(F_dff_baselinePeriod_choiceMemory,2)) = 1;
        temp_context(size(F_dff_baselinePeriod_choiceMemory,2)+...
            (1:size(F_dff_baselinePeriod_choiceOffload,2))) = 2;
        trialLabel = temp_context;
        [selectiveCellBoolIndex_baselineBin_MO,selectiveCellNum_baselineBin_MO] = ...
            fun_frameAnova1(F_dff_baselineBin_MO,trialLabel,significantThreshold,significantThreshold_lowBound);
        selectiveCellIndex_suite2p_baselineBin_MO = cellIndex_suite2p(selectiveCellBoolIndex_baselineBin_MO);
        fprintf('selectiveCellNum_baselineBin_MO = %d.\n',selectiveCellNum_baselineBin_MO);
        
        
        % selectiveCellBoolIndex_baselineBin_afterMemoryCorrectError
        F_dff_baselinePeriod_afterMemoryCorrect = F_dff_baselinePeriod(:,trialIndex_bool_afterMemoryCorrect,:);
        F_dff_baselinePeriod_afterMemoryError = F_dff_baselinePeriod(:,trialIndex_bool_afterMemoryError,:);
        F_dff_baselineBin_afterMemoryCorrectError = [F_dff_baselinePeriod_afterMemoryCorrect, F_dff_baselinePeriod_afterMemoryError];
        temp_context = zeros(1,size(F_dff_baselineBin_afterMemoryCorrectError,2));
        temp_context(1:size(F_dff_baselinePeriod_afterMemoryCorrect,2)) = 1;
        temp_context(size(F_dff_baselinePeriod_afterMemoryCorrect,2)+...
            (1:size(F_dff_baselinePeriod_afterMemoryError,2))) = 2;
        trialLabel = temp_context;
        [selectiveCellBoolIndex_baselineBin_afterMemoryCorrectError,selectiveCellNum_baselineBin_afterMemoryCorrectError] = ...
            fun_frameAnova1(F_dff_baselineBin_afterMemoryCorrectError,trialLabel,significantThreshold,significantThreshold_lowBound);
        selectiveCellIndex_suite2p_baselineBin_afterMemoryCorrectError = cellIndex_suite2p(selectiveCellBoolIndex_baselineBin_afterMemoryCorrectError);
        fprintf('selectiveCellNum_baselineBin_afterMemoryCorrectError = %d.\n',selectiveCellNum_baselineBin_afterMemoryCorrectError);
        
        
        a = 1;
        
        % selectiveCellBoolIndex_baselineBin_afterMemoryCorrect_MO
        F_dff_baselinePeriod_afterMemoryCorrect_choiceMemory = F_dff_baselinePeriod(:,trialIndex_bool_choiceMemory&trialIndex_bool_afterMemoryCorrect,:);
        F_dff_baselinePeriod_afterMemoryCorrect_choiceOffload = F_dff_baselinePeriod(:,trialIndex_bool_choiceOffload&trialIndex_bool_afterMemoryCorrect,:);
        F_dff_baselineBin_afterMemoryCorrect_MO = [F_dff_baselinePeriod_afterMemoryCorrect_choiceMemory, F_dff_baselinePeriod_afterMemoryCorrect_choiceOffload];
        temp_context = zeros(1,size(F_dff_baselineBin_afterMemoryCorrect_MO,2));
        temp_context(1:size(F_dff_baselinePeriod_afterMemoryCorrect_choiceMemory,2)) = 1;
        temp_context(size(F_dff_baselinePeriod_afterMemoryCorrect_choiceMemory,2)+...
            (1:size(F_dff_baselinePeriod_afterMemoryCorrect_choiceOffload,2))) = 2;
        trialLabel = temp_context;
        [selectiveCellBoolIndex_baselineBin_afterMemoryCorrect_MO,selectiveCellNum_baselineBin_afterMemoryCorrect_MO] = ...
            fun_frameAnova1(F_dff_baselineBin_afterMemoryCorrect_MO,trialLabel,significantThreshold,significantThreshold_lowBound);
        selectiveCellIndex_suite2p_baselineBin_afterMemoryCorrect_MO = cellIndex_suite2p(selectiveCellBoolIndex_baselineBin_afterMemoryCorrect_MO);
        fprintf('selectiveCellNum_baselineBin_afterMemoryCorrect_MO = %d.\n',selectiveCellNum_baselineBin_afterMemoryCorrect_MO);
        
        
        % selectiveCellBoolIndex_baselineBin_afterMemoryError_MO
        F_dff_baselinePeriod_afterMemoryError_choiceMemory = F_dff_baselinePeriod(:,trialIndex_bool_choiceMemory&trialIndex_bool_afterMemoryError,:);
        F_dff_baselinePeriod_afterMemoryError_choiceOffload = F_dff_baselinePeriod(:,trialIndex_bool_choiceOffload&trialIndex_bool_afterMemoryError,:);
        F_dff_baselineBin_afterMemoryError_MO = [F_dff_baselinePeriod_afterMemoryError_choiceMemory, F_dff_baselinePeriod_afterMemoryError_choiceOffload];
        temp_context = zeros(1,size(F_dff_baselineBin_afterMemoryError_MO,2));
        temp_context(1:size(F_dff_baselinePeriod_afterMemoryError_choiceMemory,2)) = 1;
        temp_context(size(F_dff_baselinePeriod_afterMemoryError_choiceMemory,2)+...
            (1:size(F_dff_baselinePeriod_afterMemoryError_choiceOffload,2))) = 2;
        trialLabel = temp_context;
        [selectiveCellBoolIndex_baselineBin_afterMemoryError_MO,selectiveCellNum_baselineBin_afterMemoryError_MO] = ...
            fun_frameAnova1(F_dff_baselineBin_afterMemoryError_MO,trialLabel,significantThreshold,significantThreshold_lowBound);
        selectiveCellIndex_suite2p_baselineBin_afterMemoryError_MO = cellIndex_suite2p(selectiveCellBoolIndex_baselineBin_afterMemoryError_MO);
        fprintf('selectiveCellNum_baselineBin_afterMemoryError_MO = %d.\n',selectiveCellNum_baselineBin_afterMemoryError_MO);
        
        
        % Get temp_trialIndex_bool_lengthxSeq in MON
        temp_pointKindsNum = 4;
        choice_BoolFlag;
        ifSelectOffloading_bool;
        sequence = trial_para.currentSequence;
        
        temp_trialIndex_bool_lengthxSeq_choiceMemory = cell(1,temp_pointKindsNum);
        temp_trialIndex_bool_lengthxSeq_choiceOffload = cell(1,temp_pointKindsNum);
        temp_trialIndex_bool_lengthxSeq_noChoice = cell(1,temp_pointKindsNum);
        for temp_lengthFlag=1:temp_pointKindsNum
            temp_trialIndex_bool_lengthxSeq_choiceMemory{temp_lengthFlag} = false(numSeq(temp_lengthFlag),trial_para.trial_count);
            temp_trialIndex_bool_lengthxSeq_choiceOffload{temp_lengthFlag} = false(numSeq(temp_lengthFlag),trial_para.trial_count);
            temp_trialIndex_bool_lengthxSeq_noChoice{temp_lengthFlag} = false(numSeq(temp_lengthFlag),trial_para.trial_count);
            for tempi=1:trial_para.trial_count
                currentSequence = trial_para.currentSequence{tempi};
                temp_seq_length = length(currentSequence);
                if temp_seq_length ~= temp_lengthFlag
                    continue
                end
                for tempj=1:numSeq(temp_lengthFlag)
                    if sum(ismember(currentSequence,target_seqSet{temp_lengthFlag}{tempj})) == temp_lengthFlag
                        if choice_BoolFlag(tempi) == true
                            if ifSelectOffloading_bool(tempi) == true
                                temp_trialIndex_bool_lengthxSeq_choiceOffload{temp_lengthFlag}(tempj,tempi) = true;
                            elseif ifSelectOffloading_bool(tempi) == false
                                temp_trialIndex_bool_lengthxSeq_choiceMemory{temp_lengthFlag}(tempj,tempi) = true;
                            end
                        elseif choice_BoolFlag(tempi) == false
                            temp_trialIndex_bool_lengthxSeq_noChoice{temp_lengthFlag}(tempj,tempi) = true;
                        end
                    end
                end
            end
        end
        
        a = 1;
        
        % Compare each lengthxSeq of F_dff_decisionPeriod in MON
        if if_compute_selectiveCellNum_decisionBin_lengthxSeq_MON == 1
            %t0 = tic;
            temp_binSize = 11;%15
            F_dff_decision_2 = F_dff_decision(:,:,1:floor(size(F_dff_decision,3)/temp_binSize)*temp_binSize);
            F_dff_decision_2_4d = reshape(F_dff_decision_2,[size(F_dff_decision_2,[1 2]) temp_binSize size(F_dff_decision_2,3)/temp_binSize]);
            F_dff_decision_2bin = squeeze(mean(F_dff_decision_2_4d,3));
            
            roiNum = size(F_dff_decisionBin,1);
            selectiveCellBoolIndex_decisionBin_lengthxSeq_MO = cell(1,temp_pointKindsNum);
            selectiveCellBoolIndex_decisionBin_lengthxSeq_MN = cell(1,temp_pointKindsNum);
            selectiveCellBoolIndex_decisionBin_lengthxSeq_ON = cell(1,temp_pointKindsNum);
            selectiveCellNum_decisionBin_lengthxSeq_MO = cell(1,temp_pointKindsNum);
            selectiveCellNum_decisionBin_lengthxSeq_MN = cell(1,temp_pointKindsNum);
            selectiveCellNum_decisionBin_lengthxSeq_ON = cell(1,temp_pointKindsNum);
            
            for temp_lengthFlag=1:temp_pointKindsNum
                selectiveCellBoolIndex_decisionBin_lengthxSeq_MO{temp_lengthFlag} = false(numSeq(temp_lengthFlag),roiNum);
                selectiveCellBoolIndex_decisionBin_lengthxSeq_MN{temp_lengthFlag} = false(numSeq(temp_lengthFlag),roiNum);
                selectiveCellBoolIndex_decisionBin_lengthxSeq_ON{temp_lengthFlag} = false(numSeq(temp_lengthFlag),roiNum);
                selectiveCellNum_decisionBin_lengthxSeq_MO{temp_lengthFlag} = zeros(numSeq(temp_lengthFlag),1);
                selectiveCellNum_decisionBin_lengthxSeq_MN{temp_lengthFlag} = zeros(numSeq(temp_lengthFlag),1);
                selectiveCellNum_decisionBin_lengthxSeq_ON{temp_lengthFlag} = zeros(numSeq(temp_lengthFlag),1);
                
                for tempi=1:numSeq(temp_lengthFlag)
                    temp_bool_lengthxSeq_M = temp_trialIndex_bool_lengthxSeq_choiceMemory{temp_lengthFlag}(tempi,:);
                    temp_bool_lengthxSeq_O = temp_trialIndex_bool_lengthxSeq_choiceOffload{temp_lengthFlag}(tempi,:);
                    temp_bool_lengthxSeq_N = temp_trialIndex_bool_lengthxSeq_noChoice{temp_lengthFlag}(tempi,:);
                    
                    temp_bool_lengthxSeq = temp_bool_lengthxSeq_M | temp_bool_lengthxSeq_O | temp_bool_lengthxSeq_N;
                    
                    temp_context_lengthxSeq = zeros(1, trial_para.trial_count);
                    temp_context_lengthxSeq(temp_bool_lengthxSeq_M) = 1;
                    temp_context_lengthxSeq(temp_bool_lengthxSeq_O) = 2;
                    temp_context_lengthxSeq(temp_bool_lengthxSeq_N) = 3;
                    temp_context_lengthxSeq_raw = temp_context_lengthxSeq;
                    temp_context_lengthxSeq = temp_context_lengthxSeq_raw(temp_context_lengthxSeq_raw > 0);
                    
                    temp_contextBool_lengthxSeq_MO = (temp_context_lengthxSeq ~= 3);
                    temp_contextBool_lengthxSeq_MN = (temp_context_lengthxSeq ~= 2);
                    temp_contextBool_lengthxSeq_ON = (temp_context_lengthxSeq ~= 1);
                    
                    %temp_dff_lengthxSeq = F_dff_decisionBin(:,temp_bool_lengthxSeq);
                    temp_dff_lengthxSeq = F_dff_decision_2bin(:,temp_bool_lengthxSeq,:);
                    
                    a = 1;
                    
                    temp_dff = temp_dff_lengthxSeq(:,temp_contextBool_lengthxSeq_MO,:);
                    trialLabel = temp_context_lengthxSeq(temp_contextBool_lengthxSeq_MO);
                    [selectiveCellBoolIndex_decisionBin_lengthxSeq_MO{temp_lengthFlag}(tempi,:),...
                        selectiveCellNum_decisionBin_lengthxSeq_MO{temp_lengthFlag}(tempi)] = ...
                        fun_frameAnova1(temp_dff,trialLabel,significantThreshold,significantThreshold_lowBound);
                    
                    temp_dff = temp_dff_lengthxSeq(:,temp_contextBool_lengthxSeq_MN,:);
                    trialLabel = temp_context_lengthxSeq(temp_contextBool_lengthxSeq_MN);
                    [selectiveCellBoolIndex_decisionBin_lengthxSeq_MN{temp_lengthFlag}(tempi,:),...
                        selectiveCellNum_decisionBin_lengthxSeq_MN{temp_lengthFlag}(tempi)] = ...
                        fun_frameAnova1(temp_dff,trialLabel,significantThreshold,significantThreshold_lowBound);
                    
                    temp_dff = temp_dff_lengthxSeq(:,temp_contextBool_lengthxSeq_ON,:);
                    trialLabel = temp_context_lengthxSeq(temp_contextBool_lengthxSeq_ON);
                    [selectiveCellBoolIndex_decisionBin_lengthxSeq_ON{temp_lengthFlag}(tempi,:),...
                        selectiveCellNum_decisionBin_lengthxSeq_ON{temp_lengthFlag}(tempi)] = ...
                        fun_frameAnova1(temp_dff,trialLabel,significantThreshold,significantThreshold_lowBound);
                    
                    %fprintf('t = %.0f secs.\n',toc(t0));
                    a = 1;
                end
            end
            
            selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MO = false(roiNum,temp_pointKindsNum);
            selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MN = false(roiNum,temp_pointKindsNum);
            selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_ON = false(roiNum,temp_pointKindsNum);
            selectiveCellNum_decisionBin_lengthxSeq_seqMerged_MO = zeros(1,temp_pointKindsNum);
            selectiveCellNum_decisionBin_lengthxSeq_seqMerged_MN = zeros(1,temp_pointKindsNum);
            selectiveCellNum_decisionBin_lengthxSeq_seqMerged_ON = zeros(1,temp_pointKindsNum);
            for temp_lengthFlag=1:temp_pointKindsNum
                selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MO(:,temp_lengthFlag) = sum(selectiveCellBoolIndex_decisionBin_lengthxSeq_MO{temp_lengthFlag},1) > 0;
                selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MN(:,temp_lengthFlag) = sum(selectiveCellBoolIndex_decisionBin_lengthxSeq_MN{temp_lengthFlag},1) > 0;
                selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_ON(:,temp_lengthFlag) = sum(selectiveCellBoolIndex_decisionBin_lengthxSeq_ON{temp_lengthFlag},1) > 0;
                selectiveCellNum_decisionBin_lengthxSeq_seqMerged_MO(temp_lengthFlag) = sum(selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MO(:,temp_lengthFlag));
                selectiveCellNum_decisionBin_lengthxSeq_seqMerged_MN(temp_lengthFlag) = sum(selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MN(:,temp_lengthFlag));
                selectiveCellNum_decisionBin_lengthxSeq_seqMerged_ON(temp_lengthFlag) = sum(selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_ON(:,temp_lengthFlag));
            end
            
            % selectiveCellBoolIndex_decisionBin_lengthxSeq_andMON
            selectiveCellBoolIndex_decisionBin_lengthxSeq_andMON = selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MO & selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MN & selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_ON;
            selectiveCellNum_decisionBin_lengthxSeq_andMON = sum(selectiveCellBoolIndex_decisionBin_lengthxSeq_andMON);
            selectiveCellIndex_suite2p_decisionBin_length1Seq_andMON = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_lengthxSeq_andMON(:,1));
            selectiveCellIndex_suite2p_decisionBin_length2Seq_andMON = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_lengthxSeq_andMON(:,2));
            selectiveCellIndex_suite2p_decisionBin_length3Seq_andMON = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_lengthxSeq_andMON(:,3));
            % fprintf('selectiveCellNum_decisionBin_lengthxSeq_andMON (choiceMemory & choiceOffload & noChoice) = [%d %d %d].\n',selectiveCellNum_decisionBin_lengthxSeq_andMON);
            
            % selectiveCellBoolIndex_decisionBin_lengthxSeq_orMON
            selectiveCellBoolIndex_decisionBin_lengthxSeq_orMON = selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MO | selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MN | selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_ON;
            selectiveCellNum_decisionBin_lengthxSeq_orMON = sum(selectiveCellBoolIndex_decisionBin_lengthxSeq_orMON);
            selectiveCellIndex_suite2p_decisionBin_length1Seq_orMON = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_lengthxSeq_orMON(:,1));
            selectiveCellIndex_suite2p_decisionBin_length2Seq_orMON = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_lengthxSeq_orMON(:,2));
            selectiveCellIndex_suite2p_decisionBin_length3Seq_orMON = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_lengthxSeq_orMON(:,3));
            % fprintf('selectiveCellNum_decisionBin_lengthxSeq_orMON (choiceMemory | choiceOffload | noChoice) = [%d %d %d].\n',selectiveCellNum_decisionBin_lengthxSeq_orMON);
            fprintf('selectiveCellNum_decisionBin_lengthxSeq_MON (and,or) = [%d %d %d],[%d %d %d].\n',selectiveCellNum_decisionBin_lengthxSeq_andMON,selectiveCellNum_decisionBin_lengthxSeq_orMON);
            
            % selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MN
            selectiveCellIndex_suite2p_decisionBin_length1Seq_seqMerged_MN = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MN(:,1));
            selectiveCellIndex_suite2p_decisionBin_length2Seq_seqMerged_MN = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MN(:,2));
            selectiveCellIndex_suite2p_decisionBin_length3Seq_seqMerged_MN = cellIndex_suite2p(selectiveCellBoolIndex_decisionBin_lengthxSeq_seqMerged_MN(:,3));
            fprintf('selectiveCellNum_decisionBin_lengthxSeq_seqMerged_MN = [%d %d %d].\n',selectiveCellNum_decisionBin_lengthxSeq_seqMerged_MN);
            
            a = 1;
            
        end
        
        % selectiveCellNum_postDecisionPeriod
        F_dff_postDecisionPeriod_choiceMemory = F_dff_decisionPeriod_choiceMemory(:,:,decisionPeriod_interval(2):decisionPeriod_interval(3));
        F_dff_postDecisionPeriod_noChoice = F_dff_decisionPeriod_noChoice(:,:,decisionPeriod_interval(2):decisionPeriod_interval(3));
        F_dff_postDecisionPeriod = [F_dff_postDecisionPeriod_choiceMemory, F_dff_postDecisionPeriod_noChoice];
        
        temp_ifChoiceCondition = false(1,size(F_dff_postDecisionPeriod,2));
        temp_ifChoiceCondition(1:size(F_dff_postDecisionPeriod_choiceMemory,2)) = true;
        trialLabel = temp_ifChoiceCondition;
        [selectiveCellBoolIndex_postDecisionPeriod,selectiveCellNum_postDecisionPeriod] = ...
            fun_frameAnova1(F_dff_postDecisionPeriod,trialLabel,significantThreshold,significantThreshold_lowBound);
        
        selectiveCellIndex_suite2p_postDecisionPeriod = cellIndex_suite2p(selectiveCellBoolIndex_postDecisionPeriod);
        % fprintf('selectiveCellNum_postDecisionPeriod (choice vs noChoice) = %d .\n',selectiveCellNum_postDecisionPeriod);
        
        selectiveCellBoolIndex_length1_T = selectiveCellBoolIndex_length1_TBin;
        selectiveCellIndex_suite2p_length1_T = cellIndex_suite2p(selectiveCellBoolIndex_length1_T);
        selectiveCellNum_length1_T = size(selectiveCellIndex_suite2p_length1_T,1);
        % fprintf('selectiveCellNum_length1_T = %d.\n',selectiveCellNum_length1_T);
        
        selectiveCellBoolIndex_length2_T = selectiveCellBoolIndex_length2_TBin;
        selectiveCellIndex_suite2p_length2_T = cellIndex_suite2p(selectiveCellBoolIndex_length2_T);
        selectiveCellNum_length2_T = size(selectiveCellIndex_suite2p_length2_T,1);
        % fprintf('selectiveCellNum_length2_T = %d.\n',selectiveCellNum_length2_T);
        
        selectiveCellBoolIndex_length3_T = selectiveCellBoolIndex_length3_TBin;
        selectiveCellIndex_suite2p_length3_T = cellIndex_suite2p(selectiveCellBoolIndex_length3_T);
        selectiveCellNum_length3_T = size(selectiveCellIndex_suite2p_length3_T,1);
        % fprintf('selectiveCellNum_length3_T = %d.\n',selectiveCellNum_length3_T);
        
        selectiveCellBoolIndex_length4_T = selectiveCellBoolIndex_length4_TBin;
        selectiveCellIndex_suite2p_length4_T = cellIndex_suite2p(selectiveCellBoolIndex_length4_T);
        selectiveCellNum_length4_T = size(selectiveCellIndex_suite2p_length4_T,1);
        % fprintf('selectiveCellNum_length4_T = %d.\n',selectiveCellNum_length4_T);        
        
        fprintf('selectiveCellNum_lengthx_T = [%d %d %d %d].\n',selectiveCellNum_length1_T,selectiveCellNum_length2_T,selectiveCellNum_length3_T,selectiveCellNum_length4_T);
        
        
        selectiveCellBoolIndex_length1_and = selectiveCellBoolIndex_length1_delay1Bin & selectiveCellBoolIndex_length1_delay2Bin;
        selectiveCellIndex_suite2p_length1_and = cellIndex_suite2p(selectiveCellBoolIndex_length1_and);
        selectiveCellNum_length1_and = size(selectiveCellIndex_suite2p_length1_and,1);
        % fprintf('selectiveCellNum_length1_and (delay1 and delay2) = %d.\n',selectiveCellNum_length1_and);
        
        selectiveCellBoolIndex_length1 = selectiveCellBoolIndex_length1_delay1Bin | selectiveCellBoolIndex_length1_delay2Bin;
        selectiveCellIndex_suite2p_length1 = cellIndex_suite2p(selectiveCellBoolIndex_length1);
        selectiveCellNum_length1 = size(selectiveCellIndex_suite2p_length1,1);
        fprintf('selectiveCellNum_length1 = %d (%d from delay1, %d from delay2, %d from both).\n',selectiveCellNum_length1,selectiveCellNum_length1_delay1Bin,selectiveCellNum_length1_delay2Bin,selectiveCellNum_length1_and);
        
        
        selectiveCellIndex_suite2p_length1_delay1Bin = cellIndex_suite2p(selectiveCellBoolIndex_length1_delay1Bin);
        selectiveCellNum_length1_delay1Bin = size(selectiveCellIndex_suite2p_length1_delay1Bin,1);
        % fprintf('selectiveCellNum_length1_delay1Bin = %d.\n',selectiveCellNum_length1_delay1Bin);
        
        selectiveCellIndex_suite2p_length1_delay2Bin = cellIndex_suite2p(selectiveCellBoolIndex_length1_delay2Bin);
        selectiveCellNum_length1_delay2Bin = size(selectiveCellIndex_suite2p_length1_delay2Bin,1);
        % fprintf('selectiveCellNum_length1_delay2Bin = %d.\n',selectiveCellNum_length1_delay2Bin);
        
        selectiveCellIndex_suite2p_length2_delay1Bin = cellIndex_suite2p(selectiveCellBoolIndex_length2_delay1Bin);
        selectiveCellIndex_suite2p_length2_delay2Bin = cellIndex_suite2p(selectiveCellBoolIndex_length2_delay2Bin);
        selectiveCellIndex_suite2p_length3_delay1Bin = cellIndex_suite2p(selectiveCellBoolIndex_length3_delay1Bin);
        selectiveCellIndex_suite2p_length3_delay2Bin = cellIndex_suite2p(selectiveCellBoolIndex_length3_delay2Bin);
        selectiveCellIndex_suite2p_length4_delay1Bin = cellIndex_suite2p(selectiveCellBoolIndex_length4_delay1Bin);
        selectiveCellIndex_suite2p_length4_delay2Bin = cellIndex_suite2p(selectiveCellBoolIndex_length4_delay2Bin);
        
        
        selectiveCellBoolIndex_length2 = selectiveCellBoolIndex_length2_delay1Bin | selectiveCellBoolIndex_length2_delay2Bin;
        selectiveCellIndex_suite2p_length2 = cellIndex_suite2p(selectiveCellBoolIndex_length2);
        selectiveCellNum_length2 = size(selectiveCellIndex_suite2p_length2,1);
        fprintf('selectiveCellNum_length2 = %d (%d from delay1, %d from delay2).\n',selectiveCellNum_length2,selectiveCellNum_length2_delay1Bin,selectiveCellNum_length2_delay2Bin);
        
        selectiveCellBoolIndex_length3 = selectiveCellBoolIndex_length3_delay1Bin | selectiveCellBoolIndex_length3_delay2Bin;
        selectiveCellIndex_suite2p_length3 = cellIndex_suite2p(selectiveCellBoolIndex_length3);
        selectiveCellNum_length3 = size(selectiveCellIndex_suite2p_length3,1);
        fprintf('selectiveCellNum_length3 = %d (%d from delay1, %d from delay2).\n',selectiveCellNum_length3,selectiveCellNum_length3_delay1Bin,selectiveCellNum_length3_delay2Bin);
        
        selectiveCellBoolIndex_length4 = selectiveCellBoolIndex_length4_delay1Bin | selectiveCellBoolIndex_length4_delay2Bin;
        selectiveCellIndex_suite2p_length4 = cellIndex_suite2p(selectiveCellBoolIndex_length4);
        selectiveCellNum_length4 = size(selectiveCellIndex_suite2p_length4,1);
        fprintf('selectiveCellNum_length4 = %d (%d from delay1, %d from delay2).\n',selectiveCellNum_length4,selectiveCellNum_length4_delay1Bin,selectiveCellNum_length4_delay2Bin);
        
        
        selectiveCellBoolIndex = ...
            selectiveCellBoolIndex_length1_delay1Bin | selectiveCellBoolIndex_length1_delay2Bin | ...
            selectiveCellBoolIndex_length2_delay1Bin | selectiveCellBoolIndex_length2_delay2Bin | ...
            selectiveCellBoolIndex_length3_delay1Bin | selectiveCellBoolIndex_length3_delay2Bin | ...
            selectiveCellBoolIndex_length4_delay1Bin | selectiveCellBoolIndex_length4_delay2Bin;
        selectiveCellIndex_suite2p = cellIndex_suite2p(selectiveCellBoolIndex);
        selectiveCellNum = size(selectiveCellIndex_suite2p,1);
        selectiveCellBoolIndex_delay1Bin = ...
            selectiveCellBoolIndex_length1_delay1Bin | ...
            selectiveCellBoolIndex_length2_delay1Bin | ...
            selectiveCellBoolIndex_length3_delay1Bin | ...            
            selectiveCellBoolIndex_length4_delay1Bin;
        selectiveCellBoolIndex_delay2Bin = ...
            selectiveCellBoolIndex_length1_delay2Bin | ...
            selectiveCellBoolIndex_length2_delay2Bin | ...
            selectiveCellBoolIndex_length3_delay2Bin | ...            
            selectiveCellBoolIndex_length4_delay2Bin;
        fprintf('selectiveCellNum = %d (%d from delay1, %d from delay2).\n',selectiveCellNum,sum(selectiveCellBoolIndex_delay1Bin),sum(selectiveCellBoolIndex_delay2Bin));
        
        selectiveCellBoolIndex_length123 = selectiveCellBoolIndex_length1 & selectiveCellBoolIndex_length2 & selectiveCellBoolIndex_length3;
        fprintf('%d cell tuning to length 1 2 3.\n',sum(selectiveCellBoolIndex_length123));
        selectiveCellIndex_length123_suite2p = cellIndex_suite2p(selectiveCellBoolIndex_length123);
        
        selectiveCellBoolIndex_onlyLength23 = ~selectiveCellBoolIndex_length1 & selectiveCellBoolIndex_length2 & selectiveCellBoolIndex_length3;
        fprintf('%d cell tuning to only length 2 3.\n',sum(selectiveCellBoolIndex_onlyLength23));
        selectiveCellIndex_onlyLength23_suite2p = cellIndex_suite2p(selectiveCellBoolIndex_onlyLength23);
        
        selectiveCellBoolIndex_only3 = ~selectiveCellBoolIndex_length1 & ~selectiveCellBoolIndex_length2 & selectiveCellBoolIndex_length3;
        fprintf('%d cell tuning to only length 3.\n',sum(selectiveCellBoolIndex_only3));
        selectiveCellIndex_onlyLength3_suite2p = cellIndex_suite2p(selectiveCellBoolIndex_only3);
        
        
        selectiveCellBoolIndex_length123_delay1Bin = selectiveCellBoolIndex_length1_delay1Bin & selectiveCellBoolIndex_length2_delay1Bin & selectiveCellBoolIndex_length3_delay1Bin;
        selectiveCellIndex_length123_delay1Bin_suite2p = cellIndex_suite2p(selectiveCellBoolIndex_length123_delay1Bin);
        
        
    end
    
    %% fun_rProb_glm
    if if_compute == 1
        if if_compute_fun_rProb_glm == 1
            fun_rProb_glm_frame_Name_v = autoGetFunName_myScripts('fun_rProb_glm_frame', [targetPATH '\functions']);
            fun_rProb_glm_frame = str2func(fun_rProb_glm_frame_Name_v);
            
            rProb_glm_options = struct;
            rProb_glm_options.offloadingProb_inOne = offloadingProb_inOne;
            rProb_glm_options.precision_inOne = precision_inOne;            
            rProb_glm_options.numSeq = numSeq;
            rProb_glm_options.target_seqSet = target_seqSet;
            rProb_glm_options.trial_para = trial_para;
                
            t0_rProb_glm = tic;
            for temptempi=1:3
                if temptempi == 1
                    temp_F_dff = F_dff_baseline;
                elseif temptempi == 2
                    temp_F_dff = F_dff_delay1;
                elseif temptempi == 3
                    temp_F_dff = F_dff_delay2;
                end                
                rProb_glm_options.temp_F_dff = temp_F_dff;
                rProb_glm_output = fun_rProb_glm_frame(rProb_glm_options);
                
                rProb_glm_output.selectiveCellBoolIndex_rProb_glm;
                temp_sum = sum(rProb_glm_output.selectiveCellBoolIndex_rProb_glm);
                
                if temptempi == 1
                    rProb_glm_output_baselineBin = rProb_glm_output;
                elseif temptempi == 2
                    rProb_glm_output_delay1Bin = rProb_glm_output;
                elseif temptempi == 3
                    rProb_glm_output_delay2Bin = rProb_glm_output;
                end                                
            end
            fprintf('selectiveCellNum_rProb_glm = [%d,%d,%d]/%d, t=%.1f secs.\n',...
                sum(rProb_glm_output_baselineBin.selectiveCellBoolIndex_rProb_glm),...
                sum(rProb_glm_output_delay1Bin.selectiveCellBoolIndex_rProb_glm),...
                sum(rProb_glm_output_delay2Bin.selectiveCellBoolIndex_rProb_glm),...
                size(temp_F_dff,1),toc(t0_rProb_glm));
            
            fprintf('selectiveCellNum_precision_glm = [%d,%d,%d]/%d, t=%.1f secs.\n',...
                sum(rProb_glm_output_baselineBin.selectiveCellBoolIndex_precision_glm),...
                sum(rProb_glm_output_delay1Bin.selectiveCellBoolIndex_precision_glm),...
                sum(rProb_glm_output_delay2Bin.selectiveCellBoolIndex_precision_glm),...
                size(temp_F_dff,1),toc(t0_rProb_glm));
            
            
            rProb_glm_output_all = struct;
            rProb_glm_output_all.baseline = rProb_glm_output_baselineBin;
            rProb_glm_output_all.delay1 = rProb_glm_output_delay1Bin;
            rProb_glm_output_all.delay2 = rProb_glm_output_delay2Bin;   
          
        end
    end
    
    selectiveCellNum_lengthx = 0;
    %% plot_lengthFlag
    if plot_lengthFlag == 1
        selectiveCellNum_lengthx = selectiveCellNum_length1;
        if if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 0
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length1_T;
        elseif if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 1
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length1_delay1Bin;
        elseif if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 2
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length1_delay2Bin;            
        end
        
        F_dff_lengthx_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memoryCorrect(1,:));
        F_dff_lengthx_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memoryCorrect(1,:));
        F_dff_lengthx_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memoryCorrect(1,:));
        sequence_lengthx_memoryCorrect_onehot = sequence_length1_memoryCorrect_onehot;
        
        F_dff_lengthx_period1 = F_dff_length1_period1;
        F_dff_lengthx_period2A = F_dff_length1_period2A;        
        F_dff_lengthx_period2 = F_dff_length1_period2;
        F_dff_lengthx_period3 = F_dff_length1_period3;
        F_dff_lengthx_period1_raw = F_dff_length1_period1_raw;
        F_dff_lengthx_period2A_raw = F_dff_length1_period2A_raw;        
        F_dff_lengthx_period2_raw = F_dff_length1_period2_raw;
        trialIndex_lengthx = trialIndex_length1;
        trialIndex_lengthx_seq = trialIndex_length1_seq;
        trialIndex_lengthx_location = trialIndex_length1_location;
        trialIndex_lengthx_memoryCorrect = trialIndex_length1_memoryCorrect;
        trialIndex_lengthx_choiceMemoryCorrect = trialIndex_length1_choiceMemoryCorrect;
        trialIndex_lengthx_noChoiceCorrect = trialIndex_length1_noChoiceCorrect;
        trialIndex_lengthx_choiceMemory = trialIndex_length1_choiceMemory;
        trialIndex_lengthx_choiceOffload = trialIndex_length1_choiceOffload;
        trialIndex_lengthx_noChoice = trialIndex_length1_noChoice;
        lengthx_period1_interval = length1_period1_interval;
        lengthx_period2A_interval = length1_period2A_interval;        
        lengthx_period2_interval = length1_period2_interval;
        lengthx_period3_interval = length1_period3_interval;
    elseif plot_lengthFlag == 2
        selectiveCellNum_lengthx = selectiveCellNum_length2;
        if if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 0
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length2_T;
        elseif if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 1
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length2_delay1Bin;
        elseif if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 2
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length2_delay2Bin;                        
        end
        
        F_dff_lengthx_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memoryCorrect(2,:));
        F_dff_lengthx_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memoryCorrect(2,:));
        F_dff_lengthx_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memoryCorrect(2,:));
        sequence_lengthx_memoryCorrect_onehot = sequence_length2_memoryCorrect_onehot;
        
        F_dff_lengthx_period1 = F_dff_length2_period1;
        F_dff_lengthx_period2A = F_dff_length2_period2A;        
        F_dff_lengthx_period2 = F_dff_length2_period2;
        F_dff_lengthx_period3 = F_dff_length2_period3;
        F_dff_lengthx_period1_raw = F_dff_length2_period1_raw;
        F_dff_lengthx_period2A_raw = F_dff_length2_period2A_raw;        
        F_dff_lengthx_period2_raw = F_dff_length2_period2_raw;
        trialIndex_lengthx = trialIndex_length2;
        trialIndex_lengthx_seq = trialIndex_length2_seq;
        trialIndex_lengthx_location = trialIndex_length2_location;
        trialIndex_lengthx_memoryCorrect = trialIndex_length2_memoryCorrect;
        trialIndex_lengthx_choiceMemoryCorrect = trialIndex_length2_choiceMemoryCorrect;
        trialIndex_lengthx_noChoiceCorrect = trialIndex_length2_noChoiceCorrect;
        trialIndex_lengthx_choiceMemory = trialIndex_length2_choiceMemory;
        trialIndex_lengthx_choiceOffload = trialIndex_length2_choiceOffload;
        trialIndex_lengthx_noChoice = trialIndex_length2_noChoice;
        lengthx_period1_interval = length2_period1_interval;
        lengthx_period2A_interval = length2_period2A_interval;                
        lengthx_period2_interval = length2_period2_interval;
        lengthx_period3_interval = length3_period3_interval;
    elseif plot_lengthFlag == 3
        selectiveCellNum_lengthx = selectiveCellNum_length3;
        if if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 0
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length3_T;
        elseif if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 1
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length3_delay1Bin;
        elseif if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 2
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length3_delay2Bin;            
        end
        
        F_dff_lengthx_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memoryCorrect(3,:));
        F_dff_lengthx_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memoryCorrect(3,:));
        F_dff_lengthx_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memoryCorrect(3,:));
        sequence_lengthx_memoryCorrect_onehot = sequence_length3_memoryCorrect_onehot;
        
        F_dff_lengthx_period1 = F_dff_length3_period1;
        F_dff_lengthx_period2A = F_dff_length3_period2A;        
        F_dff_lengthx_period2 = F_dff_length3_period2;
        F_dff_lengthx_period3 = F_dff_length3_period3;
        F_dff_lengthx_period1_raw = F_dff_length3_period1_raw;
        F_dff_lengthx_period2A_raw = F_dff_length3_period2A_raw;        
        F_dff_lengthx_period2_raw = F_dff_length3_period2_raw;
        trialIndex_lengthx = trialIndex_length3;
        trialIndex_lengthx_seq = trialIndex_length3_seq;
        trialIndex_lengthx_location = trialIndex_length3_location;
        trialIndex_lengthx_memoryCorrect = trialIndex_length3_memoryCorrect;
        trialIndex_lengthx_choiceMemoryCorrect = trialIndex_length3_choiceMemoryCorrect;
        trialIndex_lengthx_noChoiceCorrect = trialIndex_length3_noChoiceCorrect;
        trialIndex_lengthx_choiceMemory = trialIndex_length3_choiceMemory;
        trialIndex_lengthx_choiceOffload = trialIndex_length3_choiceOffload;
        trialIndex_lengthx_noChoice = trialIndex_length3_noChoice;
        lengthx_period1_interval = length3_period1_interval;
        lengthx_period2A_interval = length3_period2A_interval;                
        lengthx_period2_interval = length3_period2_interval;
        lengthx_period3_interval = length3_period3_interval;
    elseif plot_lengthFlag == 4
        selectiveCellNum_lengthx = selectiveCellNum_length4;
        if if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 0
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length4_T;
        elseif if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 1
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length4_delay1Bin;
        elseif if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 2
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length4_delay2Bin;            
        end
        
        F_dff_lengthx_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memoryCorrect(4,:));
        F_dff_lengthx_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memoryCorrect(4,:));
        F_dff_lengthx_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memoryCorrect(4,:));
        sequence_lengthx_memoryCorrect_onehot = sequence_length4_memoryCorrect_onehot;
        
        F_dff_lengthx_period1 = F_dff_length4_period1;
        F_dff_lengthx_period2A = F_dff_length4_period2A;        
        F_dff_lengthx_period2 = F_dff_length4_period2;
        F_dff_lengthx_period3 = F_dff_length4_period3;
        F_dff_lengthx_period1_raw = F_dff_length4_period1_raw;
        F_dff_lengthx_period2A_raw = F_dff_length4_period2A_raw;        
        F_dff_lengthx_period2_raw = F_dff_length4_period2_raw;
        trialIndex_lengthx = trialIndex_length4;
        trialIndex_lengthx_seq = trialIndex_length4_seq;
        trialIndex_lengthx_location = trialIndex_length4_location;
        trialIndex_lengthx_memoryCorrect = trialIndex_length4_memoryCorrect;
        trialIndex_lengthx_choiceMemoryCorrect = trialIndex_length4_choiceMemoryCorrect;
        trialIndex_lengthx_noChoiceCorrect = trialIndex_length4_noChoiceCorrect;
        trialIndex_lengthx_choiceMemory = trialIndex_length4_choiceMemory;
        trialIndex_lengthx_choiceOffload = trialIndex_length4_choiceOffload;
        trialIndex_lengthx_noChoice = trialIndex_length4_noChoice;
        lengthx_period1_interval = length4_period1_interval;
        lengthx_period2A_interval = length4_period2A_interval;                
        lengthx_period2_interval = length4_period2_interval;
        lengthx_period3_interval = length4_period3_interval;        
    end
    
    
    % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_onlyLength3_suite2p;
    % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_postDecisionPeriod;
    % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_decisionBin_andMON;
    %selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_decisionBin_orMON;
    % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_decisionBin_length1Seq_andMON;
    % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_decisionBin_length1Seq_orMON;
    % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_decisionBin_length1Seq_seqMerged_MN;
    % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length1_and;
    % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length1_delay2Bin;
    % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length1_delay1Bin;
    % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length1_T;
    %selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p(selectiveCellBoolIndex_baselineBin_MO&selectiveCellBoolIndex_decisionBin_MO);
    %selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_length123_delay1Bin_suite2p;
    
    if if_plotSelectiveCell_0T_1delay1_2delay2_3MNO == 3
        if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 0 ...
                || if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 1
            %selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_baselineBin_MO;
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_decisionBin_orMON;
        elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 2
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_baselineBin_afterMemoryCorrectError;
        elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 3
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_baselineBin_afterMemoryCorrect_MO;
        elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 4
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_baselineBin_afterMemoryError_MO;
        end
    end
    
    
    if if_plot_rProb_glm == 1
        
        rProb_glm_output = rProb_glm_output_delay1Bin;
        %rProb_glm_output = rProb_glm_output_delay2Bin;
        
        selectiveCellBoolIndex_rProb_glm = rProb_glm_output.selectiveCellBoolIndex_rProb_glm;  
        
        %selectiveCellBoolIndex_rProb_glm = selectiveCellBoolIndex_rProb_glm & tempSelectiveCellBoolIndex_length1_currentFOV;
        
        %selectiveCellBoolIndex_rProb_glm = selectiveCellBoolIndex_rProb_glm & (~tempSelectiveCellBoolIndex_length1_currentFOV);
        
        r2_rProb_glm = rProb_glm_output.r2_rProb_glm;
        
        %r2_rProb_glm(r2_rProb_glm<0.1) = -1;
        
        [r2_rProb_glm_valid_resorted,I] = sort(r2_rProb_glm(selectiveCellBoolIndex_rProb_glm),'descend');
        temp1 = r2_rProb_glm_valid_resorted > 0.1;
        I = I(temp1);
        
        %tempIndex = selectiveCellBoolIndex_rProb_glm;
        
        selectiveCellBoolIndex_rProb_glm_suite2p = cellIndex_suite2p(selectiveCellBoolIndex_rProb_glm);
        
        selectiveCellBoolIndex_rProb_glm_suite2p_resorted = selectiveCellBoolIndex_rProb_glm_suite2p(I);
        %selectiveCellBoolIndex_rProb_glm_suite2p_resorted(1:20)'

        selectiveCellIndex_suite2p_lengthx = selectiveCellBoolIndex_rProb_glm_suite2p_resorted;                
        %selectiveCellIndex_suite2p_lengthx = selectiveCellBoolIndex_rProb_glm_suite2p;    
                
    end
    
    if if_plotSelectiveCell == 1
        %plotRoiNum = min(plotRoiNum,selectiveCellNum);
        %plotRoiNum = min(plotRoiNum,selectiveCellNum_lengthx);
        %plotRoiNum = min(plotRoiNum,length(selectiveCellIndex_suite2p_lengthx));
        plotRoiNum = min(plotRoiNum,sum(selectiveCellIndex_suite2p_lengthx>=temp_id_suite2p));
        
        selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_lengthx(selectiveCellIndex_suite2p_lengthx>=temp_id_suite2p);
    end
    
    %if false
    %if true
    if if_plot_svm_beta == 1
        plotRoiNum = numFrames*svm_beta.topX;
        if plot_lengthFlag == 1
            beta_loc_lengthx_topXIndex = svm_beta.beta_loc_length1_topXIndex;
        elseif plot_lengthFlag == 2
            beta_loc_lengthx_topXIndex = svm_beta.beta_loc_length2_topXIndex;
        elseif plot_lengthFlag == 3
            beta_loc_lengthx_topXIndex = svm_beta.beta_loc_length3_topXIndex;
        end
        tempIndex = reshape(beta_loc_lengthx_topXIndex,[],1);
        selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p(tempIndex);
    end
    
    if if_plot_memoryMetaMismatch == 1
        memoryMetaMismatch;
        cellIndex_suite2p_memoryMetaMismatch = memoryMetaMismatch.cellIndex_suite2p_memoryMetaMismatch;
        cellIndex_suite2p_memoryMetaMatch = memoryMetaMismatch.cellIndex_suite2p_memoryMetaMatch;
        
        %cellIndex_suite2p_temptemp = cellIndex_suite2p_memoryMetaMismatch;
        %cellIndex_suite2p_temptemp = cellIndex_suite2p_temptemp(1);
        
        cellIndex_suite2p_temptemp = cellIndex_suite2p_memoryMetaMatch;
        
        cellIndex_suite2p_temptemp = 20;
        
        plotRoiNum = length(cellIndex_suite2p_temptemp);
        selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p_temptemp;
    end
    
    fprintf('plotCellNum = %d.\n',plotRoiNum);
    
    
    decodingData.selectiveCellBoolIndex = selectiveCellBoolIndex;
    % decodingData.selectiveCellBoolIndex = selectiveCellBoolIndex_length1;
    decodingData.F_dff_decisionPeriodA = F_dff_decisionPeriodA;    
    decodingData.decisionPeriodA_interval = decisionPeriodA_interval;    
    decodingData.F_dff_decisionPeriod = F_dff_decisionPeriod;
    decodingData.decisionPeriod_interval = decisionPeriod_interval;
    decodingData.seqIndex = seqIndex;
    decodingData.trial_para_currentSequence = trial_para.currentSequence;
    decodingData.trialIndex_bool_memoryCorrect = trialIndex_bool_memoryCorrect;
    decodingData.target_seqSet_inOne = target_seqSet_inOne;
    decodingData.target_seqSet = target_seqSet;
    
    if if_plotPCA == 1
        
        %% Compure beta from GLM, memory correct, delay1Bin
        
        a = 1;
        
        %order_glm = 1;% 0 means whole seq
        
        %plot_lengthFlag = 1;
        %order_glm = 1;
        
        order_glm_valid = min(order_glm,plot_lengthFlag);
        sequence_lengthx_memoryCorrect_onehot_oneOrder = sequence_lengthx_memoryCorrect_onehot;
        if order_glm_valid > 0
            for tempi=1:size(sequence_lengthx_memoryCorrect_onehot,1)
                a1 = find(sequence_lengthx_memoryCorrect_onehot(tempi,:)==true,order_glm_valid);
                a1 = a1(end);
                a2 = false(1,numFrames);
                a2(a1) = true;
                sequence_lengthx_memoryCorrect_onehot_oneOrder(tempi,:) = a2;
            end
        end
        
        sequence_lengthx_memoryCorrect_onehot_order = zeros(size(sequence_lengthx_memoryCorrect_onehot,1),plot_lengthFlag*numFrames);
        for tempOrder=1:plot_lengthFlag
            temp_order_glm_valid = min(tempOrder,plot_lengthFlag);
            temp_onehot = sequence_lengthx_memoryCorrect_onehot;
            if temp_order_glm_valid > 0
                for tempi=1:size(sequence_lengthx_memoryCorrect_onehot,1)
                    a1 = find(sequence_lengthx_memoryCorrect_onehot(tempi,:)==true,temp_order_glm_valid);
                    a1 = a1(end);
                    a2 = false(1,numFrames);
                    a2(a1) = true;
                    temp_onehot(tempi,:) = a2;
                end
            end
            temp_range = (1:numFrames) + numFrames*(tempOrder-1);
            sequence_lengthx_memoryCorrect_onehot_order(:,temp_range) = temp_onehot;
        end
        
        sequence_lengthx_memoryCorrect_onehot_oneOrder_raw = sequence_lengthx_memoryCorrect_onehot_oneOrder;
        
        temp_locValidBool = true(1,numFrames);
        if order_glm_valid > 0
            temp_locValidBool(numFrames-plot_lengthFlag+order_glm_valid+1:end) = false;
            temp_locValidBool(1:(order_glm_valid-1)) = false;
        end
        
        temp_locValidBool_real = true(1,numFrames);
        for tempi=1:numFrames
            if sum(sequence_lengthx_memoryCorrect_onehot_oneOrder_raw(:,tempi)) == 0
                temp_locValidBool_real(tempi) = false;
            end
        end
        
        
        %temp_locValidBool(numFrames-order_glm_valid+1)
        
        %sequence_lengthx_memoryCorrect_onehot_oneOrder = sequence_lengthx_memoryCorrect_onehot_oneOrder_raw(:,temp_locValidBool);
        sequence_lengthx_memoryCorrect_onehot_oneOrder = sequence_lengthx_memoryCorrect_onehot_oneOrder_raw;
        
        glm_mdl = cell(size(F_dff_lengthx_delay1Bin,1),1);
        if order_glm_valid == 0
            glm_beta = zeros(size(F_dff_lengthx_delay1Bin,1),numFrames);
            glm_pValue = zeros(size(F_dff_lengthx_delay1Bin,1),numFrames);
        else
            glm_beta = zeros(size(F_dff_lengthx_delay1Bin,1),numFrames*plot_lengthFlag);
            glm_pValue = zeros(size(F_dff_lengthx_delay1Bin,1),numFrames*plot_lengthFlag);
        end
        glm_r2 = zeros(size(F_dff_lengthx_delay1Bin,1),1);
        parfor tempi=1:size(F_dff_lengthx_delay1Bin,1)
            %for tempi=1:size(F_dff_lengthx_delay1Bin,1)
            warning('off');
            if order_glm_valid == 0
                x = sequence_lengthx_memoryCorrect_onehot;
            else
                x = sequence_lengthx_memoryCorrect_onehot_order;
            end
            y0 = F_dff_lengthx_delay1Bin(tempi,:)';
            y = (y0-mean(y0))/std(y0); % z-score
            % glm_mdl = fitglm(x,y,'linear','Distribution','normal');
            temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',false);
            %temp_mdl= fitglm(x,y,'linear','Distribution','poisson','Intercept',false);
            glm_mdl{tempi} = temp_mdl;
            glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
            glm_pValue(tempi,:) = temp_mdl.Coefficients.pValue;
            %glm_beta(tempi,temp_locValidBool) = temp_mdl.Coefficients.Estimate;
            %glm_pValue(tempi,temp_locValidBool) = temp_mdl.Coefficients.pValue;
            glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
            warning('on');
        end
        glm_beta;
        %glm_beta_memoryCorrect_lengthx_delay1Bin = glm_beta(:,temp_locValidBool);
        if order_glm_valid == 0
            temp_range = 1:numFrames;
        else
            temp_range = (1:numFrames) + numFrames*(order_glm_valid-1);
        end
        glm_beta_memoryCorrect_lengthx_delay1Bin = glm_beta(:,temp_range);
        
        %% Compure beta from GLM, memory correct, delay2Bin
        glm_mdl = cell(size(F_dff_lengthx_delay2Bin,1),1);
        if order_glm_valid == 0
            glm_beta = zeros(size(F_dff_lengthx_delay2Bin,1),numFrames);
            glm_pValue = zeros(size(F_dff_lengthx_delay2Bin,1),numFrames);
        else
            glm_beta = zeros(size(F_dff_lengthx_delay2Bin,1),numFrames*plot_lengthFlag);
            glm_pValue = zeros(size(F_dff_lengthx_delay2Bin,1),numFrames*plot_lengthFlag);
        end
        glm_r2 = zeros(size(F_dff_lengthx_delay2Bin,1),1);
        parfor tempi=1:size(F_dff_lengthx_delay2Bin,1)
            warning('off');
            if order_glm_valid == 0
                x = sequence_lengthx_memoryCorrect_onehot;
            else
                x = sequence_lengthx_memoryCorrect_onehot_order;
            end
            y0 = F_dff_lengthx_delay2Bin(tempi,:)';
            y = (y0-mean(y0))/std(y0); % z-score
            temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',false);
            glm_mdl{tempi} = temp_mdl;
            %glm_beta(tempi,temp_locValidBool) = temp_mdl.Coefficients.Estimate;
            %glm_pValue(tempi,temp_locValidBool) = temp_mdl.Coefficients.pValue;
            glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
            glm_pValue(tempi,:) = temp_mdl.Coefficients.pValue;
            glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
            warning('on');
        end
        %glm_beta_memoryCorrect_lengthx_delay2Bin = glm_beta(:,temp_locValidBool);
        if order_glm_valid == 0
            temp_range = 1:numFrames;
        else
            temp_range = (1:numFrames) + numFrames*(order_glm_valid-1);
        end
        glm_beta_memoryCorrect_lengthx_delay2Bin = glm_beta(:,temp_range);
        
        %% Compure beta from GLM, memory correct, baselineBin
        glm_mdl = cell(size(F_dff_lengthx_baselineBin,1),1);
        if order_glm_valid == 0
            glm_beta = zeros(size(F_dff_lengthx_baselineBin,1),numFrames);
            glm_pValue = zeros(size(F_dff_lengthx_baselineBin,1),numFrames);
        else
            glm_beta = zeros(size(F_dff_lengthx_baselineBin,1),numFrames*plot_lengthFlag);
            glm_pValue = zeros(size(F_dff_lengthx_baselineBin,1),numFrames*plot_lengthFlag);
        end
        glm_r2 = zeros(size(F_dff_lengthx_baselineBin,1),1);
        parfor tempi=1:size(F_dff_lengthx_baselineBin,1)
            warning('off');
            if order_glm_valid == 0
                x = sequence_lengthx_memoryCorrect_onehot;
            else
                x = sequence_lengthx_memoryCorrect_onehot_order;
            end
            y0 = F_dff_lengthx_baselineBin(tempi,:)';
            y = (y0-mean(y0))/std(y0); % z-score
            temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',false);
            glm_mdl{tempi} = temp_mdl;
            %glm_beta(tempi,temp_locValidBool) = temp_mdl.Coefficients.Estimate;
            %glm_pValue(tempi,temp_locValidBool) = temp_mdl.Coefficients.pValue;
            glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
            glm_pValue(tempi,:) = temp_mdl.Coefficients.pValue;
            glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
            warning('on');
        end
        %glm_beta_memoryCorrect_lengthx_baselineBin = glm_beta(:,temp_locValidBool);
        if order_glm_valid == 0
            temp_range = 1:numFrames;
        else
            temp_range = (1:numFrames) + numFrames*(order_glm_valid-1);
        end
        glm_beta_memoryCorrect_lengthx_baselineBin = glm_beta(:,temp_range);
        
        a = 1;
        
        %     %% Compure beta from GLM, choice memory, delay1Bin
        %     sequence_length1_choiceMemory = sequence_length1(choiceMemoryIndex_trialIndex_length1);
        %     sequence_length1_choiceMemory_onehot = false(size(sequence_length1_choiceMemory,2),numFrames);
        %     for tempi=1:size(sequence_length1_choiceMemory,2)
        %         sequence_length1_choiceMemory_onehot(tempi,sequence_length1_choiceMemory(tempi)) = true;
        %     end
        %
        %     F_dff_choiceMemory_length1_delay1Bin = mean(F_dff_choiceMemory_length1_period1(:,:,length1_period1_interval(3):length1_period1_interval(end)),3);
        %
        %     glm_mdl = cell(size(F_dff_choiceMemory_length1_delay1Bin,1),1);
        %     glm_beta = zeros(size(F_dff_choiceMemory_length1_delay1Bin,1),numFrames);
        %     glm_pValue = zeros(size(F_dff_choiceMemory_length1_delay1Bin,1),numFrames);
        %     glm_r2 = zeros(size(F_dff_choiceMemory_length1_delay1Bin,1),1);
        %
        %     parfor tempi=1:size(F_dff_choiceMemory_length1_delay1Bin,1)
        %         warning('off');
        %         x = sequence_length1_choiceMemory_onehot;
        %         y0 = F_dff_choiceMemory_length1_delay1Bin(tempi,:)';
        %         y = (y0-mean(y0))/std(y0); % z-score
        %         temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',false);
        %         glm_mdl{tempi} = temp_mdl;
        %         glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
        %         glm_pValue(tempi,:) = temp_mdl.Coefficients.pValue;
        %         glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
        %         warning('on');
        %     end
        %     glm_beta_choiceMemory_length1_delay1Bin = glm_beta;
        %
        %     %% Compure beta from GLM, choice memory, delay2Bin
        %     F_dff_choiceMemory_length1_delay2Bin = mean(F_dff_choiceMemory_length1_period2(:,:,length1_period2_interval(2):length1_period2_interval(end)),3);
        %
        %     glm_mdl = cell(size(F_dff_choiceMemory_length1_delay2Bin,1),1);
        %     glm_beta = zeros(size(F_dff_choiceMemory_length1_delay2Bin,1),numFrames);
        %     glm_pValue = zeros(size(F_dff_choiceMemory_length1_delay2Bin,1),numFrames);
        %     glm_r2 = zeros(size(F_dff_choiceMemory_length1_delay2Bin,1),1);
        %
        %     parfor tempi=1:size(F_dff_choiceMemory_length1_delay2Bin,1)
        %         warning('off');
        %         x = sequence_length1_choiceMemory_onehot;
        %         y0 = F_dff_choiceMemory_length1_delay2Bin(tempi,:)';
        %         y = (y0-mean(y0))/std(y0); % z-score
        %         temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',false);
        %         glm_mdl{tempi} = temp_mdl;
        %         glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
        %         glm_pValue(tempi,:) = temp_mdl.Coefficients.pValue;
        %         glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
        %         warning('on');
        %     end
        %     glm_beta_choiceMemory_length1_delay2Bin = glm_beta;
        %
        %     %% Compure beta from GLM, choice memory correct, delay1Bin
        %     sequence_length1_choiceMemoryCorrect = sequence_length1(choiceMemoryCorrectIndex_trialIndex_length1);
        %     sequence_length1_choiceMemoryCorrect_onehot = false(size(sequence_length1_choiceMemoryCorrect,2),numFrames);
        %     for tempi=1:size(sequence_length1_choiceMemoryCorrect,2)
        %         sequence_length1_choiceMemoryCorrect_onehot(tempi,sequence_length1_choiceMemoryCorrect(tempi)) = true;
        %     end
        %
        %     F_dff_choiceMemoryCorrect_length1_delay1Bin = mean(F_dff_choiceMemoryCorrect_length1_period1(:,:,length1_period1_interval(3):length1_period1_interval(end)),3);
        %
        %     glm_mdl = cell(size(F_dff_choiceMemoryCorrect_length1_delay1Bin,1),1);
        %     glm_beta = zeros(size(F_dff_choiceMemoryCorrect_length1_delay1Bin,1),numFrames);
        %     glm_pValue = zeros(size(F_dff_choiceMemoryCorrect_length1_delay1Bin,1),numFrames);
        %     glm_r2 = zeros(size(F_dff_choiceMemoryCorrect_length1_delay1Bin,1),1);
        %
        %     parfor tempi=1:size(F_dff_choiceMemoryCorrect_length1_delay1Bin,1)
        %         warning('off');
        %         x = sequence_length1_choiceMemoryCorrect_onehot;
        %         y0 = F_dff_choiceMemoryCorrect_length1_delay1Bin(tempi,:)';
        %         y = (y0-mean(y0))/std(y0); % z-score
        %         temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',false);
        %         glm_mdl{tempi} = temp_mdl;
        %         glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
        %         glm_pValue(tempi,:) = temp_mdl.Coefficients.pValue;
        %         glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
        %         warning('on');
        %     end
        %     glm_beta_choiceMemoryCorrect_length1_delay1Bin = glm_beta;
        %
        %
        %     %% Compure beta from GLM, choice memory correct, delay2Bin
        %     F_dff_choiceMemoryCorrect_length1_delay2Bin = mean(F_dff_choiceMemoryCorrect_length1_period2(:,:,length1_period2_interval(2):length1_period2_interval(end)),3);
        %
        %     glm_mdl = cell(size(F_dff_choiceMemoryCorrect_length1_delay2Bin,1),1);
        %     glm_beta = zeros(size(F_dff_choiceMemoryCorrect_length1_delay2Bin,1),numFrames);
        %     glm_pValue = zeros(size(F_dff_choiceMemoryCorrect_length1_delay2Bin,1),numFrames);
        %     glm_r2 = zeros(size(F_dff_choiceMemoryCorrect_length1_delay2Bin,1),1);
        %
        %     parfor tempi=1:size(F_dff_choiceMemoryCorrect_length1_delay2Bin,1)
        %         warning('off');
        %         x = sequence_length1_choiceMemoryCorrect_onehot;
        %         y0 = F_dff_choiceMemoryCorrect_length1_delay2Bin(tempi,:)';
        %         y = (y0-mean(y0))/std(y0); % z-score
        %         temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',false);
        %         glm_mdl{tempi} = temp_mdl;
        %         glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
        %         glm_pValue(tempi,:) = temp_mdl.Coefficients.pValue;
        %         glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
        %         warning('on');
        %     end
        %     glm_beta_choiceMemoryCorrect_length1_delay2Bin = glm_beta;
        %
        %     %% Compure beta from GLM, offload, delay1Bin
        %     sequence_length1_offload = sequence_length1(offloadIndex_trialIndex_length1);
        %     sequence_length1_offload_onehot = false(size(sequence_length1_offload,2),numFrames);
        %     for tempi=1:size(sequence_length1_offload,2)
        %         sequence_length1_offload_onehot(tempi,sequence_length1_offload(tempi)) = true;
        %     end
        %
        %     F_dff_offload_length1_delay1Bin = mean(F_dff_offload_length1_period1(:,:,length1_period1_interval(3):length1_period1_interval(end)),3);
        %
        %     glm_mdl = cell(size(F_dff_offload_length1_delay1Bin,1),1);
        %     glm_beta = zeros(size(F_dff_offload_length1_delay1Bin,1),numFrames);
        %     glm_pValue = zeros(size(F_dff_offload_length1_delay1Bin,1),numFrames);
        %     glm_r2 = zeros(size(F_dff_offload_length1_delay1Bin,1),1);
        %
        %     parfor tempi=1:size(F_dff_offload_length1_delay1Bin,1)
        %         warning('off');
        %         x = sequence_length1_offload_onehot;
        %         y0 = F_dff_offload_length1_delay1Bin(tempi,:)';
        %         y = (y0-mean(y0))/std(y0); % z-score
        %         temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',false);
        %         glm_mdl{tempi} = temp_mdl;
        %         glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
        %         glm_pValue(tempi,:) = temp_mdl.Coefficients.pValue;
        %         glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
        %         warning('on');
        %     end
        %     glm_beta_offload_length1_delay1Bin = glm_beta;
        %
        %
        %     %% Compure beta from GLM, offload, delay2Bin
        %     F_dff_offload_length1_delay2Bin = mean(F_dff_offload_length1_period2(:,:,length1_period2_interval(2):length1_period2_interval(end)),3);
        %
        %     glm_mdl = cell(size(F_dff_offload_length1_delay2Bin,1),1);
        %     glm_beta = zeros(size(F_dff_offload_length1_delay2Bin,1),numFrames);
        %     glm_pValue = zeros(size(F_dff_offload_length1_delay2Bin,1),numFrames);
        %     glm_r2 = zeros(size(F_dff_offload_length1_delay2Bin,1),1);
        %
        %     parfor tempi=1:size(F_dff_offload_length1_delay2Bin,1)
        %         warning('off');
        %         x = sequence_length1_offload_onehot;
        %         y0 = F_dff_offload_length1_delay2Bin(tempi,:)';
        %         y = (y0-mean(y0))/std(y0); % z-score
        %         temp_mdl= fitglm(x,y,'linear','Distribution','normal','Intercept',false);
        %         glm_mdl{tempi} = temp_mdl;
        %         glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
        %         glm_pValue(tempi,:) = temp_mdl.Coefficients.pValue;
        %         glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
        %         warning('on');
        %     end
        %     glm_beta_offload_length1_delay2Bin = glm_beta;
        
        
        %% Compure PCA based on GLM
        x = glm_beta_memoryCorrect_lengthx_delay1Bin';
        %     x = glm_beta_memoryCorrect_lengthx_delay2Bin';
        %     x = glm_beta_memoryCorrect_lengthx_baselineBin';
        
        % x = x(:,selectiveCellBoolIndex);
        [coeff,score,latent,tsquared,explained,mu] = pca(x);
    end
    
    if if_plotPCA == 1
        close all
        
        multi_rgbColor = ...
            [228,26,28;
            55,126,184;
            77,175,74;
            152,78,163;
            255,127,0;
            255,255,51]/255;
        
        %temp_rgbColor(1,:,:) = multi_rgbColor;
        %image(temp_rgbColor);
        %set(gca,'position',[0 0 1 1]);
        %set(gca,'XTick',[0:0:0])
        %set(gca,'YTick',[0:0:0])
        
        backgounrdColor = [1 1 1]*0.825;%0.875
        
        %plot(score(:,1),score(:,2),'k','lineWidth',1.5);
        %plot(score(1:sum(temp_locValidBool),1),score(1:sum(temp_locValidBool),2),'k','lineWidth',1.5);
        
        plot(score(temp_locValidBool,1),score(temp_locValidBool,2),'k','lineWidth',1.5);
        %plot3(score(temp_locValidBool,1),score(temp_locValidBool,2),score(temp_locValidBool,3),'k','lineWidth',1.5);
        hold on
        
        h_plot = [];
        tempj = 0;
        for tempi=1:numFrames
            %if temp_locValidBool(tempi) == false
            %    continue
            %end
            tempj = tempj + 1;
            if temp_locValidBool_real(tempi) == false
                %temph = scatter(score(tempj,1),score(tempj,2),120,multi_rgbColor(tempi,:));
            elseif temp_locValidBool_real(tempi) == true
                temph = scatter(score(tempj,1),score(tempj,2),120,multi_rgbColor(tempi,:),'filled');%36
                %temph = scatter3(score(tempi,1),score(tempi,2),score(tempi,3),120,multi_rgbColor(tempi,:),'filled');
            end
            %scatter3(score(tempi,1),score(tempi,2),score(tempi,3),36,multi_rgbColor(tempi,:),'filled');
            hold on
            h_plot = [h_plot temph];
        end
        
        
        
        %legend('location1','location2','location3','location4','location5','location6','Location','northeast','fontsize',8)    temp_labels = {'location1','location2','location3','location4','location5','location6'};
        %temp_labels = {'location1','location2','location3','location4','location5','location6'};
        temp_labels = {'1','2','3','4','5','6'};
        %legend(h_plot,temp_labels(temp_locValidBool),'Location','southeast','fontsize',20);%8
        
        set(gca,'color',backgounrdColor);
        axis equal
        xlabel(sprintf('1st Principal Component (%.1f%%)',explained(1)),'fontsize',16)
        ylabel(sprintf('2nd Principal Component (%.1f%%)',explained(2)),'fontsize',16)
        %zlabel(sprintf('3rd Principal Component (%.1f%%)',explained(3)))
    end
    
    %% Plot single cell
    if if_plotSingleCell == 1
        
        tempRoiCount = 0;
        for tempIndex=1:plotRoiNum
            %temp_id2 = find(cellIndex == temp_id);
            
            %% compute cellID_F_dff
            tempRoiCount = tempRoiCount + 1;
            tempflag = 0;
            while tempflag == 0
                if if_plotSelectiveCell == 0
                    temp_id2 = find(cellIndex == (temp_id+tempRoiCount-1));
                    if isempty(temp_id2)
                        tempRoiCount = tempRoiCount + 1;
                    else
                        tempflag = 1;
                    end
                elseif if_plotSelectiveCell == 1
                    %temp_id3 = find((selectiveCellIndex_suite2p_lengthx+1) == (temp_id+tempRoiCount-1));
                    %if isempty(temp_id3)
                    %    tempRoiCount = tempRoiCount + 1;
                    %else
                    %    tempflag = 1;
                    %    temp_id2 = find(cellIndex_suite2p == selectiveCellIndex_suite2p_lengthx(temp_id3));
                    %end
                    
                    %temp_id2 = find(cellIndex_suite2p == selectiveCellIndex_suite2p_lengthx(plotRoiNum-tempIndex+1));
                    temp_id2 = find(cellIndex_suite2p == selectiveCellIndex_suite2p_lengthx(tempIndex));
                    break
                end
                
            end
            temp_id2;
            temp_cellIndex_suite2p = cellIndex(temp_id2)-1;
            temp_cellIndex = cellIndex(temp_id2);
            
            if temp_cellIndex_suite2p == 39
                a = 1;
            end
            cellID_F_dff_lengthx_period1 = squeeze(F_dff_lengthx_period1(temp_id2,:,:));
            cellID_F_dff_lengthx_period2A = squeeze(F_dff_lengthx_period2A(temp_id2,:,:));            
            cellID_F_dff_lengthx_period2 = squeeze(F_dff_lengthx_period2(temp_id2,:,:));
            cellID_F_dff_lengthx_period3 = squeeze(F_dff_lengthx_period3(temp_id2,:,:));
            
            cellID_F_dff_lengthx_period1_location = cell(1,numFrames);
            cellID_F_dff_lengthx_period2A_location = cell(1,numFrames);            
            cellID_F_dff_lengthx_period2_location = cell(1,numFrames);
            cellID_F_dff_lengthx_period3_location = cell(1,numFrames);
            cellID_F_dff_lengthx_period1_location_mean = zeros(numFrames,size(cellID_F_dff_lengthx_period1,2));
            cellID_F_dff_lengthx_period2A_location_mean = zeros(numFrames,size(cellID_F_dff_lengthx_period2A,2));            
            cellID_F_dff_lengthx_period2_location_mean = zeros(numFrames,size(cellID_F_dff_lengthx_period2,2));
            cellID_F_dff_lengthx_period3_location_mean = zeros(numFrames,size(cellID_F_dff_lengthx_period3,2));
            cellID_F_dff_lengthx_period1_location_sem = zeros(numFrames,size(cellID_F_dff_lengthx_period1,2));
            cellID_F_dff_lengthx_period2A_location_sem = zeros(numFrames,size(cellID_F_dff_lengthx_period2A,2));            
            cellID_F_dff_lengthx_period2_location_sem = zeros(numFrames,size(cellID_F_dff_lengthx_period2,2));
            cellID_F_dff_lengthx_period3_location_sem = zeros(numFrames,size(cellID_F_dff_lengthx_period3,2));
            
            trialIndex_lengthx;
            trialIndex_lengthx_seq;
            trialIndex_lengthx_location;
            
            for tempi=1:numFrames
                %[~,Locb] = ismember(trialIndex_lengthx_seq{tempi},trialIndex_lengthx_memoryCorrect);
                [~,Locb] = ismember(trialIndex_lengthx_location{tempi},trialIndex_lengthx_memoryCorrect);
                
                Locb(Locb==0) = [];
                cellID_F_dff_lengthx_period1_location{tempi} = cellID_F_dff_lengthx_period1(Locb,:);
                cellID_F_dff_lengthx_period2A_location{tempi} = cellID_F_dff_lengthx_period2A(Locb,:);                
                cellID_F_dff_lengthx_period2_location{tempi} = cellID_F_dff_lengthx_period2(Locb,:);
                cellID_F_dff_lengthx_period3_location{tempi} = cellID_F_dff_lengthx_period3(Locb,:);
                
                cellID_F_dff_lengthx_period1_location_mean(tempi,:) = mean(cellID_F_dff_lengthx_period1_location{tempi},1);
                cellID_F_dff_lengthx_period2A_location_mean(tempi,:) = mean(cellID_F_dff_lengthx_period2A_location{tempi},1);                
                cellID_F_dff_lengthx_period2_location_mean(tempi,:) = mean(cellID_F_dff_lengthx_period2_location{tempi},1);
                cellID_F_dff_lengthx_period3_location_mean(tempi,:) = mean(cellID_F_dff_lengthx_period3_location{tempi},1);
                
                cellID_F_dff_lengthx_period1_location_sem(tempi,:) = std(cellID_F_dff_lengthx_period1_location{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period1_location{tempi},1));
                cellID_F_dff_lengthx_period2A_location_sem(tempi,:) = std(cellID_F_dff_lengthx_period2A_location{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2A_location{tempi},1));                
                cellID_F_dff_lengthx_period2_location_sem(tempi,:) = std(cellID_F_dff_lengthx_period2_location{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2_location{tempi},1));
                cellID_F_dff_lengthx_period3_location_sem(tempi,:) = std(cellID_F_dff_lengthx_period3_location{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period3_location{tempi},1));
            end
            
            a = 1;
            
            cellID_F_dff_baselinePeriod_afterMemoryCorrect = squeeze(F_dff_baselinePeriod(temp_id2,trialIndex_bool_afterMemoryCorrect,:));
            cellID_F_dff_baselinePeriod_afterMemoryError = squeeze(F_dff_baselinePeriod(temp_id2,trialIndex_bool_afterMemoryError,:));
            
            cellID_F_dff_baselinePeriod_afterMemoryCorrect_mean = mean(cellID_F_dff_baselinePeriod_afterMemoryCorrect,1);
            cellID_F_dff_baselinePeriod_afterMemoryError_mean = mean(cellID_F_dff_baselinePeriod_afterMemoryError,1);
            
            cellID_F_dff_baselinePeriod_afterMemoryCorrect_sem = std(cellID_F_dff_baselinePeriod_afterMemoryCorrect,1)...
                ./sqrt(size(cellID_F_dff_baselinePeriod_afterMemoryCorrect,1));
            cellID_F_dff_baselinePeriod_afterMemoryError_sem = std(cellID_F_dff_baselinePeriod_afterMemoryError,1)...
                ./sqrt(size(cellID_F_dff_baselinePeriod_afterMemoryError,1));
            
            cellID_F_dff_decisionPeriodA_afterMemoryCorrect = squeeze(F_dff_decisionPeriodA(temp_id2,trialIndex_bool_afterMemoryCorrect,:));
            cellID_F_dff_decisionPeriodA_afterMemoryError = squeeze(F_dff_decisionPeriodA(temp_id2,trialIndex_bool_afterMemoryError,:));
            
            cellID_F_dff_decisionPeriodA_afterMemoryCorrect_mean = mean(cellID_F_dff_decisionPeriodA_afterMemoryCorrect,1);
            cellID_F_dff_decisionPeriodA_afterMemoryError_mean = mean(cellID_F_dff_decisionPeriodA_afterMemoryError,1);
            
            cellID_F_dff_decisionPeriodA_afterMemoryCorrect_sem = std(cellID_F_dff_decisionPeriodA_afterMemoryCorrect,1)...
                ./sqrt(size(cellID_F_dff_decisionPeriodA_afterMemoryCorrect,1));
            cellID_F_dff_decisionPeriodA_afterMemoryError_sem = std(cellID_F_dff_decisionPeriodA_afterMemoryError,1)...
                ./sqrt(size(cellID_F_dff_decisionPeriodA_afterMemoryError,1));

            
            cellID_F_dff_decisionPeriod_afterMemoryCorrect = squeeze(F_dff_decisionPeriod(temp_id2,trialIndex_bool_afterMemoryCorrect,:));
            cellID_F_dff_decisionPeriod_afterMemoryError = squeeze(F_dff_decisionPeriod(temp_id2,trialIndex_bool_afterMemoryError,:));
            
            cellID_F_dff_decisionPeriod_afterMemoryCorrect_mean = mean(cellID_F_dff_decisionPeriod_afterMemoryCorrect,1);
            cellID_F_dff_decisionPeriod_afterMemoryError_mean = mean(cellID_F_dff_decisionPeriod_afterMemoryError,1);
            
            cellID_F_dff_decisionPeriod_afterMemoryCorrect_sem = std(cellID_F_dff_decisionPeriod_afterMemoryCorrect,1)...
                ./sqrt(size(cellID_F_dff_decisionPeriod_afterMemoryCorrect,1));
            cellID_F_dff_decisionPeriod_afterMemoryError_sem = std(cellID_F_dff_decisionPeriod_afterMemoryError,1)...
                ./sqrt(size(cellID_F_dff_decisionPeriod_afterMemoryError,1));
            
            
            if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 <= 2
                cellID_F_dff_baselinePeriod_choiceMemory = squeeze(F_dff_baselinePeriod(temp_id2,trialIndex_bool_choiceMemory,:));
                cellID_F_dff_baselinePeriod_choiceOffload = squeeze(F_dff_baselinePeriod(temp_id2,trialIndex_bool_choiceOffload,:));
                cellID_F_dff_baselinePeriod_noChoice = squeeze(F_dff_baselinePeriod(temp_id2,trialIndex_bool_noChoice,:));
            elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 3
                cellID_F_dff_baselinePeriod_choiceMemory = squeeze(F_dff_baselinePeriod(temp_id2,(trialIndex_bool_choiceMemory&trialIndex_bool_afterMemoryCorrect),:));
                cellID_F_dff_baselinePeriod_choiceOffload = squeeze(F_dff_baselinePeriod(temp_id2,(trialIndex_bool_choiceOffload&trialIndex_bool_afterMemoryCorrect),:));
                cellID_F_dff_baselinePeriod_noChoice = squeeze(F_dff_baselinePeriod(temp_id2,(trialIndex_bool_noChoice&trialIndex_bool_afterMemoryCorrect),:));
            elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 4
                cellID_F_dff_baselinePeriod_choiceMemory = squeeze(F_dff_baselinePeriod(temp_id2,(trialIndex_bool_choiceMemory&trialIndex_bool_afterMemoryError),:));
                cellID_F_dff_baselinePeriod_choiceOffload = squeeze(F_dff_baselinePeriod(temp_id2,(trialIndex_bool_choiceOffload&trialIndex_bool_afterMemoryError),:));
                cellID_F_dff_baselinePeriod_noChoice = squeeze(F_dff_baselinePeriod(temp_id2,(trialIndex_bool_noChoice&trialIndex_bool_afterMemoryError),:));
            end
            
            cellID_F_dff_baselinePeriod_choiceMemory_mean = mean(cellID_F_dff_baselinePeriod_choiceMemory,1);
            cellID_F_dff_baselinePeriod_choiceOffload_mean = mean(cellID_F_dff_baselinePeriod_choiceOffload,1);
            cellID_F_dff_baselinePeriod_noChoice_mean = mean(cellID_F_dff_baselinePeriod_noChoice,1);
            
            cellID_F_dff_baselinePeriod_choiceMemory_sem = std(cellID_F_dff_baselinePeriod_choiceMemory,1)...
                ./sqrt(size(cellID_F_dff_baselinePeriod_choiceMemory,1));
            
            cellID_F_dff_baselinePeriod_choiceOffload_sem = std(cellID_F_dff_baselinePeriod_choiceOffload,1)...
                ./sqrt(size(cellID_F_dff_baselinePeriod_choiceOffload,1));
            
            cellID_F_dff_baselinePeriod_noChoice_sem = std(cellID_F_dff_baselinePeriod_noChoice,1)...
                ./sqrt(size(cellID_F_dff_baselinePeriod_noChoice,1));
            
            
            if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 <= 2
                cellID_F_dff_decisionPeriod_choiceMemory = squeeze(F_dff_decisionPeriod(temp_id2,trialIndex_bool_choiceMemory,:));
                cellID_F_dff_decisionPeriod_choiceOffload = squeeze(F_dff_decisionPeriod(temp_id2,trialIndex_bool_choiceOffload,:));
                cellID_F_dff_decisionPeriod_noChoice = squeeze(F_dff_decisionPeriod(temp_id2,trialIndex_bool_noChoice,:));
                cellID_F_dff_decisionPeriodA_choiceMemory = squeeze(F_dff_decisionPeriodA(temp_id2,trialIndex_bool_choiceMemory,:));
                cellID_F_dff_decisionPeriodA_choiceOffload = squeeze(F_dff_decisionPeriodA(temp_id2,trialIndex_bool_choiceOffload,:));
                cellID_F_dff_decisionPeriodA_noChoice = squeeze(F_dff_decisionPeriodA(temp_id2,trialIndex_bool_noChoice,:));                
            elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 3
                cellID_F_dff_decisionPeriod_choiceMemory = squeeze(F_dff_decisionPeriod(temp_id2,(trialIndex_bool_choiceMemory&trialIndex_bool_afterMemoryCorrect),:));
                cellID_F_dff_decisionPeriod_choiceOffload = squeeze(F_dff_decisionPeriod(temp_id2,(trialIndex_bool_choiceOffload&trialIndex_bool_afterMemoryCorrect),:));
                cellID_F_dff_decisionPeriod_noChoice = squeeze(F_dff_decisionPeriod(temp_id2,(trialIndex_bool_noChoice&trialIndex_bool_afterMemoryCorrect),:));
                cellID_F_dff_decisionPeriodA_choiceMemory = squeeze(F_dff_decisionPeriodA(temp_id2,(trialIndex_bool_choiceMemory&trialIndex_bool_afterMemoryCorrect),:));
                cellID_F_dff_decisionPeriodA_choiceOffload = squeeze(F_dff_decisionPeriodA(temp_id2,(trialIndex_bool_choiceOffload&trialIndex_bool_afterMemoryCorrect),:));
                cellID_F_dff_decisionPeriodA_noChoice = squeeze(F_dff_decisionPeriodA(temp_id2,(trialIndex_bool_noChoice&trialIndex_bool_afterMemoryCorrect),:));                
            elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 4
                cellID_F_dff_decisionPeriod_choiceMemory = squeeze(F_dff_decisionPeriod(temp_id2,(trialIndex_bool_choiceMemory&trialIndex_bool_afterMemoryError),:));
                cellID_F_dff_decisionPeriod_choiceOffload = squeeze(F_dff_decisionPeriod(temp_id2,(trialIndex_bool_choiceOffload&trialIndex_bool_afterMemoryError),:));
                cellID_F_dff_decisionPeriod_noChoice = squeeze(F_dff_decisionPeriod(temp_id2,(trialIndex_bool_noChoice&trialIndex_bool_afterMemoryError),:));
                cellID_F_dff_decisionPeriodA_choiceMemory = squeeze(F_dff_decisionPeriodA(temp_id2,(trialIndex_bool_choiceMemory&trialIndex_bool_afterMemoryError),:));
                cellID_F_dff_decisionPeriodA_choiceOffload = squeeze(F_dff_decisionPeriodA(temp_id2,(trialIndex_bool_choiceOffload&trialIndex_bool_afterMemoryError),:));
                cellID_F_dff_decisionPeriodA_noChoice = squeeze(F_dff_decisionPeriodA(temp_id2,(trialIndex_bool_noChoice&trialIndex_bool_afterMemoryError),:));                
            end
            
            %cellID_F_dff_decisionPeriod_choiceMemory = squeeze(F_dff_decisionPeriod_choiceMemory(temp_id2,:,:));
            %cellID_F_dff_decisionPeriod_choiceOffload = squeeze(F_dff_decisionPeriod_choiceOffload(temp_id2,:,:));
            %cellID_F_dff_decisionPeriod_noChoice = squeeze(F_dff_decisionPeriod_noChoice(temp_id2,:,:));
            
            cellID_F_dff_decisionPeriod_choiceMemory_mean = mean(cellID_F_dff_decisionPeriod_choiceMemory,1);
            cellID_F_dff_decisionPeriod_choiceOffload_mean = mean(cellID_F_dff_decisionPeriod_choiceOffload,1);
            cellID_F_dff_decisionPeriod_noChoice_mean = mean(cellID_F_dff_decisionPeriod_noChoice,1);
            
            cellID_F_dff_decisionPeriod_choiceMemory_sem = std(cellID_F_dff_decisionPeriod_choiceMemory,1)...
                ./sqrt(size(cellID_F_dff_decisionPeriod_choiceMemory,1));
            
            cellID_F_dff_decisionPeriod_choiceOffload_sem = std(cellID_F_dff_decisionPeriod_choiceOffload,1)...
                ./sqrt(size(cellID_F_dff_decisionPeriod_choiceOffload,1));
            
            cellID_F_dff_decisionPeriod_noChoice_sem = std(cellID_F_dff_decisionPeriod_noChoice,1)...
                ./sqrt(size(cellID_F_dff_decisionPeriod_noChoice,1));
            
            cellID_F_dff_decisionPeriodA_choiceMemory_mean = mean(cellID_F_dff_decisionPeriodA_choiceMemory,1);
            cellID_F_dff_decisionPeriodA_choiceOffload_mean = mean(cellID_F_dff_decisionPeriodA_choiceOffload,1);
            cellID_F_dff_decisionPeriodA_noChoice_mean = mean(cellID_F_dff_decisionPeriodA_noChoice,1);
            
            cellID_F_dff_decisionPeriodA_choiceMemory_sem = std(cellID_F_dff_decisionPeriodA_choiceMemory,1)...
                ./sqrt(size(cellID_F_dff_decisionPeriodA_choiceMemory,1));
            
            cellID_F_dff_decisionPeriodA_choiceOffload_sem = std(cellID_F_dff_decisionPeriodA_choiceOffload,1)...
                ./sqrt(size(cellID_F_dff_decisionPeriodA_choiceOffload,1));
            
            cellID_F_dff_decisionPeriodA_noChoice_sem = std(cellID_F_dff_decisionPeriodA_noChoice,1)...
                ./sqrt(size(cellID_F_dff_decisionPeriodA_noChoice,1));
            
            cellID_F_dff_lengthx_period1_location_collapsed = [];
            temp_lengthx_locationIndex = zeros(1, numFrames);
            for tempi=1:numFrames
                temp_dff = cellID_F_dff_lengthx_period1_location{tempi};
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sorted = temp_dff(I,:);
                
                cellID_F_dff_lengthx_period1_location_collapsed = ...
                    [cellID_F_dff_lengthx_period1_location_collapsed;temp_dff_sorted];
                temp_lengthx_locationIndex(tempi) = ...
                    size(cellID_F_dff_lengthx_period1_location_collapsed,1) - ...
                    size(cellID_F_dff_lengthx_period1_location{tempi},1) + 1;
            end
            
            cellID_F_dff_lengthx_period2A_location_collapsed = [];
            for tempi=1:numFrames
                temp_dff1 = cellID_F_dff_lengthx_period1_location{tempi};
                temp_dff1_mean = mean(temp_dff1,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff1_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff1_mean,1);
                end
                
                temp_dff2 = cellID_F_dff_lengthx_period2A_location{tempi};                
                temp_dff2_sorted = temp_dff2(I,:);
                
                cellID_F_dff_lengthx_period2A_location_collapsed = ...
                    [cellID_F_dff_lengthx_period2A_location_collapsed;temp_dff2_sorted];
            end            
            
            cellID_F_dff_lengthx_period2_location_collapsed = [];
            for tempi=1:numFrames
                temp_dff1 = cellID_F_dff_lengthx_period1_location{tempi};
                temp_dff1_mean = mean(temp_dff1,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff1_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff1_mean,1);
                end
                
                temp_dff2 = cellID_F_dff_lengthx_period2_location{tempi};                
                temp_dff2_sorted = temp_dff2(I,:);
                
                cellID_F_dff_lengthx_period2_location_collapsed = ...
                    [cellID_F_dff_lengthx_period2_location_collapsed;temp_dff2_sorted];
            end
                        
            cellID_F_dff_lengthx_period3_location_collapsed = [];
            for tempi=1:numFrames
                temp_dff1 = cellID_F_dff_lengthx_period1_location{tempi};
                temp_dff1_mean = mean(temp_dff1,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff1_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff1_mean,1);
                end
                
                temp_dff2 = cellID_F_dff_lengthx_period3_location{tempi};
                temp_dff2_sorted = temp_dff2(I,:);
                cellID_F_dff_lengthx_period3_location_collapsed = ...
                    [cellID_F_dff_lengthx_period3_location_collapsed;temp_dff2_sorted];
            end
            
            
            if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 0 ...
                    || if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 3 ...
                    || if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 4
                temp_dff = cellID_F_dff_decisionPeriod_choiceMemory;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedA = temp_dff(I,:);
                
                temp_dff = cellID_F_dff_decisionPeriod_choiceOffload;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedB = temp_dff(I,:);
                
                temp_dff = cellID_F_dff_decisionPeriod_noChoice;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedC = temp_dff(I,:);
                
                cellID_F_dff_decisionPeriod_choice_collaped = [...
                    temp_dff_sortedA;temp_dff_sortedB;temp_dff_sortedC];
                temp_choiceIndex = [1, size(cellID_F_dff_decisionPeriod_choiceMemory,1)+1,...
                    size(cellID_F_dff_decisionPeriod_choiceMemory,1)+size(cellID_F_dff_decisionPeriod_choiceOffload,1)+1];
                
            elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 1
                cellID_F_dff_decisionPeriod_choice_collaped_raw = squeeze(F_dff_decisionPeriod(temp_id2,:,:));
                seqIndex;
                
                cellID_F_dff_decisionPeriod_choice_collaped_seq = zeros(sum(numSeq),size(F_dff_decisionPeriod,3));
                for tempi=1:sum(numSeq)
                    temp_range = find(seqIndex == tempi);
                    temp_dff = cellID_F_dff_decisionPeriod_choice_collaped_raw(temp_range,:);
                    cellID_F_dff_decisionPeriod_choice_collaped_seq(tempi,:) = mean(temp_dff,1);
                end
                cellID_F_dff_decisionPeriod_choice_collaped_seq;
                
                [temp_B,temp_I] = sort(offloadingProb_inOne);
                rProbAscendingSortedIndex_seq = temp_I;
                
                cellID_F_dff_decisionPeriod_choice_collaped = ...
                    cellID_F_dff_decisionPeriod_choice_collaped_seq(rProbAscendingSortedIndex_seq,:);
                
                
                %             [temp_B,temp_I] = sort(offloadingProb_inOne);
                %             rProbAscendingSortedIndex_seq = temp_I;
                %             rProbAscendingSortedIndex_trial = [];
                %
                %             for tempi=1:sum(numSeq)
                %                 tempj = rProbAscendingSortedIndex_seq(tempi);
                %                 tempk = find(seqIndex == tempj);
                %                 rProbAscendingSortedIndex_trial = [rProbAscendingSortedIndex_trial tempk];
                %             end
                %             cellID_F_dff_decisionPeriod_choice_collaped = ...
                %                 cellID_F_dff_decisionPeriod_choice_collaped_raw(rProbAscendingSortedIndex_trial,:);
            elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 2
                temp_dff = cellID_F_dff_decisionPeriod_afterMemoryCorrect;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedA = temp_dff(I,:);
                
                temp_dff = cellID_F_dff_decisionPeriod_afterMemoryError;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedB = temp_dff(I,:);
                
                cellID_F_dff_decisionPeriod_choice_collaped = [...
                    temp_dff_sortedA;temp_dff_sortedB];
                temp_choiceIndex = [1, size(cellID_F_dff_decisionPeriod_afterMemoryCorrect,1)+1];
            end
            
            
            if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 0 ...
                    || if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 3 ...
                    || if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 4
                temp_dff = cellID_F_dff_decisionPeriodA_choiceMemory;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedA = temp_dff(I,:);
                
                temp_dff = cellID_F_dff_decisionPeriodA_choiceOffload;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedB = temp_dff(I,:);
                
                temp_dff = cellID_F_dff_decisionPeriodA_noChoice;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedC = temp_dff(I,:);
                
                cellID_F_dff_decisionPeriodA_choice_collaped = [...
                    temp_dff_sortedA;temp_dff_sortedB;temp_dff_sortedC];
                temp_choiceIndex = [1, size(cellID_F_dff_decisionPeriodA_choiceMemory,1)+1,...
                    size(cellID_F_dff_decisionPeriodA_choiceMemory,1)+size(cellID_F_dff_decisionPeriodA_choiceOffload,1)+1];
                
            elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 1
                cellID_F_dff_decisionPeriodA_choice_collaped_raw = squeeze(F_dff_decisionPeriodA(temp_id2,:,:));
                seqIndex;
                
                cellID_F_dff_decisionPeriodA_choice_collaped_seq = zeros(sum(numSeq),size(F_dff_decisionPeriodA,3));
                for tempi=1:sum(numSeq)
                    temp_range = find(seqIndex == tempi);
                    temp_dff = cellID_F_dff_decisionPeriodA_choice_collaped_raw(temp_range,:);
                    cellID_F_dff_decisionPeriodA_choice_collaped_seq(tempi,:) = mean(temp_dff,1);
                end
                cellID_F_dff_decisionPeriodA_choice_collaped_seq;
                
                [temp_B,temp_I] = sort(offloadingProb_inOne);
                rProbAscendingSortedIndex_seq = temp_I;
                
                cellID_F_dff_decisionPeriodA_choice_collaped = ...
                    cellID_F_dff_decisionPeriodA_choice_collaped_seq(rProbAscendingSortedIndex_seq,:);
                
            elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 2
                temp_dff = cellID_F_dff_decisionPeriodA_afterMemoryCorrect;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedA = temp_dff(I,:);
                
                temp_dff = cellID_F_dff_decisionPeriodA_afterMemoryError;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedB = temp_dff(I,:);
                
                cellID_F_dff_decisionPeriodA_choice_collaped = [...
                    temp_dff_sortedA;temp_dff_sortedB];
                temp_choiceIndex = [1, size(cellID_F_dff_decisionPeriodA_afterMemoryCorrect,1)+1];
            end
            
            
            if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 0 ...
                    || if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 3 ...
                    || if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 4
                temp_dff = cellID_F_dff_baselinePeriod_choiceMemory;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedA = temp_dff(I,:);
                
                temp_dff = cellID_F_dff_baselinePeriod_choiceOffload;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedB = temp_dff(I,:);
                
                temp_dff = cellID_F_dff_baselinePeriod_noChoice;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedC = temp_dff(I,:);
                
                cellID_F_dff_baselinePeriod_choice_collaped = [...
                    temp_dff_sortedA;temp_dff_sortedB;temp_dff_sortedC];
                
            elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 1
                cellID_F_dff_baselinePeriod_choice_collaped_raw = squeeze(F_dff_baselinePeriod(temp_id2,:,:));
                seqIndex;
                
                cellID_F_dff_baselinePeriod_choice_collaped_seq = zeros(sum(numSeq),size(F_dff_baselinePeriod,3));
                for tempi=1:sum(numSeq)
                    temp_range = find(seqIndex == tempi);
                    temp_dff = cellID_F_dff_baselinePeriod_choice_collaped_raw(temp_range,:);
                    cellID_F_dff_baselinePeriod_choice_collaped_seq(tempi,:) = mean(temp_dff,1);
                end
                cellID_F_dff_baselinePeriod_choice_collaped_seq;
                
                [temp_B,temp_I] = sort(offloadingProb_inOne);
                rProbAscendingSortedIndex_seq = temp_I;
                
                cellID_F_dff_baselinePeriod_choice_collaped = ...
                    cellID_F_dff_baselinePeriod_choice_collaped_seq(rProbAscendingSortedIndex_seq,:);
                
            elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 2
                temp_dff = cellID_F_dff_baselinePeriod_afterMemoryCorrect;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedA = temp_dff(I,:);
                
                temp_dff = cellID_F_dff_baselinePeriod_afterMemoryError;
                temp_dff_mean = mean(temp_dff,2);
                if if_sortRasterPlot_delay1 == 1
                    [~,I] = sort(temp_dff_mean,'descend');
                elseif if_sortRasterPlot_delay1 == 0
                    I = 1:size(temp_dff_mean,1);
                end
                temp_dff_sortedB = temp_dff(I,:);
                
                cellID_F_dff_baselinePeriod_choice_collaped = [...
                    temp_dff_sortedA;temp_dff_sortedB];
                temp_choiceIndex = [1, size(cellID_F_dff_baselinePeriod_afterMemoryCorrect,1)+1];
            end
            
            if if_plot_location0_seq1 == 1
                
                %trialIndex_lengthx_choiceMemory = trialIndex_lengthx_memoryCorrect;
                %trialIndex_lengthx_choiceOffload = trialIndex_lengthx_memoryCorrect;
                %trialIndex_lengthx_noChoice = trialIndex_lengthx_memoryCorrect;
                
                
                %cellID_F_dff_lengthx_period1 = squeeze(F_dff_lengthx_period1(temp_id2,:,:));
                %cellID_F_dff_lengthx_period2 = squeeze(F_dff_lengthx_period2(temp_id2,:,:));
                cellID_F_dff_lengthx_period1 = squeeze(F_dff_lengthx_period1_raw(temp_id2,:,:));
                cellID_F_dff_lengthx_period2A = squeeze(F_dff_lengthx_period2A_raw(temp_id2,:,:));                
                cellID_F_dff_lengthx_period2 = squeeze(F_dff_lengthx_period2_raw(temp_id2,:,:));
                
                %% cellID_F_dff_lengthx_period12_seq
                cellID_F_dff_lengthx_period1_seq = cell(1,numSeq(plot_lengthFlag));
                cellID_F_dff_lengthx_period2A_seq = cell(1,numSeq(plot_lengthFlag));                
                cellID_F_dff_lengthx_period2_seq = cell(1,numSeq(plot_lengthFlag));
                cellID_F_dff_lengthx_period1_seq_mean = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period1,2));
                cellID_F_dff_lengthx_period2A_seq_mean = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2A,2));                
                cellID_F_dff_lengthx_period2_seq_mean = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2,2));
                cellID_F_dff_lengthx_period1_seq_sem = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period1,2));
                cellID_F_dff_lengthx_period2A_seq_sem = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2A,2));                
                cellID_F_dff_lengthx_period2_seq_sem = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2,2));
                for tempi=1:numSeq(plot_lengthFlag)
                    [~,Locb] = ismember(trialIndex_lengthx_seq{tempi},trialIndex_lengthx_memoryCorrect);
                    Locb(Locb==0) = [];
                    temp1 = trialIndex_lengthx_memoryCorrect(Locb);
                    [~,Locb2] = ismember(temp1,trialIndex_lengthx);
                    
                    cellID_F_dff_lengthx_period1_seq{tempi} = cellID_F_dff_lengthx_period1(Locb2,:);
                    cellID_F_dff_lengthx_period2A_seq{tempi} = cellID_F_dff_lengthx_period2A(Locb2,:);                    
                    cellID_F_dff_lengthx_period2_seq{tempi} = cellID_F_dff_lengthx_period2(Locb2,:);
                    cellID_F_dff_lengthx_period1_seq_mean(tempi,:) = mean(cellID_F_dff_lengthx_period1_seq{tempi},1);
                    cellID_F_dff_lengthx_period2A_seq_mean(tempi,:) = mean(cellID_F_dff_lengthx_period2A_seq{tempi},1);                    
                    cellID_F_dff_lengthx_period2_seq_mean(tempi,:) = mean(cellID_F_dff_lengthx_period2_seq{tempi},1);
                    cellID_F_dff_lengthx_period1_seq_sem(tempi,:) = std(cellID_F_dff_lengthx_period1_seq{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period1_seq{tempi},1));
                    cellID_F_dff_lengthx_period2A_seq_sem(tempi,:) = std(cellID_F_dff_lengthx_period2A_seq{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2A_seq{tempi},1));                    
                    cellID_F_dff_lengthx_period2_seq_sem(tempi,:) = std(cellID_F_dff_lengthx_period2_seq{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2_seq{tempi},1));
                end
                
                %% cellID_F_dff_lengthx_period12_seq_choiceMemory
                cellID_F_dff_lengthx_period1_seq_choiceMemory = cell(1,numSeq(plot_lengthFlag));
                cellID_F_dff_lengthx_period2A_seq_choiceMemory = cell(1,numSeq(plot_lengthFlag));                
                cellID_F_dff_lengthx_period2_seq_choiceMemory = cell(1,numSeq(plot_lengthFlag));
                cellID_F_dff_lengthx_period1_seq_mean_choiceMemory = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period1,2));
                cellID_F_dff_lengthx_period2A_seq_mean_choiceMemory = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2A,2));                
                cellID_F_dff_lengthx_period2_seq_mean_choiceMemory = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2,2));
                cellID_F_dff_lengthx_period1_seq_sem_choiceMemory = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period1,2));
                cellID_F_dff_lengthx_period2A_seq_sem_choiceMemory = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2A,2));                
                cellID_F_dff_lengthx_period2_seq_sem_choiceMemory = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2,2));
                for tempi=1:numSeq(plot_lengthFlag)
                    %[~,Locb] = ismember(trialIndex_lengthx_seq{tempi},trialIndex_lengthx_choiceMemory);
                    [~,Locb] = ismember(trialIndex_lengthx_seq{tempi},trialIndex_lengthx_choiceMemoryCorrect);
                    Locb(Locb==0) = [];
                    temp1 = trialIndex_lengthx_choiceMemoryCorrect(Locb);
                    [~,Locb2] = ismember(temp1,trialIndex_lengthx);
                    
                    cellID_F_dff_lengthx_period1_seq_choiceMemory{tempi} = cellID_F_dff_lengthx_period1(Locb2,:);
                    cellID_F_dff_lengthx_period2A_seq_choiceMemory{tempi} = cellID_F_dff_lengthx_period2A(Locb2,:);
                    cellID_F_dff_lengthx_period2_seq_choiceMemory{tempi} = cellID_F_dff_lengthx_period2(Locb2,:);                    
                    cellID_F_dff_lengthx_period1_seq_mean_choiceMemory(tempi,:) = mean(cellID_F_dff_lengthx_period1_seq_choiceMemory{tempi},1);
                    cellID_F_dff_lengthx_period2A_seq_mean_choiceMemory(tempi,:) = mean(cellID_F_dff_lengthx_period2A_seq_choiceMemory{tempi},1);                    
                    cellID_F_dff_lengthx_period2_seq_mean_choiceMemory(tempi,:) = mean(cellID_F_dff_lengthx_period2_seq_choiceMemory{tempi},1);
                    cellID_F_dff_lengthx_period1_seq_sem_choiceMemory(tempi,:) = std(cellID_F_dff_lengthx_period1_seq_choiceMemory{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period1_seq_choiceMemory{tempi},1));
                    cellID_F_dff_lengthx_period2A_seq_sem_choiceMemory(tempi,:) = std(cellID_F_dff_lengthx_period2A_seq_choiceMemory{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2A_seq_choiceMemory{tempi},1));                    
                    cellID_F_dff_lengthx_period2_seq_sem_choiceMemory(tempi,:) = std(cellID_F_dff_lengthx_period2_seq_choiceMemory{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2_seq_choiceMemory{tempi},1));
                end
                
                %% cellID_F_dff_lengthx_period12_seq_choiceOffload
                cellID_F_dff_lengthx_period1_seq_choiceOffload = cell(1,numSeq(plot_lengthFlag));
                cellID_F_dff_lengthx_period2A_seq_choiceOffload = cell(1,numSeq(plot_lengthFlag));                
                cellID_F_dff_lengthx_period2_seq_choiceOffload = cell(1,numSeq(plot_lengthFlag));
                cellID_F_dff_lengthx_period1_seq_mean_choiceOffload = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period1,2));
                cellID_F_dff_lengthx_period2A_seq_mean_choiceOffload = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2A,2));                
                cellID_F_dff_lengthx_period2_seq_mean_choiceOffload = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2,2));
                cellID_F_dff_lengthx_period1_seq_sem_choiceOffload = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period1,2));
                cellID_F_dff_lengthx_period2A_seq_sem_choiceOffload = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2A,2));                
                cellID_F_dff_lengthx_period2_seq_sem_choiceOffload = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2,2));
                for tempi=1:numSeq(plot_lengthFlag)
                    [~,Locb] = ismember(trialIndex_lengthx_seq{tempi},trialIndex_lengthx_choiceOffload);
                    Locb(Locb==0) = [];
                    temp1 = trialIndex_lengthx_choiceOffload(Locb);
                    [~,Locb2] = ismember(temp1,trialIndex_lengthx);
                    
                    cellID_F_dff_lengthx_period1_seq_choiceOffload{tempi} = cellID_F_dff_lengthx_period1(Locb2,:);
                    cellID_F_dff_lengthx_period2A_seq_choiceOffload{tempi} = cellID_F_dff_lengthx_period2A(Locb2,:);                    
                    cellID_F_dff_lengthx_period2_seq_choiceOffload{tempi} = cellID_F_dff_lengthx_period2(Locb2,:);
                    cellID_F_dff_lengthx_period1_seq_mean_choiceOffload(tempi,:) = mean(cellID_F_dff_lengthx_period1_seq_choiceOffload{tempi},1);
                    cellID_F_dff_lengthx_period2A_seq_mean_choiceOffload(tempi,:) = mean(cellID_F_dff_lengthx_period2A_seq_choiceOffload{tempi},1);                    
                    cellID_F_dff_lengthx_period2_seq_mean_choiceOffload(tempi,:) = mean(cellID_F_dff_lengthx_period2_seq_choiceOffload{tempi},1);
                    cellID_F_dff_lengthx_period1_seq_sem_choiceOffload(tempi,:) = std(cellID_F_dff_lengthx_period1_seq_choiceOffload{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period1_seq_choiceOffload{tempi},1));
                    cellID_F_dff_lengthx_period2A_seq_sem_choiceOffload(tempi,:) = std(cellID_F_dff_lengthx_period2A_seq_choiceOffload{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2A_seq_choiceOffload{tempi},1));                    
                    cellID_F_dff_lengthx_period2_seq_sem_choiceOffload(tempi,:) = std(cellID_F_dff_lengthx_period2_seq_choiceOffload{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2_seq_choiceOffload{tempi},1));
                end
                
                %% cellID_F_dff_lengthx_period12_seq_noChoice
                cellID_F_dff_lengthx_period1_seq_noChoice = cell(1,numSeq(plot_lengthFlag));
                cellID_F_dff_lengthx_period2A_seq_noChoice = cell(1,numSeq(plot_lengthFlag));                
                cellID_F_dff_lengthx_period2_seq_noChoice = cell(1,numSeq(plot_lengthFlag));
                cellID_F_dff_lengthx_period1_seq_mean_noChoice = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period1,2));
                cellID_F_dff_lengthx_period2A_seq_mean_noChoice = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2A,2));                
                cellID_F_dff_lengthx_period2_seq_mean_noChoice = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2,2));
                cellID_F_dff_lengthx_period1_seq_sem_noChoice = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period1,2));
                cellID_F_dff_lengthx_period2A_seq_sem_noChoice = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2A,2));                
                cellID_F_dff_lengthx_period2_seq_sem_noChoice = zeros(numSeq(plot_lengthFlag),size(cellID_F_dff_lengthx_period2,2));
                for tempi=1:numSeq(plot_lengthFlag)
                    %[~,Locb] = ismember(trialIndex_lengthx_seq{tempi},trialIndex_lengthx_noChoice);
                    [~,Locb] = ismember(trialIndex_lengthx_seq{tempi},trialIndex_lengthx_noChoiceCorrect);
                    Locb(Locb==0) = [];
                    temp1 = trialIndex_lengthx_noChoiceCorrect(Locb);
                    [~,Locb2] = ismember(temp1,trialIndex_lengthx);
                    
                    cellID_F_dff_lengthx_period1_seq_noChoice{tempi} = cellID_F_dff_lengthx_period1(Locb2,:);
                    cellID_F_dff_lengthx_period2A_seq_noChoice{tempi} = cellID_F_dff_lengthx_period2A(Locb2,:);                    
                    cellID_F_dff_lengthx_period2_seq_noChoice{tempi} = cellID_F_dff_lengthx_period2(Locb2,:);
                    cellID_F_dff_lengthx_period1_seq_mean_noChoice(tempi,:) = mean(cellID_F_dff_lengthx_period1_seq_noChoice{tempi},1);
                    cellID_F_dff_lengthx_period2A_seq_mean_noChoice(tempi,:) = mean(cellID_F_dff_lengthx_period2A_seq_noChoice{tempi},1);                    
                    cellID_F_dff_lengthx_period2_seq_mean_noChoice(tempi,:) = mean(cellID_F_dff_lengthx_period2_seq_noChoice{tempi},1);
                    cellID_F_dff_lengthx_period1_seq_sem_noChoice(tempi,:) = std(cellID_F_dff_lengthx_period1_seq_noChoice{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period1_seq_noChoice{tempi},1));
                    cellID_F_dff_lengthx_period2A_seq_sem_noChoice(tempi,:) = std(cellID_F_dff_lengthx_period2A_seq_noChoice{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2A_seq_noChoice{tempi},1));                    
                    cellID_F_dff_lengthx_period2_seq_sem_noChoice(tempi,:) = std(cellID_F_dff_lengthx_period2_seq_noChoice{tempi},1)./sqrt(size(cellID_F_dff_lengthx_period2_seq_noChoice{tempi},1));
                end
            end
            
            
            periodSkipInterval = 2;%3
            %temp_ticks_period1 = lengthx_period1_interval;
            %temp_ticks_period2 = lengthx_period1_interval(end)+periodSkipInterval+lengthx_period2_interval;
            %temp_ticks_period3 = lengthx_period1_interval(end)+periodSkipInterval+lengthx_period2_interval(end)+periodSkipInterval+lengthx_period3_interval;
            
            temp_ticks_period1 = lengthx_period1_interval;
            temp_ticks_period2A = temp_ticks_period1(end)+periodSkipInterval+lengthx_period2A_interval;            
            temp_ticks_period2 = temp_ticks_period2A(end)+periodSkipInterval+lengthx_period2_interval;
            temp_ticks_period3 = temp_ticks_period2(end)+periodSkipInterval+lengthx_period3_interval;
            
            
            temp_ticks_baselinePeriod = baselinePeriod_interval;
            temp_ticks_decisionPeriodA = temp_ticks_baselinePeriod(end)+periodSkipInterval+decisionPeriodA_interval;            
            temp_ticks_decisionPeriod = temp_ticks_decisionPeriodA(end)+periodSkipInterval+decisionPeriod_interval;            
            %temp_ticks_decisionPeriod = temp_ticks_baselinePeriod(end)+periodSkipInterval+decisionPeriod_interval;
            
            % if plot_lengthFlag == 1
            if if_plot_location0_seq1 == 0
                xlim_padSize1 = 2;
                xlim_padSize2 = 20;%9-->15-->18-->20
                
                multi_rgbColor = ...
                    [228,26,28;
                    55,126,184;
                    77,175,74;
                    152,78,163;
                    255,127,0;
                    255,255,51]/255;
                
                backgounrdColor = [1 1 1]*0.825;%0.875
                
                [min_period1,max_period1] = bounds(cellID_F_dff_lengthx_period1_location_mean,'all');                
                [min_period2,max_period2] = bounds(cellID_F_dff_lengthx_period2_location_mean,'all');
                [min_period3,max_period3] = bounds(cellID_F_dff_lengthx_period3_location_mean,'all');
                
                min_period12 = min(min_period1,min_period2) - max([cellID_F_dff_lengthx_period1_location_sem cellID_F_dff_lengthx_period2_location_sem],[],'all');
                max_period12 = max(max_period1,max_period2) + max([cellID_F_dff_lengthx_period1_location_sem cellID_F_dff_lengthx_period2_location_sem],[],'all');
                
                if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 ~= 2
                    [min_decisionPeriodA,max_decisionPeriodA] = bounds([cellID_F_dff_decisionPeriodA_choiceMemory_mean;cellID_F_dff_decisionPeriodA_choiceOffload_mean],'all');                    
                    [min_decisionPeriod,max_decisionPeriod] = bounds([cellID_F_dff_decisionPeriod_choiceMemory_mean;cellID_F_dff_decisionPeriod_choiceOffload_mean],'all');
                    [min_baselinePeriod,max_baselinePeriod] = bounds([cellID_F_dff_baselinePeriod_choiceMemory_mean;cellID_F_dff_baselinePeriod_choiceOffload_mean],'all');
                    min_baselineDecisionPeriod = min([min_decisionPeriodA,min_decisionPeriod,min_baselinePeriod]) - max([cellID_F_dff_decisionPeriodA_choiceMemory_sem cellID_F_dff_decisionPeriod_choiceMemory_sem cellID_F_dff_decisionPeriod_choiceOffload_sem],[],'all');
                    max_baselineDecisionPeriod = max([max_decisionPeriodA,max_decisionPeriod,max_baselinePeriod]) + max([cellID_F_dff_decisionPeriodA_choiceMemory_sem cellID_F_dff_decisionPeriod_choiceMemory_sem cellID_F_dff_decisionPeriod_choiceOffload_sem],[],'all');
                elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 2
                    [min_decisionPeriodA,max_decisionPeriodA] = bounds([cellID_F_dff_decisionPeriodA_afterMemoryCorrect_mean;cellID_F_dff_decisionPeriodA_afterMemoryError_mean],'all');                    
                    [min_decisionPeriod,max_decisionPeriod] = bounds([cellID_F_dff_decisionPeriod_afterMemoryCorrect_mean;cellID_F_dff_decisionPeriod_afterMemoryError_mean],'all');
                    [min_baselinePeriod,max_baselinePeriod] = bounds([cellID_F_dff_baselinePeriod_afterMemoryCorrect_mean;cellID_F_dff_baselinePeriod_afterMemoryError_mean],'all');
                    min_baselineDecisionPeriod = min([min_decisionPeriodA,min_decisionPeriod,min_baselinePeriod]) - max([cellID_F_dff_decisionPeriodA_afterMemoryCorrect_sem cellID_F_dff_decisionPeriod_afterMemoryCorrect_sem cellID_F_dff_decisionPeriod_afterMemoryError_sem],[],'all');
                    max_baselineDecisionPeriod = max([max_decisionPeriodA,max_decisionPeriod,max_baselinePeriod]) + max([cellID_F_dff_decisionPeriodA_afterMemoryCorrect_sem cellID_F_dff_decisionPeriod_afterMemoryCorrect_sem cellID_F_dff_decisionPeriod_afterMemoryError_sem],[],'all');
                end
                
                %fig1 = figure('Name','Fig1','NumberTitle','off');
                %fig1 = figure('Name',sprintf('fig%d, cell id %d',tempIndex,temp_cellIndex_suite2p),'NumberTitle','off');
                fig1 = figure('Name',sprintf('fig%d, cell id %d, FOV id %s',tempIndex,temp_cellIndex_suite2p,currentSession(18:22)),'NumberTitle','off');                
                %set(gcf,'Position',[35+260 35+0 1100 850]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                %set(gcf,'Position',[35+260 35+0 1100 1000]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                %set(gcf,'Position',[35+0 35+0 1630 1000]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                set(gcf,'Position',[35+0 42+0 1630 840]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point                
                %t = tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
                %t = tiledlayout(4,11,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
                t = tiledlayout(4,11,'TileSpacing','tight','Padding','tight'); %#ok<*NASGU>
                
                t.Title.String = sprintf('suite2p cell id=%d                      %s',...
                    temp_cellIndex_suite2p,currentSession(18:end));
                t.Title.FontSize = 16;
                t.Title.Interpreter = 'none';
                
                %set(gcf,'color',backgounrdColor);
                
                %% subplot 1, raster plot in period1
                %nexttile
                nexttile([1,11]);
                x = [temp_ticks_period1(1) temp_ticks_period1(end)];
                y = [1,size(cellID_F_dff_lengthx_period1_location_collapsed,1)];
                C = cellID_F_dff_lengthx_period1_location_collapsed;
                imagesc(x,y,C);
                hold on
                
                x = [temp_ticks_period2A(1) temp_ticks_period2A(end)];
                y = [1,size(cellID_F_dff_lengthx_period1_location_collapsed,1)];
                C = cellID_F_dff_lengthx_period2A_location_collapsed;
                imagesc(x,y,C);
                hold on
                
                x = [temp_ticks_period2(1) temp_ticks_period2(end)];
                y = [1,size(cellID_F_dff_lengthx_period1_location_collapsed,1)];
                C = cellID_F_dff_lengthx_period2_location_collapsed;
                imagesc(x,y,C);
                hold on
                
                x = [temp_ticks_period3(1) temp_ticks_period3(end)];
                y = [1,size(cellID_F_dff_lengthx_period1_location_collapsed,1)];
                C = cellID_F_dff_lengthx_period3_location_collapsed;
                imagesc(x,y,C);
                hold on
                
                
                for tempi=1:length(lengthx_period1_interval)
                    plot(temp_ticks_period1(tempi)*[1 1],[1 size(cellID_F_dff_lengthx_period1_location_collapsed,1)],...
                        '-','LineWidth',1,'Color',[1 1 1])
                    hold on
                end
                
                for tempi=1:length(lengthx_period2A_interval)
                    plot(temp_ticks_period2A(tempi)*[1 1],[1 size(cellID_F_dff_lengthx_period2A_location_collapsed,1)],...
                        '-','LineWidth',1,'Color',[1 1 1]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period2_interval)
                    plot(temp_ticks_period2(tempi)*[1 1],[1 size(cellID_F_dff_lengthx_period2_location_collapsed,1)],...
                        '-','LineWidth',1,'Color',[1 1 1]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period3_interval)
                    plot(temp_ticks_period3(tempi)*[1 1],[1 size(cellID_F_dff_lengthx_period3_location_collapsed,1)],...
                        '-','LineWidth',1,'Color',[1 1 1]);
                    hold on
                end
                
                for tempi=1:numFrames
                    plot(temp_ticks_period3(end)+[0.5 5],(temp_lengthx_locationIndex(tempi)-0.5)*[1 1],...
                        '-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                end
                
                for tempi=1:numFrames
                    %plot([0.5 temp_ticks_period3(end)+0.5],(temp_lengthx_locationIndex(tempi)-0.5)*[1 1],...
                    %   '-','LineWidth',1,'Color',[0 0 0]);
                    plot([0.5 temp_ticks_period3(end)+0.5],(temp_lengthx_locationIndex(tempi)-0.5)*[1 1],...
                        ':','LineWidth',0.5,'Color',[1 1 1]);
                    hold on
                end
                
                for tempi=1:plot_lengthFlag
                    x = [lengthx_period1_interval(3+tempi-1) 1;...
                        lengthx_period1_interval(3+tempi-1) size(cellID_F_dff_lengthx_period1_location_collapsed,1);...
                        lengthx_period1_interval(3+tempi-1)+6 size(cellID_F_dff_lengthx_period1_location_collapsed,1);...
                        lengthx_period1_interval(3+tempi-1)+6 1];
                    y = [1 2 3 4];
                    patch('Faces',y,'Vertices',x,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.2,'EdgeColor','none');
                    hold on
                end
                
                set(gca,'color',backgounrdColor);
                xlim([1-xlim_padSize1 temp_ticks_period3(end)+xlim_padSize2]);
                ylim([1-1 size(cellID_F_dff_lengthx_period1_location_collapsed,1)+1]);
                set(gca,'xticklabel',[])
                set(gca,'xtick',[])
                
                xticks([temp_ticks_period1 temp_ticks_period2A temp_ticks_period2 temp_ticks_period3]);
                if plot_lengthFlag == 1
                    %xticklabels({'PreFixation','Fixation','T1','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
                    %xticklabels({'PreFixation','Fixation','T1','Delay1',' ',' ','Decision','','','Submit',''}); 
                    xticklabels({'PreFixation','Fixation','T1','Delay1','','','ChoiceCue','','','Decision','','','Submit',''}); 
                elseif plot_lengthFlag == 2
                    xticklabels({'PreFixation','Fixation','T1','T2','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 3
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 4
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','T4','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});                    
                end
                
                text(temp_ticks_period1(end),size(cellID_F_dff_lengthx_period1_location_collapsed,1),'/','fontsize',15);
                text(temp_ticks_period2A(1),size(cellID_F_dff_lengthx_period1_location_collapsed,1),'/','fontsize',15);
                
                text(temp_ticks_period2A(end),size(cellID_F_dff_lengthx_period1_location_collapsed,1),'/','fontsize',15);
                text(temp_ticks_period2(1),size(cellID_F_dff_lengthx_period1_location_collapsed,1),'/','fontsize',15);
                
                text(temp_ticks_period2(end),size(cellID_F_dff_lengthx_period1_location_collapsed,1),'/','fontsize',15);
                text(temp_ticks_period3(1),size(cellID_F_dff_lengthx_period1_location_collapsed,1),'/','fontsize',15);
                
                set(gca, 'FontSize', 11)
                set(gca,'box','off');% 取消右、上边框
                %ylabel(sprintf('Trials\ncellID=%d',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
                ylabel(sprintf('Trials'), 'FontSize', 14, 'FontWeight', 'bold');
                
                colorbar;
                
                %% subplot 2, trial average in period1
                %nexttile
                nexttile([1,11]);
                h_line = [];
                for tempi=1:numFrames
                    x = temp_ticks_period1(1):temp_ticks_period1(end);
                    y = cellID_F_dff_lengthx_period1_location_mean(tempi,:);
                    h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))];
                    hold on
                    y_sem = cellID_F_dff_lengthx_period1_location_sem(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                end
                
                for tempi=1:numFrames
                    x = temp_ticks_period2A(1):temp_ticks_period2A(end);
                    y = cellID_F_dff_lengthx_period2A_location_mean(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period2A_location_sem(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                end
                
                for tempi=1:numFrames
                    x = temp_ticks_period2(1):temp_ticks_period2(end);
                    y = cellID_F_dff_lengthx_period2_location_mean(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period2_location_sem(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                end
                
                for tempi=1:numFrames
                    x = temp_ticks_period3(1):temp_ticks_period3(end);
                    y = cellID_F_dff_lengthx_period3_location_mean(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period3_location_sem(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period1_interval)
                    plot(temp_ticks_period1(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period2A_interval)
                    plot(temp_ticks_period2A(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);                    
                    hold on
                end
                
                for tempi=1:length(lengthx_period2_interval)
                    plot(temp_ticks_period2(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);                    
                    hold on
                end
                
                for tempi=1:length(lengthx_period3_interval)
                    plot(temp_ticks_period3(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:plot_lengthFlag
                    x = [lengthx_period1_interval(3+tempi-1) min_period12;...
                        lengthx_period1_interval(3+tempi-1) max_period12;...
                        lengthx_period1_interval(3+tempi-1)+6 max_period12;...
                        lengthx_period1_interval(3+tempi-1)+6 min_period12];
                    y = [1 2 3 4];
                    patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.1,'EdgeColor','none');
                    hold on
                end
                
                le = legend(h_line,'location1','location2','location3','location4','location5','location6','Location','northeast','fontsize',10);%10                
                le.ItemTokenSize = ones(1,6)*20;
                le.Color = backgounrdColor;
                a = 1;
                
                set(gca,'color',backgounrdColor);
                xlim([1-xlim_padSize1 temp_ticks_period3(end)+xlim_padSize2]);
                ylim([min_period12,max_period12]);
                set(gca,'xticklabel',[])
                set(gca,'xtick',[])
                
                xticks([temp_ticks_period1 temp_ticks_period2A temp_ticks_period2 temp_ticks_period3]);
                if plot_lengthFlag == 1
                    %xticklabels({'PreFixation','Fixation','T1','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
                    %xticklabels({'PreFixation','Fixation','T1','Delay1',' ',' ','Decision','','','Submit',''}); 
                    xticklabels({'PreFixation','Fixation','T1','Delay1','','','ChoiceCue','','','Decision','','','Submit',''}); 
                elseif plot_lengthFlag == 2
                    xticklabels({'PreFixation','Fixation','T1','T2','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 3
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 4
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','T4','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});                    
                end
                
                text(temp_ticks_period1(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period2A(1),min_period12,'/','fontsize',15);
                
                text(temp_ticks_period2A(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period2(1),min_period12,'/','fontsize',15);
                
                text(temp_ticks_period2(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period3(1),min_period12,'/','fontsize',15);
                
                set(gca, 'FontSize', 11)
                set(gca,'box','off');% 取消右、上边框
                %ylabel(sprintf('dF / F\ncellID=%d',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
                ylabel(sprintf('dF / F'), 'FontSize', 14, 'FontWeight', 'bold');
                
                
                %% subplot 3, raster plot in decision-making
                %nexttile
                nexttile([1,7]);
                x = [temp_ticks_baselinePeriod(1) temp_ticks_baselinePeriod(end)];
                y = [1,size(cellID_F_dff_baselinePeriod_choice_collaped,1)];
                C = cellID_F_dff_baselinePeriod_choice_collaped;
                imagesc(x,y,C);
                hold on
                
                x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];
                y = [1,size(cellID_F_dff_decisionPeriodA_choice_collaped,1)];
                C = cellID_F_dff_decisionPeriodA_choice_collaped;
                imagesc(x,y,C);
                hold on
                
                x = [temp_ticks_decisionPeriod(1) temp_ticks_decisionPeriod(end)];
                y = [1,size(cellID_F_dff_decisionPeriod_choice_collaped,1)];
                C = cellID_F_dff_decisionPeriod_choice_collaped;
                imagesc(x,y,C);
                hold on
                
                
                for tempi=1:length(baselinePeriod_interval)
                    plot(temp_ticks_baselinePeriod(tempi)*[1 1],[1 size(cellID_F_dff_baselinePeriod_choice_collaped,1)],...
                        '-','LineWidth',1,'Color',[1 1 1])
                    hold on
                end
                
                for tempi=1:length(decisionPeriodA_interval)
                    plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[1 size(cellID_F_dff_decisionPeriodA_choice_collaped,1)],...
                        '-','LineWidth',1,'Color',[1 1 1]);
                    hold on
                end
                
                for tempi=1:length(decisionPeriod_interval)
                    plot(temp_ticks_decisionPeriod(tempi)*[1 1],[1 size(cellID_F_dff_decisionPeriod_choice_collaped,1)],...
                        '-','LineWidth',1,'Color',[1 1 1]);
                    hold on
                end
                
                
                if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 0 ...
                        || if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 3 ...
                        || if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 4
                    tempColor = [...
                        0.25,0.75,0.25;
                        0.75,0.25,0.25;
                        0.25,0.25,0.25];
                    %for tempi=1:2
                    for tempi=1:3
                        plot([temp_ticks_decisionPeriod(end)+0.5, temp_ticks_decisionPeriod(end)+5]...
                            ,(temp_choiceIndex(tempi)-0.5)*[1 1],...
                            '-','LineWidth',1,'Color',tempColor(tempi,:));
                        hold on
                        
                        plot([0.5 temp_ticks_decisionPeriod(end)+0.5],(temp_choiceIndex(tempi)-0.5)*[1 1],...
                            ':','LineWidth',0.5,'Color',[1 1 1]);%0.5
                        hold on
                    end
                elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 2
                    tempColor = [...
                        0.25,0.25,0.25;
                        0.75,0.25,0.25];
                    for tempi=1:2
                        plot([temp_ticks_decisionPeriod(end)+0.5, temp_ticks_decisionPeriod(end)+5]...
                            ,(temp_choiceIndex(tempi)-0.5)*[1 1],...
                            '-','LineWidth',1,'Color',tempColor(tempi,:));
                        hold on
                        
                        plot([0.5 temp_ticks_decisionPeriod(end)+0.5],(temp_choiceIndex(tempi)-0.5)*[1 1],...
                            ':','LineWidth',0.5,'Color',[1 1 1]);%0.5
                        hold on
                    end
                end
                
                set(gca,'color',backgounrdColor);
                xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
                ylim([1-1 size(cellID_F_dff_decisionPeriod_choice_collaped,1)+1]);
                set(gca,'xticklabel',[])
                set(gca,'xtick',[])
                
                xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
                %xticklabels({'PreFixation','Fixation','','','Decision',''});
                xticklabels({'PreFixation','Fixation','','','ChoiceCue','','','Decision',''});
                
                
                text(temp_ticks_baselinePeriod(end),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',15);
                text(temp_ticks_decisionPeriodA(1),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',15);
                
                text(temp_ticks_decisionPeriodA(end),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',15);
                text(temp_ticks_decisionPeriod(1),size(cellID_F_dff_baselinePeriod_choice_collaped,1),'/','fontsize',15);
                
                set(gca, 'FontSize', 11)
                set(gca,'box','off');% 取消右、上边框
                if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 ~= 1
                    %ylabel(sprintf('Trials\ncellID=%d',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
                    ylabel(sprintf('Trials'), 'FontSize', 14, 'FontWeight', 'bold');                    
                elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 1
                    %ylabel(sprintf('Seqs (offload sort)\ncellID=%d',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
                    ylabel(sprintf('Seqs (offload sort)'), 'FontSize', 14, 'FontWeight', 'bold');                    
                end
                
                colorbar;
                %% subplot 4, roi anatomy info
                ax4 = nexttile([2,4]);
                temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
                
                temp_roi_npix = double(temp_roi_stat.npix);
                
                temp_roi_xpix = double(temp_roi_stat.xpix);
                temp_roi_ypix = double(temp_roi_stat.ypix);
                
                [temp_x_min,temp_x_max] = bounds(temp_roi_xpix);
                [temp_y_min,temp_y_max] = bounds(temp_roi_ypix);
                
                temp_delta = 50;%20
                temp_y_min = max(temp_y_min - temp_delta,1);
                temp_y_max = min(temp_y_max + temp_delta,512);
                temp_x_min = max(temp_x_min - temp_delta,1);
                temp_x_max = min(temp_x_max + temp_delta,512);
                
                
                temp_template = template(temp_y_min:temp_y_max,temp_x_min:temp_x_max);
                temp_template = temp_template / max(temp_template,[],'all');
                
                %temp_template2 = (temp_template.^0.7)*255;
                temp_template2 = (temp_template.^0.55)*255;
                
                image(temp_template2);
                hold on
                colormap(ax4,gray);
                
                temp_I = false(512,512);
                for tempi=1:length(temp_roi_xpix)
                    temp_I(temp_roi_ypix(tempi),temp_roi_xpix(tempi)) = true;
                end
                
                %se = strel('square',2);
                se = strel('disk',3);%1,2
                temp_I = imdilate(temp_I,se);                
                
                temp_I2 = temp_I(temp_y_min:temp_y_max,temp_x_min:temp_x_max);
                B = bwboundaries(temp_I2,'noholes');
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:,2)+1,boundary(:,1)+1,'color',[0.4660 0.6740 0.1880], 'LineWidth', 1.5);
                    hold on
                end
                
                axis equal;
                xticks([1 (temp_x_max-temp_x_min)]);
                xticklabels({sprintf('%d',temp_x_min),sprintf('%d',temp_x_max)});
                yticks([1 (temp_y_max-temp_y_min)]);
                yticklabels({sprintf('%d',temp_y_min),sprintf('%d',temp_y_max)})
                %set(gca,'xtick',[]);
                %set(gca,'ytick',[]);
                %set(gca,'Visible','off');
                title(sprintf('npix=%d',temp_roi_npix),'FontSize',12);
                
                %% subplot 5, trial average in decision-making
                %nexttile
                nexttile([1,7]);
                color_choiceMemory = [0.4660 0.6740 0.1880];
                color_choiceOffload = [0.6350 0.0780 0.1840];
                
                if_plot_noChoice = 0;
                
                if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 ~= 2
                    x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
                    y = cellID_F_dff_baselinePeriod_choiceMemory_mean;
                    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemory);
                    hold on
                    y_sem = cellID_F_dff_baselinePeriod_choiceMemory_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
                    y = cellID_F_dff_baselinePeriod_choiceOffload_mean;
                    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffload);
                    hold on
                    y_sem = cellID_F_dff_baselinePeriod_choiceOffload_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    if if_plot_noChoice == 1
                    x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
                    y = cellID_F_dff_baselinePeriod_noChoice_mean;
                    plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.25,0.25]);
                    hold on
                    y_sem = cellID_F_dff_baselinePeriod_noChoice_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    end
                    
                    
                    x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                    y = cellID_F_dff_decisionPeriodA_choiceMemory_mean;
                    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemory);
                    hold on
                    y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                    y = cellID_F_dff_decisionPeriodA_choiceOffload_mean;
                    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffload);
                    hold on
                    y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    if if_plot_noChoice == 1
                    x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                    y = cellID_F_dff_decisionPeriodA_noChoice_mean;
                    plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.25,0.25]);
                    hold on
                    y_sem = cellID_F_dff_decisionPeriodA_noChoice_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    end
                    
                    h_line = [];
                    
                    x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
                    y = cellID_F_dff_decisionPeriod_choiceMemory_mean;
                    h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemory)];
                    hold on
                    y_sem = cellID_F_dff_decisionPeriod_choiceMemory_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceMemory,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
                    y = cellID_F_dff_decisionPeriod_choiceOffload_mean;
                    h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffload)];
                    hold on
                    y_sem = cellID_F_dff_decisionPeriod_choiceOffload_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    if if_plot_noChoice == 1
                    x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
                    y = cellID_F_dff_decisionPeriod_noChoice_mean;
                    h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])];
                    hold on
                    y_sem = cellID_F_dff_decisionPeriod_noChoice_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    end
                    
                    
                    if if_plot_memoryMetaMismatch == 1
                        test_temp_plot_memoryMetaMismatch_singleNeuron_Name_v = autoGetFunName_myScripts('test_temp_plot_memoryMetaMismatch_singleNeuron',targetPATH);
                        script_plot_memoryMetaMismatch_singleNeuron = str2func(test_temp_plot_memoryMetaMismatch_singleNeuron_Name_v);
                        
                        script_plot_memoryMetaMismatch_singleNeuron();
                    end
                    
                elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 2
                    
                    x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
                    y = cellID_F_dff_baselinePeriod_afterMemoryCorrect_mean;
                    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemory);
                    hold on
                    y_sem = cellID_F_dff_baselinePeriod_afterMemoryCorrect_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
                    y = cellID_F_dff_baselinePeriod_afterMemoryError_mean;
                    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffload);
                    hold on
                    y_sem = cellID_F_dff_baselinePeriod_afterMemoryError_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                    y = cellID_F_dff_decisionPeriodA_afterMemoryCorrect_mean;
                    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemory);
                    hold on
                    y_sem = cellID_F_dff_decisionPeriodA_afterMemoryCorrect_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                    y = cellID_F_dff_decisionPeriodA_afterMemoryError_mean;
                    plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffload);
                    hold on
                    y_sem = cellID_F_dff_decisionPeriodA_afterMemoryError_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    
                    h_line = [];
                    
                    x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
                    y = cellID_F_dff_decisionPeriod_afterMemoryCorrect_mean;
                    h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',color_choiceMemory)];
                    hold on
                    y_sem = cellID_F_dff_decisionPeriod_afterMemoryCorrect_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
                    y = cellID_F_dff_decisionPeriod_afterMemoryError_mean;
                    h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',color_choiceOffload)];
                    hold on
                    y_sem = cellID_F_dff_decisionPeriod_afterMemoryError_sem;
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',color_choiceOffload,'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                end
                
                for tempi=1:length(baselinePeriod_interval)
                    plot(temp_ticks_baselinePeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:length(decisionPeriodA_interval)
                    plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                        '-','LineWidth',1,'Color',[0 0 0]);                    
                    hold on
                end
                
                for tempi=1:length(decisionPeriod_interval)
                    plot(temp_ticks_decisionPeriod(tempi)*[1 1],[min_baselineDecisionPeriod max_baselineDecisionPeriod],...
                        '-','LineWidth',1,'Color',[0 0 0]);                    
                    hold on
                end
                
                if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 <= 1
                    if if_plot_noChoice == 1
                        le = legend(h_line,'choiceMemory','choiceOffload','noChoice','Location','northeast','fontsize',10);%10
                    elseif if_plot_noChoice == 0
                        %le = legend(h_line,'choiceMemory','choiceOffload','Location','northeast','fontsize',10,'NumColumns',2);%10
                        le = legend(h_line,'choiceMemory','choiceOffload','Location','northeast','fontsize',9);%10                        
                    end
                    le.ItemTokenSize = ones(1,3)*10;%20
                    le.Color = backgounrdColor;
                elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 2
                    le = legend(h_line,'After correct','After error','Location','northeast','fontsize',10);
                    le.ItemTokenSize = ones(1,2)*20;
                    le.Color = backgounrdColor;                    
                elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 3
                    le = legend(h_line,'AfterCorrectMemory','AfterCorrectOffload','AfterCorrectNoChoice','Location','northeast','fontsize',10);
                    le.ItemTokenSize = ones(1,3)*20;
                    le.Color = backgounrdColor;                    
                elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 4
                    le = legend(h_line,'AfterErrorMemory','AfterErrorOffload','AfterErrorNoChoice','Location','northeast','fontsize',10);
                    le.ItemTokenSize = ones(1,3)*20;
                    le.Color = backgounrdColor;                    
                end
                
                set(gca,'color',backgounrdColor);
                xlim([1-xlim_padSize1 temp_ticks_decisionPeriod(end)+xlim_padSize2]);
                ylim([min_baselineDecisionPeriod,max_baselineDecisionPeriod]);
                
                set(gca,'xticklabel',[])
                set(gca,'xtick',[])
                
                xticks([temp_ticks_baselinePeriod temp_ticks_decisionPeriodA temp_ticks_decisionPeriod]);
                %xticklabels({'PreFixation','Fixation','','','Decision',''});
                xticklabels({'PreFixation','Fixation','','','ChoiceCue','','','Decision',''});
                
                text(temp_ticks_baselinePeriod(end),min_baselineDecisionPeriod,'/','fontsize',15);
                text(temp_ticks_decisionPeriodA(1),min_baselineDecisionPeriod,'/','fontsize',15);
                
                text(temp_ticks_decisionPeriodA(end),min_baselineDecisionPeriod,'/','fontsize',15);
                text(temp_ticks_decisionPeriod(1),min_baselineDecisionPeriod,'/','fontsize',15);
                
                set(gca, 'FontSize', 11)
                set(gca,'box','off');% 取消右、上边框
                %ylabel('dF / F', 'FontSize', 14, 'FontWeight', 'bold');
                %ylabel(sprintf('dF / F\ncellID=%d',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
                ylabel(sprintf('dF / F'), 'FontSize', 14, 'FontWeight', 'bold');                
                %title(sprintf('cell id (suite2p) = %d',temp_cellIndex_suite2p),'FontSize',14);
                
                a = 1;
                %             % Save fig
                %             currentFigName = sprintf('cellID%d',temp_cellIndex_suite2p);
                %             fileName_fig = [output_path_singleCell '\' currentFigName '.png'];%emf,pdf,png,jpg,tiff
                %             % to judge whether save figure or not
                %             if ifSave_fig == 1
                %                 %exportgraphics(fig1,fileName_fig,'BackgroundColor','none','ContentType','vector');
                %                 exportgraphics(fig1,fileName_fig,'Resolution',300);
                %                 close all
                %             end
                
                % elseif plot_lengthFlag ~= 1
                
                
                %                 %% subplot 6, roi anatomy info
                %                 nexttile([1,3]);
                %                 temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
                %
                %                 temp_roi_npix = double(temp_roi_stat.npix);
                %
                %                 temp_roi_xpix = double(temp_roi_stat.xpix);
                %                 temp_roi_ypix = double(temp_roi_stat.ypix);
                %
                %                 [temp_x_min,temp_x_max] = bounds(temp_roi_xpix);
                %                 [temp_y_min,temp_y_max] = bounds(temp_roi_ypix);
                %
                %
                %                 temp_diff = temp_x_max-temp_x_min;
                %                 if temp_diff < 5
                %                     temp_x_min = temp_x_min - 2;
                %                     temp_x_max = temp_x_max + 2;
                %                 end
                %                 temp_diff = temp_y_max-temp_y_min;
                %                 if temp_diff < 5
                %                     temp_y_min = temp_y_min - 2;
                %                     temp_y_max = temp_y_max + 2;
                %                 end
                %
                %                 a = 1;
                %
                %                 %fig1 = figure('Name',sprintf('fig1'),'NumberTitle','off');
                %                 %set(gcf,'Position',[35+0 35+0 200 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                %
                %                 scatter(temp_roi_xpix,temp_roi_ypix,60,'s','filled','MarkerFaceColor',[0.4660 0.6740 0.1880]);
                %
                %                 xlim([temp_x_min,temp_x_max]);
                %                 ylim([temp_y_min,temp_y_max]);
                %
                %                 axis equal;
                %                 set(gca, 'YDir','reverse');
                %                 title(sprintf('npix=%d',temp_roi_npix),'FontSize',12);
                
                
                
            elseif if_plot_location0_seq1 == 1
                xlim_padSize1 = 2;
                xlim_padSize2 = 12;
                
                
                multi_hslColor = zeros(numSeq(plot_lengthFlag), 3);
                multi_hslColor(:,2) = 0.5;
                multi_hslColor(:,3) = 0.5;
                for tempi=1:numSeq(plot_lengthFlag)
                    multi_hslColor(tempi,1) = tempi/numSeq(plot_lengthFlag);
                end
                multi_rgbColor = hsl2rgb(multi_hslColor);
                
                
                backgounrdColor = [1 1 1]*0.825;%0.875
                
                [min_period1,max_period1] = bounds(cellID_F_dff_lengthx_period1_seq_mean,'all');
                [min_period2,max_period2] = bounds(cellID_F_dff_lengthx_period2_seq_mean,'all');
                min_period12 = min(min_period1,min_period2) - max([cellID_F_dff_lengthx_period1_seq_sem cellID_F_dff_lengthx_period2_seq_sem],[],'all');
                max_period12 = max(max_period1,max_period2) + max([cellID_F_dff_lengthx_period1_seq_sem cellID_F_dff_lengthx_period2_seq_sem],[],'all');
                
                [min_period1_choiceMemory,max_period1_choiceMemory] = bounds(cellID_F_dff_lengthx_period1_seq_mean_choiceMemory,'all');
                [min_period2_choiceMemory,max_period2_choiceMemory] = bounds(cellID_F_dff_lengthx_period2_seq_mean_choiceMemory,'all');
                min_period12_choiceMemory = min(min_period1_choiceMemory,min_period2_choiceMemory) - max([cellID_F_dff_lengthx_period1_seq_sem_choiceMemory cellID_F_dff_lengthx_period2_seq_sem_choiceMemory],[],'all');
                max_period12_choiceMemory = max(max_period1_choiceMemory,max_period2_choiceMemory) + max([cellID_F_dff_lengthx_period1_seq_sem_choiceMemory cellID_F_dff_lengthx_period2_seq_sem_choiceMemory],[],'all');
                
                [min_period1_choiceOffload,max_period1_choiceOffload] = bounds(cellID_F_dff_lengthx_period1_seq_mean_choiceOffload,'all');
                [min_period2_choiceOffload,max_period2_choiceOffload] = bounds(cellID_F_dff_lengthx_period2_seq_mean_choiceOffload,'all');
                min_period12_choiceOffload = min(min_period1_choiceOffload,min_period2_choiceOffload) - max([cellID_F_dff_lengthx_period1_seq_sem_choiceOffload cellID_F_dff_lengthx_period2_seq_sem_choiceOffload],[],'all');
                max_period12_choiceOffload = max(max_period1_choiceOffload,max_period2_choiceOffload) + max([cellID_F_dff_lengthx_period1_seq_sem_choiceOffload cellID_F_dff_lengthx_period2_seq_sem_choiceOffload],[],'all');
                
                [min_period1_noChoice,max_period1_noChoice] = bounds(cellID_F_dff_lengthx_period1_seq_mean_noChoice,'all');
                [min_period2_noChoice,max_period2_noChoice] = bounds(cellID_F_dff_lengthx_period2_seq_mean_noChoice,'all');
                min_period12_noChoice = min(min_period1_noChoice,min_period2_noChoice) - max([cellID_F_dff_lengthx_period1_seq_sem_noChoice cellID_F_dff_lengthx_period2_seq_sem_noChoice],[],'all');
                max_period12_noChoice = max(max_period1_noChoice,max_period2_noChoice) + max([cellID_F_dff_lengthx_period1_seq_sem_noChoice cellID_F_dff_lengthx_period2_seq_sem_noChoice],[],'all');
                
                min_period12 = min([min_period12,min_period12_choiceMemory,min_period12_choiceOffload,min_period12_noChoice]);
                max_period12 = max([max_period12,max_period12_choiceMemory,max_period12_choiceOffload,max_period12_noChoice]);
                
                fig1 = figure('Name','Fig1','NumberTitle','off');
                set(gcf,'Position',[35+0 35+0 1630 1000]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                t = tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
                
                
                %% subplot 1, trial average in period12, memoryCorrect condition
                nexttile
                labels = [];
                h_line = [];
                for tempi=1:numSeq(plot_lengthFlag)
                    if size(cellID_F_dff_lengthx_period1_seq{tempi},1) <= 1
                        continue
                    end
                    temp_str = [];
                    for tempj=1:plot_lengthFlag
                        temp_str = [temp_str num2str(target_seqSet{plot_lengthFlag}{tempi}(tempj))]; %#ok<*AGROW>
                    end
                    temp_str = string(temp_str);
                    
                    labels = [labels temp_str];
                    
                    
                    x = temp_ticks_period1(1):temp_ticks_period1(end);
                    y = cellID_F_dff_lengthx_period1_seq_mean(tempi,:);
                    h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))];
                    hold on
                    y_sem = cellID_F_dff_lengthx_period1_seq_sem(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_period2A(1):temp_ticks_period2A(end);
                    y = cellID_F_dff_lengthx_period2A_seq_mean(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period2A_seq_sem(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_period2(1):temp_ticks_period2(end);
                    y = cellID_F_dff_lengthx_period2_seq_mean(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period2_seq_sem(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                end
                
                for tempi=1:length(lengthx_period1_interval)
                    plot(temp_ticks_period1(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period2A_interval)
                    plot(temp_ticks_period2A(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period2_interval)
                    plot(temp_ticks_period2(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:plot_lengthFlag
                    x = [lengthx_period1_interval(3+tempi-1) min_period12;...
                        lengthx_period1_interval(3+tempi-1) max_period12;...
                        lengthx_period1_interval(3+tempi-1)+6 max_period12;...
                        lengthx_period1_interval(3+tempi-1)+6 min_period12];
                    y = [1 2 3 4];
                    patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.1,'EdgeColor','none');
                    hold on
                end
                
                legend(h_line,labels,'Location','northeast','fontsize',6)
                
                set(gca,'color',backgounrdColor);
                %xlim([1-xlim_padSize1 lengthx_period1_interval(end)+periodSkipInterval+lengthx_period2_interval(end)+xlim_padSize2]);
                xlim([temp_ticks_period1(1)-xlim_padSize1 temp_ticks_period2(end)+xlim_padSize2]);
                ylim([min_period12,max_period12]);
                set(gca,'xticklabel',[])
                set(gca,'xtick',[])
                
                xticks([temp_ticks_period1 temp_ticks_period2A temp_ticks_period2 temp_ticks_period3]);
                if plot_lengthFlag == 1
                    %xticklabels({'PreFixation','Fixation','T1','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
                    %xticklabels({'PreFixation','Fixation','T1','Delay1',' ',' ','Decision','','','Submit',''}); 
                    xticklabels({'PreFixation','Fixation','T1','Delay1','','','ChoiceCue','','','Decision','','','Submit',''}); 
                elseif plot_lengthFlag == 2
                    xticklabels({'PreFixation','Fixation','T1','T2','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 3
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 4
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','T4','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});                    
                end                
                                
                text(temp_ticks_period1(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period2A(1),min_period12,'/','fontsize',15);
                
                text(temp_ticks_period2A(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period2(1),min_period12,'/','fontsize',15);
                
                set(gca, 'FontSize', 11)
                set(gca,'box','off');% 取消右、上边框
                %ylabel(sprintf('dF / F, memoryCorrect\ncellID=%d',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
                ylabel(sprintf('dF / F, cellID=%d\nmemoryCorrect',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
                
                %% subplot 2, trial average in period12, noChoice condition
                nexttile
                labels = [];
                h_line = [];
                for tempi=1:numSeq(plot_lengthFlag)
                    if size(cellID_F_dff_lengthx_period1_seq_noChoice{tempi},1) <= 1
                        continue
                    end
                    temp_str = [];
                    for tempj=1:plot_lengthFlag
                        temp_str = [temp_str num2str(target_seqSet{plot_lengthFlag}{tempi}(tempj))]; %#ok<*AGROW>
                    end
                    temp_str = string(temp_str);
                    
                    labels = [labels temp_str];
                    
                    
                    x = temp_ticks_period1(1):temp_ticks_period1(end);
                    y = cellID_F_dff_lengthx_period1_seq_mean_noChoice(tempi,:);
                    h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))];
                    hold on
                    y_sem = cellID_F_dff_lengthx_period1_seq_sem_noChoice(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_period2A(1):temp_ticks_period2A(end);
                    y = cellID_F_dff_lengthx_period2A_seq_mean_noChoice(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period2A_seq_sem_noChoice(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_period2(1):temp_ticks_period2(end);
                    y = cellID_F_dff_lengthx_period2_seq_mean_noChoice(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period2_seq_sem_noChoice(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                end
                
                for tempi=1:length(lengthx_period1_interval)
                    plot(temp_ticks_period1(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period2A_interval)
                    plot(temp_ticks_period2A(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period2_interval)
                    plot(temp_ticks_period2(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:plot_lengthFlag
                    x = [lengthx_period1_interval(3+tempi-1) min_period12;...
                        lengthx_period1_interval(3+tempi-1) max_period12;...
                        lengthx_period1_interval(3+tempi-1)+6 max_period12;...
                        lengthx_period1_interval(3+tempi-1)+6 min_period12];
                    y = [1 2 3 4];
                    patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.1,'EdgeColor','none');
                    hold on
                end
                
                legend(h_line,labels,'Location','northeast','fontsize',6)
                
                set(gca,'color',backgounrdColor);
                xlim([temp_ticks_period1(1)-xlim_padSize1 temp_ticks_period2(end)+xlim_padSize2]);
                ylim([min_period12,max_period12]);
                set(gca,'xticklabel',[])
                set(gca,'xtick',[])
                
                xticks([temp_ticks_period1 temp_ticks_period2A temp_ticks_period2 temp_ticks_period3]);
                if plot_lengthFlag == 1
                    %xticklabels({'PreFixation','Fixation','T1','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
                    %xticklabels({'PreFixation','Fixation','T1','Delay1',' ',' ','Decision','','','Submit',''}); 
                    xticklabels({'PreFixation','Fixation','T1','Delay1','','','ChoiceCue','','','Decision','','','Submit',''}); 
                elseif plot_lengthFlag == 2
                    xticklabels({'PreFixation','Fixation','T1','T2','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 3
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 4
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','T4','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});                    
                end                
                                
                text(temp_ticks_period1(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period2A(1),min_period12,'/','fontsize',15);
                
                text(temp_ticks_period2A(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period2(1),min_period12,'/','fontsize',15);
                
                set(gca, 'FontSize', 11)
                set(gca,'box','off');% 取消右、上边框
                ylabel(sprintf('dF / F, cellID=%d\nnoChoice',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
                
                %% subplot 3, trial average in period12, choiceMemory condition
                nexttile
                labels = [];
                h_line = [];
                for tempi=1:numSeq(plot_lengthFlag)
                    if size(cellID_F_dff_lengthx_period1_seq_choiceMemory{tempi},1) <= 1
                        continue
                    end
                    temp_str = [];
                    for tempj=1:plot_lengthFlag
                        temp_str = [temp_str num2str(target_seqSet{plot_lengthFlag}{tempi}(tempj))]; %#ok<*AGROW>
                    end
                    temp_str = string(temp_str);
                    
                    labels = [labels temp_str];
                    
                    
                    x = temp_ticks_period1(1):temp_ticks_period1(end);
                    y = cellID_F_dff_lengthx_period1_seq_mean_choiceMemory(tempi,:);
                    h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))];
                    hold on
                    y_sem = cellID_F_dff_lengthx_period1_seq_sem_choiceMemory(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_period2A(1):temp_ticks_period2A(end);
                    y = cellID_F_dff_lengthx_period2A_seq_mean_choiceMemory(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period2A_seq_sem_choiceMemory(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_period2(1):temp_ticks_period2(end);
                    y = cellID_F_dff_lengthx_period2_seq_mean_choiceMemory(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period2_seq_sem_choiceMemory(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                end
                
                for tempi=1:length(lengthx_period1_interval)
                    plot(temp_ticks_period1(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period2A_interval)
                    plot(temp_ticks_period2A(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period2_interval)
                    plot(temp_ticks_period2(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:plot_lengthFlag
                    x = [lengthx_period1_interval(3+tempi-1) min_period12;...
                        lengthx_period1_interval(3+tempi-1) max_period12;...
                        lengthx_period1_interval(3+tempi-1)+6 max_period12;...
                        lengthx_period1_interval(3+tempi-1)+6 min_period12];
                    y = [1 2 3 4];
                    patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.1,'EdgeColor','none');
                    hold on
                end
                
                legend(h_line,labels,'Location','northeast','fontsize',6)
                
                set(gca,'color',backgounrdColor);
                xlim([temp_ticks_period1(1)-xlim_padSize1 temp_ticks_period2(end)+xlim_padSize2]);
                ylim([min_period12,max_period12]);
                set(gca,'xticklabel',[])
                set(gca,'xtick',[])
                
                xticks([temp_ticks_period1 temp_ticks_period2A temp_ticks_period2 temp_ticks_period3]);
                if plot_lengthFlag == 1
                    %xticklabels({'PreFixation','Fixation','T1','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
                    %xticklabels({'PreFixation','Fixation','T1','Delay1',' ',' ','Decision','','','Submit',''}); 
                    xticklabels({'PreFixation','Fixation','T1','Delay1','','','ChoiceCue','','','Decision','','','Submit',''}); 
                elseif plot_lengthFlag == 2
                    xticklabels({'PreFixation','Fixation','T1','T2','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 3
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 4
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','T4','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});                                        
                end                
                                
                text(temp_ticks_period1(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period2A(1),min_period12,'/','fontsize',15);
                
                text(temp_ticks_period2A(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period2(1),min_period12,'/','fontsize',15);
                
                set(gca, 'FontSize', 11)
                set(gca,'box','off');% 取消右、上边框
                ylabel(sprintf('dF / F, cellID=%d\nchoiceMemory',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');
                
                %% subplot 4, trial average in period12, choiceOffload condition
                nexttile
                labels = [];
                h_line = [];
                for tempi=1:numSeq(plot_lengthFlag)
                    if size(cellID_F_dff_lengthx_period1_seq_choiceOffload{tempi},1) <= 1
                        continue
                    end
                    temp_str = [];
                    for tempj=1:plot_lengthFlag
                        temp_str = [temp_str num2str(target_seqSet{plot_lengthFlag}{tempi}(tempj))]; %#ok<*AGROW>
                    end
                    temp_str = string(temp_str);
                    
                    labels = [labels temp_str];
                    
                    
                    x = temp_ticks_period1(1):temp_ticks_period1(end);
                    y = cellID_F_dff_lengthx_period1_seq_mean_choiceOffload(tempi,:);
                    h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:))];
                    hold on
                    y_sem = cellID_F_dff_lengthx_period1_seq_sem_choiceOffload(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_period2A(1):temp_ticks_period2A(end);
                    y = cellID_F_dff_lengthx_period2A_seq_mean_choiceOffload(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period2A_seq_sem_choiceOffload(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                    x = temp_ticks_period2(1):temp_ticks_period2(end);
                    y = cellID_F_dff_lengthx_period2_seq_mean_choiceOffload(tempi,:);
                    plot(x,y,'-','LineWidth',1.5,'Color',multi_rgbColor(tempi,:));
                    hold on
                    y_sem = cellID_F_dff_lengthx_period2_seq_sem_choiceOffload(tempi,:);
                    patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',multi_rgbColor(tempi,:),'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                    hold on
                    
                end
                
                for tempi=1:length(lengthx_period1_interval)
                    plot(temp_ticks_period1(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period2A_interval)
                    plot(temp_ticks_period2A(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:length(lengthx_period2_interval)
                    plot(temp_ticks_period2(tempi)*[1 1],[min_period12 max_period12],...
                        '-','LineWidth',1,'Color',[0 0 0]);
                    hold on
                end
                
                for tempi=1:plot_lengthFlag
                    x = [lengthx_period1_interval(3+tempi-1) min_period12;...
                        lengthx_period1_interval(3+tempi-1) max_period12;...
                        lengthx_period1_interval(3+tempi-1)+6 max_period12;...
                        lengthx_period1_interval(3+tempi-1)+6 min_period12];
                    y = [1 2 3 4];
                    patch('Faces',y,'Vertices',x,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.1,'EdgeColor','none');
                    hold on
                end
                
                legend(h_line,labels,'Location','northeast','fontsize',6)
                
                set(gca,'color',backgounrdColor);
                xlim([temp_ticks_period1(1)-xlim_padSize1 temp_ticks_period2(end)+xlim_padSize2]);
                ylim([min_period12,max_period12]);
                set(gca,'xticklabel',[])
                set(gca,'xtick',[])
                
                xticks([temp_ticks_period1 temp_ticks_period2A temp_ticks_period2 temp_ticks_period3]);
                if plot_lengthFlag == 1
                    %xticklabels({'PreFixation','Fixation','T1','Delay1','LayoutOff',' ',' ','Decision','','','Submit',''});
                    %xticklabels({'PreFixation','Fixation','T1','Delay1',' ',' ','Decision','','','Submit',''}); 
                    xticklabels({'PreFixation','Fixation','T1','Delay1','','','ChoiceCue','','','Decision','','','Submit',''}); 
                elseif plot_lengthFlag == 2
                    xticklabels({'PreFixation','Fixation','T1','T2','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 3
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});
                elseif plot_lengthFlag == 4
                    xticklabels({'PreFixation','Fixation','T1','T2','T3','T4','Delay1','','','ChoiceCue','','','Decision','','','Submit',''});                                        
                end                
                                
                text(temp_ticks_period1(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period2A(1),min_period12,'/','fontsize',15);
                
                text(temp_ticks_period2A(end),min_period12,'/','fontsize',15);
                text(temp_ticks_period2(1),min_period12,'/','fontsize',15);
                
                set(gca, 'FontSize', 11)
                set(gca,'box','off');% 取消右、上边框
                ylabel(sprintf('dF / F, cellID=%d\nchoiceOffload',temp_cellIndex_suite2p), 'FontSize', 14, 'FontWeight', 'bold');                                               
                
            end
            
            % Save fig
            %currentFigName = sprintf('cellID%d',temp_cellIndex_suite2p);
            %currentFigName = sprintf('fig%d_cellID%d',tempIndex,temp_cellIndex_suite2p);
            if if_plot_location0_seq1 == 0
                currentFigName = sprintf('FOVID%s_length%d_fig%d_cellID%d',currentSession(18:22),plot_lengthFlag,tempIndex,temp_cellIndex_suite2p);
            elseif if_plot_location0_seq1 == 1
                currentFigName = sprintf('FOVID%s_seq_cellID%d_length%d_fig%d',currentSession(18:22),temp_cellIndex_suite2p,plot_lengthFlag,tempIndex);                
            end
            fileName_fig = [output_path_singleCell '\' currentFigName '.png'];%emf,pdf,png,jpg,tiff
            % to judge whether save figure or not
            if ifSave_fig == 1
                %exportgraphics(fig1,fileName_fig,'BackgroundColor','none','ContentType','vector');
                exportgraphics(fig1,fileName_fig,'Resolution',300);
                close all
            end
            
        end
    end
    
    a = 1;
    
    %% F_dff_lengthx_sample_raw_allTrial
    % F_dff_length1_period1_raw;
    % length1_period1_interval;
    temp0 = length1_period1_interval(3:4);
    temp0(end) = temp0(end)-1;
    temp_range = temp0(1):temp0(end);
    length1_sample_interval = temp0-temp0(1)+1;
    F_dff_length1_sample_raw = F_dff_length1_period1_raw(:,:,temp_range);
    
    % F_dff_length2_period1_raw;
    % length2_period1_interval;
    temp0 = length2_period1_interval(3:5);
    temp0(end) = temp0(end)-1;
    temp_range = temp0(1):temp0(end);
    length2_sample_interval = temp0-temp0(1)+1;
    F_dff_length2_sample_raw = F_dff_length2_period1_raw(:,:,temp_range);
    
    
    % F_dff_length3_period1_raw;
    % length3_period1_interval;
    temp0 = length3_period1_interval(3:6);
    temp0(end) = temp0(end)-1;
    temp_range = temp0(1):temp0(end);
    length3_sample_interval = temp0-temp0(1)+1;
    F_dff_length3_sample_raw = F_dff_length3_period1_raw(:,:,temp_range);
    
    F_dff_length1_sample_raw;
    F_dff_length2_sample_raw;
    F_dff_length3_sample_raw;
    
    F_dff_length1_sample_raw_allTrial = nan(size(F_dff_length1_sample_raw,1),trial_para.trial_count,size(F_dff_length1_sample_raw,3));
    F_dff_length2_sample_raw_allTrial = nan(size(F_dff_length2_sample_raw,1),trial_para.trial_count,size(F_dff_length2_sample_raw,3));
    F_dff_length3_sample_raw_allTrial = nan(size(F_dff_length3_sample_raw,1),trial_para.trial_count,size(F_dff_length3_sample_raw,3));
    
    F_dff_length1_sample_raw_allTrial(:,trialIndex_length1_bool,:) = F_dff_length1_sample_raw;
    F_dff_length2_sample_raw_allTrial(:,trialIndex_length2_bool,:) = F_dff_length2_sample_raw;
    F_dff_length3_sample_raw_allTrial(:,trialIndex_length3_bool,:) = F_dff_length3_sample_raw;
    
    
    a1 = squeeze(F_dff_length1_sample_raw_allTrial(:,1,:));
    a2 = squeeze(F_dff_length1_sample_raw(:,1,:));
    a11 = squeeze(F_dff_length1_sample_raw_allTrial(:,3,:));
    
    F_dff_length1_sample = F_dff_length1_sample_raw_allTrial;
    F_dff_length2_sample = F_dff_length2_sample_raw_allTrial;
    F_dff_length3_sample = F_dff_length3_sample_raw_allTrial;
    length1_sample_interval;
    length2_sample_interval;
    length3_sample_interval;

    %% decodingDataSimplified
    decodingDataSimplified = struct;
    decodingData;
    
    decodingDataSimplified.F_dff_baselinePeriod = F_dff_baselinePeriod;
    decodingDataSimplified.baselinePeriod_interval = baselinePeriod_interval;        
    decodingDataSimplified.F_dff_length1_sample = F_dff_length1_sample;
    decodingDataSimplified.length1_sample_interval = length1_sample_interval;    
    decodingDataSimplified.F_dff_length2_sample = F_dff_length2_sample;
    decodingDataSimplified.length2_sample_interval = length2_sample_interval;    
    decodingDataSimplified.F_dff_length3_sample = F_dff_length3_sample;
    decodingDataSimplified.length3_sample_interval = length3_sample_interval; 
    decodingDataSimplified.F_dff_decisionPeriodB = F_dff_decisionPeriodB;
    decodingDataSimplified.decisionPeriodB_interval = decisionPeriodB_interval;    
    decodingDataSimplified.F_dff_decisionPeriodA = F_dff_decisionPeriodA;
    decodingDataSimplified.decisionPeriodA_interval = decisionPeriodA_interval;
    decodingDataSimplified.F_dff_decisionPeriod = F_dff_decisionPeriod;
    decodingDataSimplified.decisionPeriod_interval = decisionPeriod_interval;
    decodingDataSimplified.seqIndex = seqIndex;
    decodingDataSimplified.trial_para_currentSequence = trial_para.currentSequence;
    decodingDataSimplified.trialIndex_bool_memoryCorrect = trialIndex_bool_memoryCorrect;
    decodingDataSimplified.target_seqSet_inOne = target_seqSet_inOne;
    decodingDataSimplified.target_seqSet = target_seqSet;
    decodingDataSimplified.trial_para_ifSelectOffloading = trial_para.ifSelectOffloading;
    decodingDataSimplified.trial_para_isFilled = trial_para.isFilled;    
    decodingDataSimplified.trial_para_choiceCondition_flag = trial_para.choiceCondition_flag;
    decodingDataSimplified.trial_para_isCorrect = trial_para.isCorrect;

    
    if if_save_decoding_and_PTB_data == 1
        save([output_path,'\','decodingData.mat'],'decodingData');
        save([output_path,'\','trial_para.mat'],'trial_para');
        
        save([output_path,'\','decodingDataSimplified.mat'],'decodingDataSimplified');        
    end
    
    % if if_save_glmData == 1
    %     fprintf('length = [%d, %d, %d].\n',size(sequence_length1_memoryCorrect_onehot,1),size(sequence_length2_memoryCorrect_onehot,1),size(sequence_length3_memoryCorrect_onehot,1));
    %     sequence_length1_memoryCorrect_onehot;
    %     sequence_length2_memoryCorrect_onehot;
    %     sequence_length3_memoryCorrect_onehot;
    %
    %     F_dff_length1_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memoryCorrect(1,:));
    %     F_dff_length1_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memoryCorrect(1,:));
    %     F_dff_length1_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memoryCorrect(1,:));
    %
    %     F_dff_length2_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memoryCorrect(2,:));
    %     F_dff_length2_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memoryCorrect(2,:));
    %     F_dff_length2_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memoryCorrect(2,:));
    %
    %     F_dff_length3_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memoryCorrect(3,:));
    %     F_dff_length3_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memoryCorrect(3,:));
    %     F_dff_length3_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memoryCorrect(3,:));
    %
    %     glmData = struct;
    %     glmData.sequence_length1_memoryCorrect_onehot = sequence_length1_memoryCorrect_onehot;
    %     glmData.sequence_length2_memoryCorrect_onehot = sequence_length2_memoryCorrect_onehot;
    %     glmData.sequence_length3_memoryCorrect_onehot = sequence_length3_memoryCorrect_onehot;
    %     glmData.F_dff_length1_delay1Bin = F_dff_length1_delay1Bin;
    %     glmData.F_dff_length1_delay2Bin = F_dff_length1_delay2Bin;
    %     glmData.F_dff_length1_baselineBin = F_dff_length1_baselineBin;
    %     glmData.F_dff_length2_delay1Bin = F_dff_length2_delay1Bin;
    %     glmData.F_dff_length2_delay2Bin = F_dff_length2_delay2Bin;
    %     glmData.F_dff_length2_baselineBin = F_dff_length2_baselineBin;
    %     glmData.F_dff_length3_delay1Bin = F_dff_length3_delay1Bin;
    %     glmData.F_dff_length3_delay2Bin = F_dff_length3_delay2Bin;
    %     glmData.F_dff_length3_baselineBin = F_dff_length3_baselineBin;
    %
    %     %tempName = ['glm',currentSession(14:22),'.mat']; %#ok<*UNRCH>
    %     tempName = ['glm',currentSession(14:22),'_memoryCorrect.mat']; %#ok<*UNRCH>
    %     %save(['C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts','\data\',tempName],'glmData');
    %     save(['C:\ASDROOT\STUDY\temp','\data\',tempName],'glmData');
    % end
    
    
    
    % sequence_length1_onehot
    sequence_length1_onehot = false(size(sequence_length1,2),numFrames);
    for tempi=1:size(sequence_length1,2)
        sequence_length1_onehot(tempi,sequence_length1(tempi)) = true;
    end
    % sequence_length2_onehot
    sequence_length2_onehot = false(size(sequence_length2,2),numFrames);
    for tempi=1:size(sequence_length2,2)
        sequence_length2_onehot(tempi,sequence_length2(:,tempi)) = true;
    end
    % sequence_length3_onehot
    sequence_length3_onehot = false(size(sequence_length3,2),numFrames);
    for tempi=1:size(sequence_length3,2)
        sequence_length3_onehot(tempi,sequence_length3(:,tempi)) = true;
    end
    % sequence_length4_onehot
    sequence_length4_onehot = false(size(sequence_length4,2),numFrames);
    for tempi=1:size(sequence_length4,2)
        sequence_length4_onehot(tempi,sequence_length4(:,tempi)) = true;
    end
    responseSeq_length1_onehot;
    responseSeq_length2_onehot;
    responseSeq_length3_onehot;
    responseSeq_length4_onehot;
    
    
    
    trialIndex_length1;
    ifSelectOffloading_bool;
    
    tempBoolIndex = ifSelectOffloading_bool(trialIndex_length1)==false;
    sequence_length1_onehot_memory = sequence_length1_onehot(tempBoolIndex,:);
    responseSeq_length1_onehot_memory = responseSeq_length1_onehot(tempBoolIndex,:);
    
    tempBoolIndex = ifSelectOffloading_bool(trialIndex_length2)==false;
    sequence_length2_onehot_memory = sequence_length2_onehot(tempBoolIndex,:);
    responseSeq_length2_onehot_memory = responseSeq_length2_onehot(tempBoolIndex,:);
    
    tempBoolIndex = ifSelectOffloading_bool(trialIndex_length3)==false;
    sequence_length3_onehot_memory = sequence_length3_onehot(tempBoolIndex,:);
    responseSeq_length3_onehot_memory = responseSeq_length3_onehot(tempBoolIndex,:);
    
    tempBoolIndex = ifSelectOffloading_bool(trialIndex_length4)==false;
    sequence_length4_onehot_memory = sequence_length4_onehot(tempBoolIndex,:);
    responseSeq_length4_onehot_memory = responseSeq_length4_onehot(tempBoolIndex,:);
    
    
    if if_save_glmData == 1
        
        fprintf('length = [%d, %d, %d, %d].\n',size(sequence_length1_onehot_memory,1),size(sequence_length2_onehot_memory,1),...
            size(sequence_length3_onehot_memory,1),size(sequence_length4_onehot_memory,1));
        
        F_dff_length1_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memory(1,:));
        F_dff_length1_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memory(1,:));
        F_dff_length1_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memory(1,:));
        
        F_dff_length2_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memory(2,:));
        F_dff_length2_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memory(2,:));
        F_dff_length2_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memory(2,:));
        
        F_dff_length3_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memory(3,:));
        F_dff_length3_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memory(3,:));
        F_dff_length3_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memory(3,:));
        
        F_dff_length4_delay1Bin = F_dff_delay1Bin(:,trialIndex_lengthx_bool_memory(4,:));
        F_dff_length4_delay2Bin = F_dff_delay2Bin(:,trialIndex_lengthx_bool_memory(4,:));
        F_dff_length4_baselineBin = F_dff_baselineBin(:,trialIndex_lengthx_bool_memory(4,:));
        
        glmData = struct;
        glmData.sequence_length1_onehot = sequence_length1_onehot_memory;
        glmData.sequence_length2_onehot = sequence_length2_onehot_memory;
        glmData.sequence_length3_onehot = sequence_length3_onehot_memory;
        glmData.sequence_length4_onehot = sequence_length4_onehot_memory;
        glmData.responseSeq_length1_onehot = responseSeq_length1_onehot_memory;
        glmData.responseSeq_length2_onehot = responseSeq_length2_onehot_memory;
        glmData.responseSeq_length3_onehot = responseSeq_length3_onehot_memory;
        glmData.responseSeq_length4_onehot = responseSeq_length4_onehot_memory;
        glmData.F_dff_length1_delay1Bin = F_dff_length1_delay1Bin;
        glmData.F_dff_length1_delay2Bin = F_dff_length1_delay2Bin;
        glmData.F_dff_length1_baselineBin = F_dff_length1_baselineBin;
        glmData.F_dff_length2_delay1Bin = F_dff_length2_delay1Bin;
        glmData.F_dff_length2_delay2Bin = F_dff_length2_delay2Bin;
        glmData.F_dff_length2_baselineBin = F_dff_length2_baselineBin;
        glmData.F_dff_length3_delay1Bin = F_dff_length3_delay1Bin;
        glmData.F_dff_length3_delay2Bin = F_dff_length3_delay2Bin;
        glmData.F_dff_length3_baselineBin = F_dff_length3_baselineBin;
        glmData.F_dff_length4_delay1Bin = F_dff_length4_delay1Bin;
        glmData.F_dff_length4_delay2Bin = F_dff_length4_delay2Bin;
        glmData.F_dff_length4_baselineBin = F_dff_length4_baselineBin;
        
        glmData.cellIndex = cellIndex;
        glmData.cellIndex_suite2p = cellIndex_suite2p;
        glmData.selectiveCellBoolIndex_length1 = selectiveCellBoolIndex_length1;
        glmData.selectiveCellBoolIndex_length1_delay1Bin = selectiveCellBoolIndex_length1_delay1Bin;
        glmData.selectiveCellBoolIndex_length1_delay2Bin = selectiveCellBoolIndex_length1_delay2Bin;
        glmData.selectiveCellBoolIndex_length1_T = selectiveCellBoolIndex_length1_T;
        glmData.selectiveCellBoolIndex_length2 = selectiveCellBoolIndex_length2;
        glmData.selectiveCellBoolIndex_length2_delay1Bin = selectiveCellBoolIndex_length2_delay1Bin;
        glmData.selectiveCellBoolIndex_length2_delay2Bin = selectiveCellBoolIndex_length2_delay2Bin;        
        glmData.selectiveCellBoolIndex_length2_T = selectiveCellBoolIndex_length2_T;
        glmData.selectiveCellBoolIndex_length3 = selectiveCellBoolIndex_length3;
        glmData.selectiveCellBoolIndex_length3_delay1Bin = selectiveCellBoolIndex_length3_delay1Bin;
        glmData.selectiveCellBoolIndex_length3_delay2Bin = selectiveCellBoolIndex_length3_delay2Bin;        
        glmData.selectiveCellBoolIndex_length3_T = selectiveCellBoolIndex_length3_T;
        glmData.selectiveCellBoolIndex_length4 = selectiveCellBoolIndex_length4;
        glmData.selectiveCellBoolIndex_length4_delay1Bin = selectiveCellBoolIndex_length4_delay1Bin;
        glmData.selectiveCellBoolIndex_length4_delay2Bin = selectiveCellBoolIndex_length4_delay2Bin;        
        glmData.selectiveCellBoolIndex_length4_T = selectiveCellBoolIndex_length4_T;
        
        glmData.rProb_glm_output_all = rProb_glm_output_all;
        
        %tempName = ['glm',currentSession(14:22),'.mat']; %#ok<*UNRCH>
        tempName = ['glm',currentSession(14:22),'_allTrial.mat']; %#ok<*UNRCH>
        %save(['C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts','\data\',tempName],'glmData');
        save(['C:\ASDROOT\STUDY\temp','\data\',tempName],'glmData');
    end
    
    if if_profile == 1
        profile viewer
    end
    
    fprintf('Time is %.1f secs (roiNum = %d).\n',toc(t0),size(F_dff,1));
    
end

fprintf('Time of the stepB_analysis is %.1f secs.\n',toc(time0));
fprintf('num_FOV=%d.\n',num_FOV);

%% End
