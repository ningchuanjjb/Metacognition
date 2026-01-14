%% Initialization
% clear
close all
% home;% To scroll down in command window

% plot_lengthFlag = 1;
% valid_lengthRange_forPlot = [1 2 3 4];
% valid_lengthRange_forPlot = 4:4;
valid_lengthRange_forPlot = 1:1;
% valid_lengthRange_forPlot = 1:3;

if_load = 0;
preLoad_trigger = 0;
preLoad_status = 1;

if_plotSingleCell = 1;
ifSave_fig = 0;

if_singleFOV = 1;

if_0normal_1length123_2length1234 = 0;

order_glm = 0;% 0 means whole seq

if_compute_selectiveCellNum_decisionBin_lengthxSeq_MON = 0;


if_plot_mappingTwoSessions = 0;

if_plot_svm_beta = 0;

if_plot_location0_seq1 = 0;

if_substractBaseLine = 0;%1

if_plotSelectiveCell = 1;
%selectiveCellIndex_suite2p_lengthx

if_sortRasterPlot_delay1 = 0;
if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 = 1;%2
if_plotSelectiveCell_T0_delay1_MNO2 = 1;

plotRoiNum = 1;
temp_id_suite2p = 0;

temp_id = temp_id_suite2p + 1;

windowSize = 5;%3-->5

compensateImageDelay = 0;

if_profile = 0;

if if_profile == 1
    profile on
end

currentSession_multi = string;


currentSession_multi = [currentSession_multi; '113Recording_20230426A_Ding_Site16'];
currentSession_multi = [currentSession_multi; '113Recording_20230427A_Ding_Site16_sameFOV0426'];
currentSession_multi = [currentSession_multi; '113Recording_20230502A_Ding_Site13'];
currentSession_multi = [currentSession_multi; '113Recording_20230503A_Ding_Site13_sameFOV0502'];
currentSession_multi = [currentSession_multi; '113Recording_20230504A_Ding_Site02'];
currentSession_multi = [currentSession_multi; '113Recording_20230508A_Ding_Site02_sameFOV0509'];
currentSession_multi = [currentSession_multi; '113Recording_20230509A_Ding_Site02'];% 660000 frames, easy to crash

currentSession_multi = [currentSession_multi; '113Recording_20230510A_Ding_Site05_sameFOV0511'];
currentSession_multi = [currentSession_multi; '113Recording_20230510B_Ding_Site05_sameFOV0511'];
currentSession_multi = [currentSession_multi; '113Recording_20230511A_Ding_Site05'];
currentSession_multi = [currentSession_multi; '113Recording_20230512A_Ding_Site09'];
currentSession_multi = [currentSession_multi; '113Recording_20230513A_Ding_Site09_sameFOV0512'];

currentSession_multi = [currentSession_multi; '113Recording_20230515A_Ding_Site24_sameFOV0516'];
currentSession_multi = [currentSession_multi; '113Recording_20230516A_Ding_Site24'];
currentSession_multi = [currentSession_multi; '113Recording_20230517A_Ding_Site16B'];
currentSession_multi = [currentSession_multi; '113Recording_20230522A_Ding_Site05B'];
currentSession_multi = [currentSession_multi; '113Recording_20230525A_Ding_Site09B'];
currentSession_multi = [currentSession_multi; '113Recording_20230526A_Ding_Site09B_sameFOV0525'];

currentSession_multi = [currentSession_multi; '113Recording_20230527A_Ding_Site02B'];
currentSession_multi = [currentSession_multi; '113Recording_20230528A_Ding_Site02B_sameFOV0527'];
currentSession_multi = [currentSession_multi; '113Recording_20230530A_Ding_Site05C'];
currentSession_multi = [currentSession_multi; '113Recording_20230531A_Ding_Site05C_sameFOV0530'];
currentSession_multi = [currentSession_multi; '113Recording_20230601A_Ding_Site13B'];
currentSession_multi = [currentSession_multi; '113Recording_20230602A_Ding_Site13B_sameFOV0601'];

currentSession_multi = [currentSession_multi; '113Recording_20230604A_Ding_Site07'];
currentSession_multi = [currentSession_multi; '113Recording_20230605A_Ding_Site07_sameFOV0604'];
currentSession_multi = [currentSession_multi; '113Recording_20230612A_Ding_Site14'];
currentSession_multi = [currentSession_multi; '113Recording_20230614A_Ding_Site05D'];
currentSession_multi = [currentSession_multi; '113Recording_20230615A_Ding_Site05D_sameFOV0614'];
currentSession_multi = [currentSession_multi; '113Recording_20230619A_Ding_Site02C'];
currentSession_multi = [currentSession_multi; '113Recording_20230620A_Ding_Site05E'];

currentSession_multi(1) = [];
num_FOV = length(currentSession_multi);

glm_r2_lengthx_multiFOV = glm_r2_lengthx_delay1Bin_multiFOV;

tempSelectiveCellBoolIndex_length1 = tempBoolIndex1;
tempSelectiveCellBoolIndex_length2 = tempBoolIndex2;
tempSelectiveCellBoolIndex_length3 = tempBoolIndex3;
tempSelectiveCellBoolIndex_length4 = tempBoolIndex4;

tempSelectiveCellBoolIndex_length12_and = tempSelectiveCellBoolIndex_length1 & tempSelectiveCellBoolIndex_length2;
tempSelectiveCellBoolIndex_length13_and = tempSelectiveCellBoolIndex_length1 & tempSelectiveCellBoolIndex_length3;
tempSelectiveCellBoolIndex_length23_and = tempSelectiveCellBoolIndex_length2 & tempSelectiveCellBoolIndex_length3;
tempSelectiveCellBoolIndex_length123_and = tempSelectiveCellBoolIndex_length1 & tempSelectiveCellBoolIndex_length2 & tempSelectiveCellBoolIndex_length3;
tempSelectiveCellBoolIndex_length1234_and = tempSelectiveCellBoolIndex_length1 & tempSelectiveCellBoolIndex_length2 & tempSelectiveCellBoolIndex_length3 & tempSelectiveCellBoolIndex_length4;


% tempSelectiveCellBoolIndex_length123_and = tempSelectiveCellBoolIndex_length23_and & (~tempSelectiveCellBoolIndex_length123_and);

if if_plot_mappingTwoSessions == 1
    tempMappingCellIndex_suite2p_withCorr_sorted = decodingDataSimplified_AB.extraForMerged.tempMappingCellIndex_suite2p_withCorr_sorted;
    if_load = 0;
    %if_preLoad = 1;
    
    temp_sortedStartIndex = 196;    
end



% if if_preLoad == 1
if preLoad_trigger == 1
    t0_preLoad = tic;
    
    markerParse_trialLevel_preLoad = cell(1,num_FOV);
    basic_para_preLoad = cell(1,num_FOV);
    trial_para_preLoad = cell(1,num_FOV);
    iscell_preLoad = cell(1,num_FOV);
    roi_stats_preLoad = cell(1,num_FOV);
    F_dff_preLoad = cell(1,num_FOV);
    F0_preLoad = cell(1,num_FOV);
    template_preLoad = cell(1,num_FOV);
    
    for tempSessionIndex=1:num_FOV
        currentSession = currentSession_multi{tempSessionIndex};
        
        %searchName_gAcc = 'from23-04-26to23-06-12_D_gAcc_1';
        %searchName_rProb = 'from23-04-26to23-06-12_D_offloadingProb_1';
        searchName_gAcc = 'from23-04-26to23-06-20_D_gAcc_1';
        searchName_rProb = 'from23-04-26to23-06-20_D_offloadingProb_1';
        
        targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
        cd(targetPATH)
        
        output_shortPath = 'D:\twoPhotonData_motionCorrected';
        temp_currentSession_path = [output_shortPath '\' currentSession];
        temp_if_max0_min1 = 0;
        output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
        path_plane = [output_path,'\plane0'];
        
        rawData_path = 'R:\twoPhotonRawData\ToBeProcessed';
        currentSession_path = [rawData_path '\' currentSession];
        PTB_name = '.mat';
        PTB_fullName = autoGetMarkerFileName(['2*',PTB_name], currentSession_path);
        
        
        temp_load = load([output_path,'\','markerParse_trialLevel.mat'],'markerParse_trialLevel');
        markerParse_trialLevel_preLoad{tempSessionIndex} = temp_load.markerParse_trialLevel;
        
        basic_para = [];
        trial_para = [];
        temp_load = load(PTB_fullName,'basic_para','trial_para');
        basic_para_preLoad{tempSessionIndex} = temp_load.basic_para;
        trial_para_preLoad{tempSessionIndex} = temp_load.trial_para;
        
        path_behavior = [output_shortPath '\behavior'];
        
        % Load other processed results
        load_gAcc = loadMat_single(searchName_gAcc, path_behavior);
        gAcc_noChoice_collapsed_inOne = load_gAcc.gAcc_noChoice_collapsed_inOne;
        
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
        
        temp_load = load(fullFileName_Fall,'F_dff_raw');
        F_dff_raw = temp_load.F_dff_raw;
        
        temp_load = load(fullFileName_Fall,'F0_raw');
        F0_raw = temp_load.F0_raw;
        
        iscell = readNPY(fullFileName_iscell);
        iscell_preLoad{tempSessionIndex} = iscell;
        
        F_dff = F_dff_raw(iscell(:,1)==1, :);
        F0 = F0_raw(iscell(:,1)==1, :);
        clear F_dff_raw
        clear F0_raw
        
        F_dff = smoothdata(F_dff,2,'gaussian',windowSize);
        F_dff(:,1:end-compensateImageDelay) = F_dff(:,1+compensateImageDelay:end); % to compensate imaging delay
        F_dff_preLoad{tempSessionIndex} = F_dff;
        
        F0 = smoothdata(F0,2,'gaussian',windowSize);
        F0(:,1:end-compensateImageDelay) = F0(:,1+compensateImageDelay:end); % to compensate imaging delay
        F0_preLoad{tempSessionIndex} = F0;
        
        
        s = load(fullFileName_Fall,'stat');
        roi_stats_raw = s.stat;
        temp_cellIndex = find(iscell(:,1)==1);
        roi_stats = roi_stats_raw(temp_cellIndex);
        roi_stats_preLoad{tempSessionIndex} = roi_stats;
        
        
        temp_if_max0_min1 = 0;
        template_path = autoGetFileName_general('maxProjection*.tif', output_path,temp_if_max0_min1);
        template = double(loadtiff(template_path));
        template_preLoad{tempSessionIndex} = template;
    end
    preLoad_status = 1;
    fprintf('Time of pre-load is %.1f secs.\n',toc(t0_preLoad));
end




time0 = tic;
for tempSessionIndex=1:num_FOV
    % for tempSessionIndex=num_FOV:-1:1
    
    t0 = tic;
    
    if if_singleFOV == 1
        tempSessionIndex = 20; %#ok<*FXSET>
    end
    
    currentSession = currentSession_multi{tempSessionIndex};
    fprintf('currentSession = %s.\n',currentSession);
    
    FOVName = selevtivity_multiFOV.FOVName_multiFOV(tempSessionIndex,:);
    tempCellIndex_currentFOV = find(selevtivity_multiFOV.cell_FOVIndex_multiFOV == tempSessionIndex);
    
    %glm_r2_currentFOV = glm_r2_lengthx_multiFOV(tempCellIndex_currentFOV,1:3);
    glm_r2_currentFOV = glm_r2_lengthx_multiFOV(tempCellIndex_currentFOV,1:4);
    
    if if_0normal_1length123_2length1234 == 0
        tempSelectiveCellBoolIndex_length1_currentFOV = tempSelectiveCellBoolIndex_length1(tempCellIndex_currentFOV);
        tempSelectiveCellBoolIndex_length2_currentFOV = tempSelectiveCellBoolIndex_length2(tempCellIndex_currentFOV);
        tempSelectiveCellBoolIndex_length3_currentFOV = tempSelectiveCellBoolIndex_length3(tempCellIndex_currentFOV);
        tempSelectiveCellBoolIndex_length4_currentFOV = tempSelectiveCellBoolIndex_length4(tempCellIndex_currentFOV);
    elseif if_0normal_1length123_2length1234 == 1
        tempSelectiveCellBoolIndex_length1_currentFOV = tempSelectiveCellBoolIndex_length123_and(tempCellIndex_currentFOV);
        tempSelectiveCellBoolIndex_length2_currentFOV = tempSelectiveCellBoolIndex_length123_and(tempCellIndex_currentFOV);
        tempSelectiveCellBoolIndex_length3_currentFOV = tempSelectiveCellBoolIndex_length123_and(tempCellIndex_currentFOV);
        tempSelectiveCellBoolIndex_length4_currentFOV = tempSelectiveCellBoolIndex_length123_and(tempCellIndex_currentFOV);
    elseif if_0normal_1length123_2length1234 == 2
        tempSelectiveCellBoolIndex_length1_currentFOV = tempSelectiveCellBoolIndex_length1234_and(tempCellIndex_currentFOV);
        tempSelectiveCellBoolIndex_length2_currentFOV = tempSelectiveCellBoolIndex_length1234_and(tempCellIndex_currentFOV);
        tempSelectiveCellBoolIndex_length3_currentFOV = tempSelectiveCellBoolIndex_length1234_and(tempCellIndex_currentFOV);
        tempSelectiveCellBoolIndex_length4_currentFOV = tempSelectiveCellBoolIndex_length1234_and(tempCellIndex_currentFOV);
    end
    
    a = 1;
    temp_sum = 0;
    for tempi=1:length(valid_lengthRange_forPlot)
        temp_length = valid_lengthRange_forPlot(tempi);
        if temp_length == 1
            temp_sum = temp_sum + sum(tempSelectiveCellBoolIndex_length1_currentFOV);
        elseif temp_length == 2
            temp_sum = temp_sum + sum(tempSelectiveCellBoolIndex_length2_currentFOV);
        elseif temp_length == 3
            temp_sum = temp_sum + sum(tempSelectiveCellBoolIndex_length3_currentFOV);
        elseif temp_length == 4
            temp_sum = temp_sum + sum(tempSelectiveCellBoolIndex_length4_currentFOV);
        end
    end
    if temp_sum == 0
        fprintf('No valid roi to be ploted in current FOV, skip it.\n');
        continue
    end
    
    
    
    %searchName_gAcc = 'from23-04-26to23-06-12_D_gAcc_1';
    %searchName_rProb = 'from23-04-26to23-06-12_D_offloadingProb_1';
    searchName_gAcc = 'from23-04-26to23-06-20_D_gAcc_1';
    searchName_rProb = 'from23-04-26to23-06-20_D_offloadingProb_1';
    
    targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
    cd(targetPATH)
    
    output_shortPath = 'D:\twoPhotonData_motionCorrected';
    temp_currentSession_path = [output_shortPath '\' currentSession];
    temp_if_max0_min1 = 0;
    output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
    path_plane = [output_path,'\plane0'];
    
    output_path_singleCell = [output_shortPath,'\singleCell'];
    if ~exist(output_path_singleCell, 'dir')
        mkdir(output_path_singleCell);
    end
    
    
    rawData_path = 'R:\twoPhotonRawData\ToBeProcessed';
    currentSession_path = [rawData_path '\' currentSession];
    PTB_name = '.mat';
    PTB_fullName = autoGetMarkerFileName(['2*',PTB_name], currentSession_path);
    
    if if_load == 1
        load([output_path,'\','markerParse_trialLevel.mat'],'markerParse_trialLevel');
        
        basic_para = [];
        trial_para = [];
        load(PTB_fullName,'basic_para','trial_para');
        a = 1;
        
        path_behavior = [output_shortPath '\behavior'];
        
        % Load other processed results
        load_gAcc = loadMat_single(searchName_gAcc, path_behavior);
        gAcc_noChoice_collapsed_inOne = load_gAcc.gAcc_noChoice_collapsed_inOne;
        
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
        
        load(fullFileName_Fall,'F_dff_raw');
        load(fullFileName_Fall,'F0_raw');
        iscell = readNPY(fullFileName_iscell);
        
        F_dff = F_dff_raw(iscell(:,1)==1, :);
        F0 = F0_raw(iscell(:,1)==1, :);
        clear F_dff_raw
        clear F0_raw
        
        s = load(fullFileName_Fall,'stat');
        roi_stats_raw = s.stat;
        temp_cellIndex = find(iscell(:,1)==1);
        roi_stats = roi_stats_raw(temp_cellIndex);
        
        x = F_dff;
        F_dff = smoothdata(x,2,'gaussian',windowSize);
        F_dff(:,1:end-compensateImageDelay) = F_dff(:,1+compensateImageDelay:end); % to compensate imaging delay
        
        F0 = smoothdata(F0,2,'gaussian',windowSize);
        F0(:,1:end-compensateImageDelay) = F0(:,1+compensateImageDelay:end); % to compensate imaging delay
        
        temp_if_max0_min1 = 0;
        template_path = autoGetFileName_general('maxProjection*.tif', output_path,temp_if_max0_min1);
        template = double(loadtiff(template_path));
    end
    
    if preLoad_status == 1
        markerParse_trialLevel = markerParse_trialLevel_preLoad{tempSessionIndex};
        basic_para = basic_para_preLoad{tempSessionIndex};
        trial_para = trial_para_preLoad{tempSessionIndex};
        iscell = iscell_preLoad{tempSessionIndex};
        roi_stats = roi_stats_preLoad{tempSessionIndex};
        F_dff = F_dff_preLoad{tempSessionIndex};
        F0 = F0_preLoad{tempSessionIndex};
        template = template_preLoad{tempSessionIndex};
    end
    
    
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
    end
    baselinePeriod_interval = baselinePeriod_frameIndex(1,:)-baselinePeriod_frameIndex(1,1)+1;
    
    [baselinePeriod_frameRangeMin,~] = bounds(baselinePeriod_frameIndex,2);
    temp_frameRange = repmat(baselinePeriod_frameRangeMin,1,baselinePeriod_interval(end)-baselinePeriod_interval(1)+1);
    baselinePeriod_frameRange = temp_frameRange + ((baselinePeriod_interval(1):baselinePeriod_interval(end))-1);
    
    temp_dff = reshape(F_dff(:,baselinePeriod_frameRange'),[size(F_dff,1),size(baselinePeriod_frameRange')]);
    F_dff_baselinePeriod = permute(temp_dff,[1 3 2]);
    
    F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,baselinePeriod_interval(1):baselinePeriod_interval(3)),3);
    
    F_dff_baselineInit = mean(F_dff_baselinePeriod(:,:,1),2);
    if if_substractBaseLine == 1
        F0_baseline = F_dff_baselineInit.*F0;
        
        F = (F_dff.*F0)+F0;
        F_dff = (F-F0-F0_baseline)./(F0+F0_baseline);
        
        temp_dff = reshape(F_dff(:,baselinePeriod_frameRange'),[size(F_dff,1),size(baselinePeriod_frameRange')]);
        F_dff_baselinePeriod = permute(temp_dff,[1 3 2]);
        F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,baselinePeriod_interval(1):baselinePeriod_interval(3)),3);
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
    
    %% Extract decision
    % Decision period : SELECTING AND DELAY 2 ON - 1500 ms--> SELECTING AND DELAY 2 ON --> SELECTING AND DELAY 2 ON + 1500 ms
    %F_dff
    %F_dff_decisionPeriod_choiceMemory
    %F_dff_decisionPeriod_choiceOffload
    
    
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
    
    
    for plot_lengthFlag=valid_lengthRange_forPlot
        %for plot_lengthFlag=1:3
        
        
        %% plot_lengthFlag
        if plot_lengthFlag == 1
            if if_plotSelectiveCell_T0_delay1_MNO2 == 0
                % abandon now
            elseif if_plotSelectiveCell_T0_delay1_MNO2 == 1
                %selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p(tempSelectiveCellBoolIndex_length1_currentFOV);
                
                [M,I] = sort(glm_r2_currentFOV(tempSelectiveCellBoolIndex_length1_currentFOV,plot_lengthFlag),'descend');
                selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p(tempSelectiveCellBoolIndex_length1_currentFOV);
                selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_lengthx(I);
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
            if if_plotSelectiveCell_T0_delay1_MNO2 == 0
                % abandon now
            elseif if_plotSelectiveCell_T0_delay1_MNO2 == 1
                %selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p(tempSelectiveCellBoolIndex_length2_currentFOV);
                
                if if_0normal_1length123_2length1234 == 0
                    [M,I] = sort(glm_r2_currentFOV(tempSelectiveCellBoolIndex_length2_currentFOV,plot_lengthFlag),'descend');
                elseif if_0normal_1length123_2length1234 == 1
                    [M,I] = sort(glm_r2_currentFOV(tempSelectiveCellBoolIndex_length2_currentFOV,1),'descend');
                end
                selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p(tempSelectiveCellBoolIndex_length2_currentFOV);
                selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_lengthx(I);
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
            if if_plotSelectiveCell_T0_delay1_MNO2 == 0
                % abandon now
            elseif if_plotSelectiveCell_T0_delay1_MNO2 == 1
                %selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p(tempSelectiveCellBoolIndex_length3_currentFOV);
                
                if if_0normal_1length123_2length1234 == 0
                    [M,I] = sort(glm_r2_currentFOV(tempSelectiveCellBoolIndex_length3_currentFOV,plot_lengthFlag),'descend');
                elseif if_0normal_1length123_2length1234 == 1
                    [M,I] = sort(glm_r2_currentFOV(tempSelectiveCellBoolIndex_length3_currentFOV,1),'descend');
                end
                selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p(tempSelectiveCellBoolIndex_length3_currentFOV);
                selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_lengthx(I);
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
            if if_plotSelectiveCell_T0_delay1_MNO2 == 0
                % abandon now
            elseif if_plotSelectiveCell_T0_delay1_MNO2 == 1
                %selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p(tempSelectiveCellBoolIndex_length4_currentFOV);
                
                if if_0normal_1length123_2length1234 == 0
                    [M,I] = sort(glm_r2_currentFOV(tempSelectiveCellBoolIndex_length4_currentFOV,plot_lengthFlag),'descend');
                elseif if_0normal_1length123_2length1234 == 1
                    [M,I] = sort(glm_r2_currentFOV(tempSelectiveCellBoolIndex_length4_currentFOV,1),'descend');
                end
                selectiveCellIndex_suite2p_lengthx = cellIndex_suite2p(tempSelectiveCellBoolIndex_length4_currentFOV);
                selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_lengthx(I);
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
        
        
        % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length1_delay1Bin;
        % selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_length1_T;
        
        if if_plotSelectiveCell_T0_delay1_MNO2 == 2
            % abandon now
        end
        
        if if_plotSelectiveCell == 1
            temp_plotRoiNum = min(plotRoiNum,sum(selectiveCellIndex_suite2p_lengthx>=temp_id_suite2p));
            selectiveCellIndex_suite2p_lengthx = selectiveCellIndex_suite2p_lengthx(selectiveCellIndex_suite2p_lengthx>=temp_id_suite2p);
        else
            temp_plotRoiNum = plotRoiNum;
        end
        
        %if false
        %if true
        if if_plot_svm_beta == 1
            % abandon now
        end
        
        
        if if_plot_mappingTwoSessions == 1
            tempMappingCellIndex_suite2p_withCorr_sorted;
            temp_sortedStartIndex;
            plotRoiNum;
            
            temp_plotRoiNum = min(size(tempMappingCellIndex_suite2p_withCorr_sorted,1)-temp_sortedStartIndex+1,plotRoiNum);
            
            temptempCellIndex_suite2p = tempMappingCellIndex_suite2p_withCorr_sorted(:,tempSessionIndex);
            
            selectiveCellIndex_suite2p_lengthx = temptempCellIndex_suite2p(temp_sortedStartIndex:(temp_sortedStartIndex+temp_plotRoiNum-1));
        end
        
        
        fprintf('plotCellNum = %d.\n',temp_plotRoiNum);
        
        
        %% Plot single cell
        if if_plotSingleCell == 1
            
            tempRoiCount = 0;
            for tempIndex=1:temp_plotRoiNum
                
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
                    % abandon now
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
                    
                    %fig1 = figure('Name',sprintf('fig%d, cell id %d, FOV id %s',tempIndex,temp_cellIndex_suite2p,currentSession(18:22)),'NumberTitle','off');
                    %fig1 = figure('Name',sprintf('length%d, fig%d, cell id %d',plot_lengthFlag,tempIndex,temp_cellIndex_suite2p),'NumberTitle','off');
                    
                    if if_plot_mappingTwoSessions == 0
                        tempIndex2 = tempIndex;
                    elseif if_plot_mappingTwoSessions == 1
                        tempIndex2 = temp_sortedStartIndex+tempIndex-1;
                    end
                    
                    fig1 = figure('Name',sprintf('length%d, fig%d, cell id %d',plot_lengthFlag,tempIndex2,temp_cellIndex_suite2p),'NumberTitle','off');                    
                    set(gcf,'Position',[35+0 42+0 1630 840]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                    t = tiledlayout(4,11,'TileSpacing','tight','Padding','tight'); %#ok<*NASGU>
                    
                    %t.Title.String = sprintf('suite2p cell id=%d                      %s',...
                    %    temp_cellIndex_suite2p,currentSession(18:end));
                    t.Title.String = sprintf('suite2p cell id=%d (valid id=%d)                      %s',...
                       temp_cellIndex_suite2p,temp_id2,currentSession(18:end));                    
                    t.Title.FontSize = 16;
                    t.Title.Interpreter = 'none';
                    
                    
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
                    if if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 ~= 2
                        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
                        y = cellID_F_dff_baselinePeriod_choiceMemory_mean;
                        plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.75,0.25]);
                        hold on
                        y_sem = cellID_F_dff_baselinePeriod_choiceMemory_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.75,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
                        y = cellID_F_dff_baselinePeriod_choiceOffload_mean;
                        plot(x,y,'-','LineWidth',1.5,'Color',[0.75,0.25,0.25]);
                        hold on
                        y_sem = cellID_F_dff_baselinePeriod_choiceOffload_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.75,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
                        y = cellID_F_dff_baselinePeriod_noChoice_mean;
                        plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.25,0.25]);
                        hold on
                        y_sem = cellID_F_dff_baselinePeriod_noChoice_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        
                        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                        y = cellID_F_dff_decisionPeriodA_choiceMemory_mean;
                        plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.75,0.25]);
                        hold on
                        y_sem = cellID_F_dff_decisionPeriodA_choiceMemory_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.75,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                        y = cellID_F_dff_decisionPeriodA_choiceOffload_mean;
                        plot(x,y,'-','LineWidth',1.5,'Color',[0.75,0.25,0.25]);
                        hold on
                        y_sem = cellID_F_dff_decisionPeriodA_choiceOffload_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.75,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                        y = cellID_F_dff_decisionPeriodA_noChoice_mean;
                        plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.25,0.25]);
                        hold on
                        y_sem = cellID_F_dff_decisionPeriodA_noChoice_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        h_line = [];
                        
                        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
                        y = cellID_F_dff_decisionPeriod_choiceMemory_mean;
                        h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.75,0.25])];
                        hold on
                        y_sem = cellID_F_dff_decisionPeriod_choiceMemory_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.75,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
                        y = cellID_F_dff_decisionPeriod_choiceOffload_mean;
                        h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',[0.75,0.25,0.25])];
                        hold on
                        y_sem = cellID_F_dff_decisionPeriod_choiceOffload_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.75,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
                        y = cellID_F_dff_decisionPeriod_noChoice_mean;
                        h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])];
                        hold on
                        y_sem = cellID_F_dff_decisionPeriod_noChoice_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                    elseif if_plot_MNOtrial0_MNOseq1_afterMtrial2_aftrMCrrMO3_aftrMErrMO4 == 2
                        
                        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
                        y = cellID_F_dff_baselinePeriod_afterMemoryCorrect_mean;
                        plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.75,0.25]);
                        hold on
                        y_sem = cellID_F_dff_baselinePeriod_afterMemoryCorrect_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        x = temp_ticks_baselinePeriod(1):temp_ticks_baselinePeriod(end);
                        y = cellID_F_dff_baselinePeriod_afterMemoryError_mean;
                        plot(x,y,'-','LineWidth',1.5,'Color',[0.75,0.25,0.25]);
                        hold on
                        y_sem = cellID_F_dff_baselinePeriod_afterMemoryError_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.75,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                        y = cellID_F_dff_decisionPeriodA_afterMemoryCorrect_mean;
                        plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.75,0.25]);
                        hold on
                        y_sem = cellID_F_dff_decisionPeriodA_afterMemoryCorrect_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        x = temp_ticks_decisionPeriodA(1):temp_ticks_decisionPeriodA(end);
                        y = cellID_F_dff_decisionPeriodA_afterMemoryError_mean;
                        plot(x,y,'-','LineWidth',1.5,'Color',[0.75,0.25,0.25]);
                        hold on
                        y_sem = cellID_F_dff_decisionPeriodA_afterMemoryError_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.75,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        
                        h_line = [];
                        
                        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
                        y = cellID_F_dff_decisionPeriod_afterMemoryCorrect_mean;
                        h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',[0.25,0.75,0.25])];
                        hold on
                        y_sem = cellID_F_dff_decisionPeriod_afterMemoryCorrect_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
                        hold on
                        
                        x = temp_ticks_decisionPeriod(1):temp_ticks_decisionPeriod(end);
                        y = cellID_F_dff_decisionPeriod_afterMemoryError_mean;
                        h_line = [h_line plot(x,y,'-','LineWidth',1.5,'Color',[0.75,0.25,0.25])];
                        hold on
                        y_sem = cellID_F_dff_decisionPeriod_afterMemoryError_sem;
                        patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.75,0.25,0.25],'FaceAlpha',0.2,'EdgeColor',[0.62 0.62 0.62]);
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
                        le = legend(h_line,'Memory','Offload','noChoice','Location','northeast','fontsize',10);%10
                        le.ItemTokenSize = ones(1,3)*20;
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
                    
                elseif if_plot_location0_seq1 == 1
                    % abandon now
                end
                
                % Save fig
                if if_0normal_1length123_2length1234 == 0
                    currentFigName = sprintf('FOVID%s_length%d_fig%d_cellID%d',currentSession(18:22),plot_lengthFlag,tempIndex,temp_cellIndex_suite2p);
                    temp_output_path_singleCell_lengthx = [output_path_singleCell '\' sprintf('length%d',plot_lengthFlag)];
                elseif if_0normal_1length123_2length1234 == 1
                    currentFigName = sprintf('FOVID%s_cellID%d_length%d_fig%d',currentSession(18:22),temp_cellIndex_suite2p,plot_lengthFlag,tempIndex);
                    temp_output_path_singleCell_lengthx = [output_path_singleCell '\' sprintf('length123')];
                end
                
                if if_plot_mappingTwoSessions == 1
                    %currentFigName = sprintf('length%d_fig%d_cellID%d_FOVID%s',plot_lengthFlag,tempIndex,temp_cellIndex_suite2p,currentSession(18:22));                    
                    
                    tempIndex2 = temp_sortedStartIndex+tempIndex-1;

                    currentFigName = sprintf('length%d_fig%d_FOVID%s_cellID%d',plot_lengthFlag,tempIndex2,currentSession(18:22),temp_cellIndex_suite2p);
                end
                
                temp_file = dir(temp_output_path_singleCell_lengthx);
                if isempty(temp_file) == 1
                    mkdir(temp_output_path_singleCell_lengthx);
                end
                fileName_fig = [temp_output_path_singleCell_lengthx '\' currentFigName '.png'];%emf,pdf,png,jpg,tiff
                % to judge whether save figure or not
                if ifSave_fig == 1
                    exportgraphics(fig1,fileName_fig,'Resolution',300);
                    close all
                end
                
            end
        end
        
    end
    fprintf('Time is %.1f secs (roiNum = %d).\n',toc(t0),size(F_dff,1));
    
    if if_singleFOV == 1
        break
    end
end

fprintf('Time of the stepB2_analysis is %.1f secs.\n',toc(time0));
fprintf('num_FOV=%d.\n',num_FOV);


if if_profile == 1
    profile viewer
end

%% End
