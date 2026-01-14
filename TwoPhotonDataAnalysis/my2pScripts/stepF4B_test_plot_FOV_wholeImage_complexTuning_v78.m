% Chuan's 9th script (20251214)
% This script: To conduct spatial organzation analysis in one FOV, related to figure 6.
close all


% if_tuning_location0_precision1 = 0;
if_tuning_location0_precision1 = if_tuning_location0_precision1;

if if_dff_singleSession1_twoSession2_allSession3 == 1
    temp_currentSessionIndex = currentSessionIndex;
elseif if_dff_singleSession1_twoSession2_allSession3 == 2
    temp_currentSessionIndex = currentSessionIndex_B;
end

% if_trialTuning0_seqTuning1_mix2 = 2;%0
if_trialTuning0_seqTuning1_mix2 = if_plotBeta_trialTuning0_seqTuning1_mix2;%2

if_computeShuffle = 1;%1

temp_shuffleNum = 10000;%10000
temp_shuffleNum_A = 100;%100
temp_shuffleNum_B = 1000;%100

threshold_r2_prctile = 0;%0,50,75

% threshold_prctileA = 5;%5-->0.1
% threshold_prctileB = 95;%95-->99.9
threshold_prctileA = 2.5;%5-->0.1-->2.5
threshold_prctileB = 100-threshold_prctileA;%95-->99.9


threshold_prctile_clusterA = 95;

threshold_prctile_clusterB_high = 95;%95,99
threshold_prctile_clusterB_low = 100-threshold_prctile_clusterB_high;



temp_resampleNum = 100;%10000

if_load = 1;%1

if_plot_twoFeatureMerge = 0;
if_plot_contour0_r2Circle1 = 1;

if_plot_twoFeatureSplit = 0;
if_plot_twoFeatureSplitB = 1;

if_plot_centriodComparison = 0;

if_plot_clusterA = 0;
if_plot_clusterB = 1;%1

if_sessionB = 1;%1


if_r2_binary = 1;%1
if_r2Binary_unitView0_subspaceView1 = 0;%0,1

if_plot_distance_pairwise0_centroid1 = 0;%0

if if_sessionB == 1
    cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow = cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryLow; %#ok<*UNRCH>
    cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh = cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryHigh;
    cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh = cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryHigh;
    cellIndex_suite2p_memoryPrecisionLow_choiceMemoryLow = cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryLow;
    cellIndex_suite2p_choiceMemoryBaselineHigh = cellIndex_suite2p_B_choiceMemoryBaselineHigh;
    
    %tempMappingCellIndex_suite2p = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,2);
    tempMappingCellIndex_suite2p = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end);
else
    tempMappingCellIndex_suite2p = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,1);
end
cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow;
cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh;
cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh;
cellIndex_suite2p_memoryPrecisionLow_choiceMemoryLow;
cellIndex_suite2p_choiceMemoryBaselineHigh;

r2_memoryPrecision = r2_memoryPrecision; %#ok<*ASGSL>
r2_choiceMemory = r2_choiceMemory;

tempMappingCellIndex_suite2p;


temptemp_cellIndex_suite2p = tempMappingCellIndex_suite2p;
% temptemp_cellIndex_suite2p = [cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow;...
%     cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh;cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh;cellIndex_suite2p_memoryPrecisionLow_choiceMemoryLow];

% temptemp_cellIndex_suite2p = [cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow;...
%     cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh;cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh];

% temptemp_cellIndex_suite2p = [cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow;cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh;...
%     cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh;cellIndex_suite2p_choiceMemoryBaselineHigh];


r2_memoryPrecision_all = nan(size(r2_memoryPrecision));
r2_choiceMemory_all = nan(size(r2_choiceMemory));
r2_choiceMemory_baseline_all = nan(size(r2_choiceMemory));
r2_memoryPrecision_valid = nan(size(r2_memoryPrecision));
r2_choiceMemory_valid = nan(size(r2_choiceMemory));
r2_choiceMemory_baseline_valid = nan(size(r2_choiceMemory));
for tempi=1:length(temptemp_cellIndex_suite2p)
    tempIndex = find(tempMappingCellIndex_suite2p==temptemp_cellIndex_suite2p(tempi));
    if isempty(tempIndex)
        continue
    end
    
    if if_tuning_location0_precision1 == 1
        if if_trialTuning0_seqTuning1_mix2 == 0
            r2_memoryPrecision_valid(tempIndex) = r2_memoryPrecision(tempIndex);
            r2_choiceMemory_valid(tempIndex) = r2_choiceMemory(tempIndex);
            r2_choiceMemory_baseline_valid(tempIndex) = r2_choiceMemory_baseline(tempIndex);
        elseif if_trialTuning0_seqTuning1_mix2 == 1
            r2_memoryPrecision_valid(tempIndex) = r2_seqPrecision(tempIndex);
            r2_choiceMemory_valid(tempIndex) = r2_rProb(tempIndex);
            r2_choiceMemory_baseline_valid(tempIndex) = r2_rProb_baseline(tempIndex);
        elseif if_trialTuning0_seqTuning1_mix2 == 2
            r2_memoryPrecision_valid(tempIndex) = r2_seqPrecision(tempIndex);
            r2_choiceMemory_valid(tempIndex) = r2_rProb(tempIndex);
            r2_choiceMemory_baseline_valid(tempIndex) = r2_choiceMemory_baseline(tempIndex);
        end
        
    elseif if_tuning_location0_precision1== 0
        if if_trialTuning0_seqTuning1_mix2 == 0
            r2_memoryPrecision_valid(tempIndex) = r2_6loc(tempIndex);
            r2_choiceMemory_valid(tempIndex) = r2_choiceMemory(tempIndex);
            r2_choiceMemory_baseline_valid(tempIndex) = r2_choiceMemory_baseline(tempIndex);
        elseif if_trialTuning0_seqTuning1_mix2 == 1
            r2_memoryPrecision_valid(tempIndex) = r2_6loc(tempIndex);
            r2_choiceMemory_valid(tempIndex) = r2_rProb(tempIndex);
            r2_choiceMemory_baseline_valid(tempIndex) = r2_rProb_baseline(tempIndex);
        elseif if_trialTuning0_seqTuning1_mix2 == 2
            r2_memoryPrecision_valid(tempIndex) = r2_6loc(tempIndex);
            r2_choiceMemory_valid(tempIndex) = r2_rProb(tempIndex);
            r2_choiceMemory_baseline_valid(tempIndex) = r2_choiceMemory_baseline(tempIndex);
        end
    end
    
end

r2_memoryPrecision_valid;
r2_choiceMemory_valid;
r2_choiceMemory_baseline_valid;

threshold_r2_prctile;

temp_threshold = prctile(r2_memoryPrecision_valid,threshold_r2_prctile);
r2_memoryPrecision_valid(r2_memoryPrecision_valid<temp_threshold) = 0;

temp_threshold = prctile(r2_choiceMemory_valid,threshold_r2_prctile);
r2_choiceMemory_valid(r2_choiceMemory_valid<temp_threshold) = 0;

temp_threshold = prctile(r2_choiceMemory_baseline_valid,threshold_r2_prctile);
r2_choiceMemory_baseline_valid(r2_choiceMemory_baseline_valid<temp_threshold) = 0;

if if_tuning_location0_precision1 == 1
    if if_trialTuning0_seqTuning1_mix2 == 0
        r2_memoryPrecision_all = r2_memoryPrecision;
        r2_choiceMemory_all = r2_choiceMemory;
        r2_choiceMemory_baseline_all = r2_choiceMemory_baseline;
    elseif if_trialTuning0_seqTuning1_mix2 == 1
        r2_memoryPrecision_all = r2_seqPrecision;
        r2_choiceMemory_all = r2_rProb;
        r2_choiceMemory_baseline_all = r2_rProb_baseline;
    elseif if_trialTuning0_seqTuning1_mix2 == 2
        r2_memoryPrecision_all = r2_seqPrecision;
        r2_choiceMemory_all = r2_rProb;
        r2_choiceMemory_baseline_all = r2_choiceMemory_baseline;
    end
elseif if_tuning_location0_precision1 == 0
    if if_trialTuning0_seqTuning1_mix2 == 0
        r2_memoryPrecision_all = r2_6loc;
        r2_choiceMemory_all = r2_choiceMemory;
        r2_choiceMemory_baseline_all = r2_choiceMemory_baseline;
    elseif if_trialTuning0_seqTuning1_mix2 == 1
        r2_memoryPrecision_all = r2_6loc;
        r2_choiceMemory_all = r2_rProb;
        r2_choiceMemory_baseline_all = r2_rProb_baseline;
    elseif if_trialTuning0_seqTuning1_mix2 == 2
        r2_memoryPrecision_all = r2_6loc;
        r2_choiceMemory_all = r2_rProb;
        r2_choiceMemory_baseline_all = r2_choiceMemory_baseline;
    end
end

if if_r2_binary == 1
    
    temp_p1 = p_6loc;
    temp_p2 = p_rProb;
    temp_p3 = p_choiceMemory_baseline;
    
    temp1 = zeros(size(r2_memoryPrecision_all));
    temp2 = zeros(size(r2_memoryPrecision_all));
    temp3 = zeros(size(r2_memoryPrecision_all));
    
    temp1(temp_p1<0.01) = 1;
    temp2(temp_p2<0.01) = 1;
    temp3(temp_p3<0.01) = 1;
    
   if if_r2Binary_unitView0_subspaceView1 == 1
       if exist('tempSelectiveBoolIndex_memorySubspace','var') == 1
           temp1 = double(tempSelectiveBoolIndex_memorySubspace);
           temp2 = double(tempSelectiveBoolIndex_metaSubspace);
           temp3 = double(tempSelectiveBoolIndex_priorSubspace);
                      
           if true % new in 20250417
               r2_memoryPrecision_all = abs(beta_memory_linearAxis);
               r2_choiceMemory_all = abs(beta_meta_linearAxis);
               r2_choiceMemory_baseline_all = abs(beta_prior_linearAxis);
               
               r2_memoryPrecision_all(~tempSelectiveBoolIndex_memorySubspace) = 0;
               r2_choiceMemory_all(~tempSelectiveBoolIndex_metaSubspace) = 0;
               r2_choiceMemory_baseline_all(~tempSelectiveBoolIndex_priorSubspace) = 0;
           end
           
       else
           if_r2Binary_unitView0_subspaceView1 = 0;
       end
   end
   
   if if_r2Binary_unitView0_subspaceView1 == 0
       r2_memoryPrecision_all(~temp1) = 0;
       r2_choiceMemory_all(~temp2) = 0;
       r2_choiceMemory_baseline_all(~temp3) = 0;
   end
   fprintf('if_r2Binary_unitView0_subspaceView1 = %d.\n',if_r2Binary_unitView0_subspaceView1);
   
    r2_memoryPrecision_raw = r2_memoryPrecision_all;
    r2_choiceMemory_raw = r2_choiceMemory_all;
    r2_choiceMemory_baseline_raw = r2_choiceMemory_baseline_all;
    
    r2n11n_memoryPrecision_raw = rescale(r2_memoryPrecision_raw,0,1);%0.2
    r2n11n_choiceMemory_raw = rescale(r2_choiceMemory_raw,0,1);
    r2n11n_choiceMemory_baseline_raw = rescale(r2_choiceMemory_baseline_raw,0,1);
    
    %r2n11n_memoryPrecision_raw2 = r2n11n_memoryPrecision_raw .* temp1;
    %r2n11n_choiceMemory_raw2  = r2n11n_choiceMemory_raw .* temp2;
    %r2n11n_choiceMemory_baseline_raw2  = r2n11n_choiceMemory_baseline_raw .* temp3;
    
    r2n11n_memoryPrecision_raw2 = r2n11n_memoryPrecision_raw;
    r2n11n_choiceMemory_raw2  = r2n11n_choiceMemory_raw;
    r2n11n_choiceMemory_baseline_raw2  = r2n11n_choiceMemory_baseline_raw;
    
    r2_memoryPrecision_valid = temp1;
    r2_choiceMemory_valid = temp2;
    r2_choiceMemory_baseline_valid = temp3;
    
    r2_memoryPrecision_all = temp1;
    r2_choiceMemory_all = temp2;
    r2_choiceMemory_baseline_all = temp3;
    
end

r2n11n_memoryPrecision_valid = rescale(r2_memoryPrecision_all,0,1);%0.2
r2n11n_choiceMemory_valid = rescale(r2_choiceMemory_all,0,1);
r2n11n_choiceMemory_baseline_valid = rescale(r2_choiceMemory_baseline_all,0,1);

max_r2_memoryPrecision_valid = max(r2_memoryPrecision_all);
max_r2_choiceMemory_valid = max(r2_choiceMemory_all);
max_r2_choiceMemory_baseline_valid = max(r2_choiceMemory_baseline_all);

cellType = nan(size(r2_memoryPrecision));
cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow;% cellType=1
cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh;% cellType=2
cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh;% cellType=3

for tempi=1:length(cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow)
    tempIndex = find(tempMappingCellIndex_suite2p==cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow(tempi));
    cellType(tempIndex) = 1;
end
for tempi=1:length(cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh)
    tempIndex = find(tempMappingCellIndex_suite2p==cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh(tempi));
    cellType(tempIndex) = 2;
end
for tempi=1:length(cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh)
    tempIndex = find(tempMappingCellIndex_suite2p==cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh(tempi));
    cellType(tempIndex) = 3;
end


tempMappingCellIndex_suite2p;
r2n11n_memoryPrecision_valid;
r2n11n_choiceMemory_valid;
cellType;



%%
temp_range = FOVAllCellRange_multiFOV(temp_currentSessionIndex,1):FOVAllCellRange_multiFOV(temp_currentSessionIndex,2);
cellIndex_suite2p = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range);

% cellIndex_suite2p = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end);
% mark in 20240401, need to modify here to make it proper in 3 sessions FOV



currentSession_multi = string;

if if_monkey_D0_Z1 == 0
    
    currentSession_multi = [currentSession_multi; '113Recording_20230426A_Ding_Site16'];%1
    currentSession_multi = [currentSession_multi; '113Recording_20230427A_Ding_Site16_sameFOV0426'];%2
    currentSession_multi = [currentSession_multi; '113Recording_20230502A_Ding_Site13'];%3
    currentSession_multi = [currentSession_multi; '113Recording_20230503A_Ding_Site13_sameFOV0502'];%4
    currentSession_multi = [currentSession_multi; '113Recording_20230504A_Ding_Site02'];%5
    currentSession_multi = [currentSession_multi; '113Recording_20230508A_Ding_Site02_sameFOV0509'];%6
    currentSession_multi = [currentSession_multi; '113Recording_20230509A_Ding_Site02'];% 660000 frames, easy to crash, 7
    %
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
    %
    currentSession_multi = [currentSession_multi; '113Recording_20230527A_Ding_Site02B'];%19
    currentSession_multi = [currentSession_multi; '113Recording_20230528A_Ding_Site02B_sameFOV0527'];%20
    currentSession_multi = [currentSession_multi; '113Recording_20230530A_Ding_Site05C'];%21
    currentSession_multi = [currentSession_multi; '113Recording_20230531A_Ding_Site05C_sameFOV0530'];%22
    currentSession_multi = [currentSession_multi; '113Recording_20230601A_Ding_Site13B'];%23
    currentSession_multi = [currentSession_multi; '113Recording_20230602A_Ding_Site13B_sameFOV0601'];%24
    %
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



targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)


temp_currentSession = currentSession_multi{temp_currentSessionIndex};
fprintf('temp_currentSession = %s.\n',temp_currentSession);


output_shortPath = 'D:\twoPhotonData_motionCorrected';
temp_currentSession_path = [output_shortPath '\' temp_currentSession];
temp_if_max0_min1 = 0;
output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
path_plane = [output_path,'\plane0'];

temp_if_max0_min1 = 0;
%template_path = autoGetFileName_general('template*.tif', output_path,temp_if_max0_min1);
template_path = autoGetFileName_general('maxProjection*.tif', output_path,temp_if_max0_min1);
template = double(loadtiff(template_path));

if if_load== 1
    fileName_Fall = 'Fall.mat';
    fileName_iscell = 'iscell.npy';
    fullFileName_Fall = [path_plane,'\',fileName_Fall];
    fullFileName_iscell = [path_plane,'\',fileName_iscell];
    
    iscell = readNPY(fullFileName_iscell);
    
    s = load(fullFileName_Fall,'stat');
    roi_stats_raw = s.stat;
    temp_cellIndex = find(iscell(:,1)==1);
    roi_stats = roi_stats_raw(temp_cellIndex);
end


temp_roiNum = roiNum;

roi_med = nan(length(r2n11n_memoryPrecision_valid),2); % (y,x)
for tempj=1:temp_roiNum
    if isnan(r2n11n_memoryPrecision_valid(tempj))
        continue
    end
    
    temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
    
    temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
    temp_roi_xpix = double(temp_roi_stat.xpix);
    temp_roi_ypix = double(temp_roi_stat.ypix);
    temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
    
    roi_med(tempj,1) = temp_roi_med(1);
    roi_med(tempj,2) = temp_roi_med(2);
end


inmin = 0;%0
inmax = 1500;%843-->1500

outmin = 55;%30-->130
outmax = 250;%200

% temp_template2 = rescale(template,30,200,'InputMin',inmin,'InputMax',inmax);
temp_template2 = rescale(template,outmin,outmax,'InputMin',inmin,'InputMax',inmax);



%% Shulffle

% temptempBoolIndex = ~isnan(r2_memoryPrecision_all);
% tempSum1 = sum(temptempBoolIndex);
%
% temptempIndex = find(temptempBoolIndex==true);
%
% r2_memoryPrecision_valid_shuffled = nan(size(r2_memoryPrecision_all,1),temp_shuffleNum);
% r2_choiceMemory_valid_shuffled = nan(size(r2_choiceMemory_all,1),temp_shuffleNum);
% for tempi=1:temp_shuffleNum
%     temptempIndex_shuffled = temptempIndex(randperm(length(temptempIndex))');
%
%     temp_r2_memoryPrecision_valid_shuffled = r2_memoryPrecision_all;
%     temp_r2_memoryPrecision_valid_shuffled(temptempIndex) = r2_memoryPrecision_all(temptempIndex_shuffled);
%
%     temp_r2_choiceMemory_valid_shuffled = r2_choiceMemory_all;
%     temp_r2_choiceMemory_valid_shuffled(temptempIndex) = r2_choiceMemory_all(temptempIndex_shuffled);
%
%     r2_memoryPrecision_valid_shuffled(:,tempi) = temp_r2_memoryPrecision_valid_shuffled;
%     r2_choiceMemory_valid_shuffled(:,tempi) = temp_r2_choiceMemory_valid_shuffled;
% end
% % r2_memoryPrecision_valid_shuffled_mean = mean(r2_memoryPrecision_valid_shuffled,2);
% % r2_choiceMemory_valid_shuffled_mean = mean(r2_choiceMemory_valid_shuffled,2);



temptempIndex = (1:size(roi_med,1))';

roi_med_y_shuffled = nan(size(roi_med,1),temp_shuffleNum);
roi_med_x_shuffled = nan(size(roi_med,1),temp_shuffleNum);
for tempi=1:temp_shuffleNum
    temptempIndex_shuffled = temptempIndex(randperm(length(temptempIndex))');
    
    temp_roi_med_y_shuffled = roi_med(:,1);
    temp_roi_med_y_shuffled(temptempIndex) = roi_med(temptempIndex_shuffled,1);
    
    temp_roi_med_x_shuffled = roi_med(:,2);
    temp_roi_med_x_shuffled(temptempIndex) = roi_med(temptempIndex_shuffled,2);
    
    roi_med_y_shuffled(:,tempi) = temp_roi_med_y_shuffled;
    roi_med_x_shuffled(:,tempi) = temp_roi_med_x_shuffled;
end
roi_med_shuffled = nan(size(roi_med,1),2,temp_shuffleNum);
roi_med_shuffled(:,1,:) = roi_med_y_shuffled;
roi_med_shuffled(:,2,:) = roi_med_x_shuffled;


%% Centriod, offloading-tuning roi
% fprintf('----------------\n');
% fprintf('All offloading-tuning roi: \n');

% temp_r2A = r2_memoryPrecision_all;
% temp_r2B = r2_choiceMemory_all;
temp_r2A = r2_memoryPrecision_valid;
temp_r2B = r2_choiceMemory_valid;
temp_r2C = r2_choiceMemory_baseline_valid;

temp1 = roi_med(:,2) .* temp_r2A ./ sum(temp_r2A,'omitnan');
centriod_memoryPrecision_x = sum(temp1,'omitnan');
temp1 = roi_med(:,1) .* temp_r2A ./ sum(temp_r2A,'omitnan');
centriod_memoryPrecision_y = sum(temp1,'omitnan');
% fprintf('centriod_memoryPrecision: %.1f, %.1f \n',centriod_memoryPrecision_x,centriod_memoryPrecision_y);

temp1 = roi_med(:,2) .* temp_r2B ./ sum(temp_r2B,'omitnan');
centriod_choiceMemory_x = sum(temp1,'omitnan');
temp1 = roi_med(:,1) .* temp_r2B ./ sum(temp_r2B,'omitnan');
centriod_choiceMemory_y = sum(temp1,'omitnan');
% fprintf('centriod_choiceMemory:    %.1f, %.1f \n',centriod_choiceMemory_x,centriod_choiceMemory_y);

temp1 = roi_med(:,2) .* temp_r2C ./ sum(temp_r2C,'omitnan');
centriod_choiceMemory_baseline_x = sum(temp1,'omitnan');
temp1 = roi_med(:,1) .* temp_r2C ./ sum(temp_r2C,'omitnan');
centriod_choiceMemory_baseline_y = sum(temp1,'omitnan');

centriod_memoryPrecision_all = [centriod_memoryPrecision_x,centriod_memoryPrecision_y];
centriod_meta_all = [centriod_choiceMemory_x,centriod_choiceMemory_y];
centriod_meta_baseline_all = [centriod_choiceMemory_baseline_x,centriod_choiceMemory_baseline_y];


% temp_r2A_shuffled = r2_memoryPrecision_valid_shuffled;
% temp_r2B_shuffled = r2_choiceMemory_valid_shuffled;
% temp_r2A_shuffled = temp_r2A;
% temp_r2B_shuffled = temp_r2B;
temp_r2A_shuffled = r2_memoryPrecision_all;
temp_r2B_shuffled = r2_choiceMemory_all;
temp_r2C_shuffled = r2_choiceMemory_baseline_all;

temp_roi_med_x_shuffled = roi_med_x_shuffled;
temp_roi_med_y_shuffled = roi_med_y_shuffled;

temp1 = temp_roi_med_x_shuffled .* temp_r2A_shuffled ./ sum(temp_r2A_shuffled,1,'omitnan');
centriod_memoryPrecision_x = sum(temp1,1,'omitnan');
temp1 = temp_roi_med_y_shuffled .* temp_r2A_shuffled ./ sum(temp_r2A_shuffled,1,'omitnan');
centriod_memoryPrecision_y = sum(temp1,1,'omitnan');
centriod_memoryPrecision_all_shuffled = [centriod_memoryPrecision_x;centriod_memoryPrecision_y]';

temp1 = temp_roi_med_x_shuffled .* temp_r2B_shuffled ./ sum(temp_r2B_shuffled,1,'omitnan');
centriod_choiceMemory_x = sum(temp1,1,'omitnan');
temp1 = temp_roi_med_y_shuffled .* temp_r2B_shuffled ./ sum(temp_r2B_shuffled,1,'omitnan');
centriod_choiceMemory_y = sum(temp1,1,'omitnan');
centriod_meta_all_shuffled = [centriod_choiceMemory_x;centriod_choiceMemory_y]';

temp1 = temp_roi_med_x_shuffled .* temp_r2C_shuffled ./ sum(temp_r2C_shuffled,1,'omitnan');
centriod_choiceMemory_baseline_x = sum(temp1,1,'omitnan');
temp1 = temp_roi_med_y_shuffled .* temp_r2C_shuffled ./ sum(temp_r2C_shuffled,1,'omitnan');
centriod_choiceMemory_baseline_y = sum(temp1,1,'omitnan');
centriod_meta_baseline_all_shuffled = [centriod_choiceMemory_baseline_x;centriod_choiceMemory_baseline_y]';


%% Centriod, offloading-tuning & memoryPrecisionHigh_choiceMemoryHigh
% fprintf('----------------\n');
% fprintf('Offloading-tuning & memoryPrecisionHigh_choiceMemoryHigh roi: \n');
temptempBoolIndex = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);

temp_r2A = r2_memoryPrecision_all;
temp_r2A(~temptempBoolIndex) = nan;
temp_r2B = r2_choiceMemory_all;
temp_r2B(~temptempBoolIndex) = nan;


temp1 = roi_med(:,2) .* temp_r2A ./ sum(temp_r2A,'omitnan');
centriod_memoryPrecision_x = sum(temp1,'omitnan');
temp1 = roi_med(:,1) .* temp_r2A ./ sum(temp_r2A,'omitnan');
centriod_memoryPrecision_y = sum(temp1,'omitnan');
% fprintf('centriod_memoryPrecision: %.1f, %.1f \n',centriod_memoryPrecision_x,centriod_memoryPrecision_y);

temp1 = roi_med(:,2) .* temp_r2B ./ sum(temp_r2B,'omitnan');
centriod_choiceMemory_x = sum(temp1,'omitnan');
temp1 = roi_med(:,1) .* temp_r2B ./ sum(temp_r2B,'omitnan');
centriod_choiceMemory_y = sum(temp1,'omitnan');
% fprintf('centriod_choiceMemory:    %.1f, %.1f \n',centriod_choiceMemory_x,centriod_choiceMemory_y);

centriod_memoryPrecision_mix = [centriod_memoryPrecision_x,centriod_memoryPrecision_y];
centriod_meta_mix = [centriod_choiceMemory_x,centriod_choiceMemory_y];


% temp_r2A_shuffled = r2_memoryPrecision_valid_shuffled;
% temp_r2A_shuffled(~temptempBoolIndex,:) = nan;
% temp_r2B_shuffled = r2_choiceMemory_valid_shuffled;
% temp_r2B_shuffled(~temptempBoolIndex,:) = nan;
temp_r2A_shuffled = temp_r2A;
temp_r2B_shuffled = temp_r2B;
temp_roi_med_x_shuffled = roi_med_x_shuffled;
temp_roi_med_y_shuffled = roi_med_y_shuffled;

temp1 = temp_roi_med_x_shuffled .* temp_r2A_shuffled ./ sum(temp_r2A_shuffled,1,'omitnan');
centriod_memoryPrecision_x = sum(temp1,1,'omitnan');
temp1 = temp_roi_med_y_shuffled .* temp_r2A_shuffled ./ sum(temp_r2A_shuffled,1,'omitnan');
centriod_memoryPrecision_y = sum(temp1,1,'omitnan');
centriod_memoryPrecision_mix_shuffled = [centriod_memoryPrecision_x;centriod_memoryPrecision_y]';

temp1 = temp_roi_med_x_shuffled .* temp_r2B_shuffled ./ sum(temp_r2B_shuffled,1,'omitnan');
centriod_choiceMemory_x = sum(temp1,1,'omitnan');
temp1 = temp_roi_med_y_shuffled .* temp_r2B_shuffled ./ sum(temp_r2B_shuffled,1,'omitnan');
centriod_choiceMemory_y = sum(temp1,1,'omitnan');
centriod_meta_mix_shuffled = [centriod_choiceMemory_x;centriod_choiceMemory_y]';


%% Centriod, offloading-tuning & memoryPrecisionHigh_choiceMemoryLow
% fprintf('----------------\n');
% fprintf('Offloading-tuning & memoryPrecisionHigh_choiceMemoryLow roi: \n');
temptempBoolIndex = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow);

temp_r2A = r2_memoryPrecision_all;
temp_r2A(~temptempBoolIndex) = nan;

temp1 = roi_med(:,2) .* temp_r2A ./ sum(temp_r2A,'omitnan');
centriod_memoryPrecision_x = sum(temp1,'omitnan');
temp1 = roi_med(:,1) .* temp_r2A ./ sum(temp_r2A,'omitnan');
centriod_memoryPrecision_y = sum(temp1,'omitnan');
% fprintf('centriod_memoryPrecision: %.1f, %.1f \n',centriod_memoryPrecision_x,centriod_memoryPrecision_y);

centriod_memoryPrecision_pure = [centriod_memoryPrecision_x,centriod_memoryPrecision_y];


% temp_r2A_shuffled = r2_memoryPrecision_valid_shuffled;
% temp_r2A_shuffled(~temptempBoolIndex,:) = nan;
temp_r2A_shuffled = temp_r2A;
temp_roi_med_x_shuffled = roi_med_x_shuffled;
temp_roi_med_y_shuffled = roi_med_y_shuffled;

temp1 = temp_roi_med_x_shuffled .* temp_r2A_shuffled ./ sum(temp_r2A_shuffled,1,'omitnan');
centriod_memoryPrecision_x = sum(temp1,1,'omitnan');
temp1 = temp_roi_med_y_shuffled .* temp_r2A_shuffled ./ sum(temp_r2A_shuffled,1,'omitnan');
centriod_memoryPrecision_y = sum(temp1,1,'omitnan');
centriod_memoryPrecision_pure_shuffled = [centriod_memoryPrecision_x;centriod_memoryPrecision_y]';



%% Centriod, offloading-tuning & memoryPrecisionLow_choiceMemoryHigh
% fprintf('----------------\n');
% fprintf('Offloading-tuning & memoryPrecisionLow_choiceMemoryHigh roi: \n');
temptempBoolIndex = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh);

temp_r2B = r2_choiceMemory_all;
temp_r2B(~temptempBoolIndex) = nan;

temp1 = roi_med(:,2) .* temp_r2B ./ sum(temp_r2B,'omitnan');
centriod_choiceMemory_x = sum(temp1,'omitnan');
temp1 = roi_med(:,1) .* temp_r2B ./ sum(temp_r2B,'omitnan');
centriod_choiceMemory_y = sum(temp1,'omitnan');
% fprintf('centriod_choiceMemory:    %.1f, %.1f \n',centriod_choiceMemory_x,centriod_choiceMemory_y);

centriod_meta_pure = [centriod_choiceMemory_x,centriod_choiceMemory_y];


% temp_r2B_shuffled = r2_choiceMemory_valid_shuffled;
% temp_r2B_shuffled(~temptempBoolIndex,:) = nan;
temp_r2B_shuffled = temp_r2B;
temp_roi_med_x_shuffled = roi_med_x_shuffled;
temp_roi_med_y_shuffled = roi_med_y_shuffled;

temp1 = temp_roi_med_x_shuffled .* temp_r2B_shuffled ./ sum(temp_r2B_shuffled,1,'omitnan');
centriod_choiceMemory_x = sum(temp1,1,'omitnan');
temp1 = temp_roi_med_y_shuffled .* temp_r2B_shuffled ./ sum(temp_r2B_shuffled,1,'omitnan');
centriod_choiceMemory_y = sum(temp1,1,'omitnan');
centriod_meta_pure_shuffled = [centriod_choiceMemory_x;centriod_choiceMemory_y]';


%% Centriod distance
centriod_memoryPrecision_all;
centriod_meta_all;
centriod_memoryPrecision_mix;
centriod_meta_mix;
centriod_memoryPrecision_pure;
centriod_meta_pure;

centriod_memoryPrecision_all_shuffled;
centriod_meta_all_shuffled;
centriod_memoryPrecision_mix_shuffled;
centriod_meta_mix_shuffled;
centriod_memoryPrecision_pure_shuffled;
centriod_meta_pure_shuffled;

x1 = centriod_memoryPrecision_all_shuffled;
x2 = centriod_meta_all_shuffled;
y = sum(abs(x1-x2).^2,2).^0.5;
centriodDis_all_shuffled = y;

x1 = centriod_memoryPrecision_mix_shuffled;
x2 = centriod_meta_mix_shuffled;
y = sum(abs(x1-x2).^2,2).^0.5;
centriodDis_mix_shuffled = y;

x1 = centriod_memoryPrecision_pure_shuffled;
x2 = centriod_meta_pure_shuffled;
y = sum(abs(x1-x2).^2,2).^0.5;
centriodDis_pure_shuffled = y;

x1 = centriod_memoryPrecision_all_shuffled;
x2 = centriod_meta_baseline_all_shuffled;
y = sum(abs(x1-x2).^2,2).^0.5;
centriodDis_all_AC_shuffled = y;

x1 = centriod_meta_all_shuffled;
x2 = centriod_meta_baseline_all_shuffled;
y = sum(abs(x1-x2).^2,2).^0.5;
centriodDis_all_BC_shuffled = y;

centriodDis_all_shuffled;
centriodDis_mix_shuffled;
centriodDis_pure_shuffled;


centriodDis_all_shuffled_prctileA = prctile(centriodDis_all_shuffled,threshold_prctileA);
centriodDis_mix_shuffled_prctileA = prctile(centriodDis_mix_shuffled,threshold_prctileA);
centriodDis_pure_shuffled_prctileA = prctile(centriodDis_pure_shuffled,threshold_prctileA);
centriodDis_all_AC_shuffled_prctileA = prctile(centriodDis_all_AC_shuffled,threshold_prctileA);
centriodDis_all_BC_shuffled_prctileA = prctile(centriodDis_all_BC_shuffled,threshold_prctileA);

centriodDis_all_shuffled_prctileB = prctile(centriodDis_all_shuffled,threshold_prctileB);
centriodDis_mix_shuffled_prctileB = prctile(centriodDis_mix_shuffled,threshold_prctileB);
centriodDis_pure_shuffled_prctileB = prctile(centriodDis_pure_shuffled,threshold_prctileB);
centriodDis_all_AC_shuffled_prctileB = prctile(centriodDis_all_AC_shuffled,threshold_prctileB);
centriodDis_all_BC_shuffled_prctileB = prctile(centriodDis_all_BC_shuffled,threshold_prctileB);


x1 = centriod_memoryPrecision_all;
x2 = centriod_meta_all;
y = sum(abs(x1-x2).^2,2).^0.5;
centriodDis_all = y;

x1 = centriod_memoryPrecision_mix;
x2 = centriod_meta_mix;
y = sum(abs(x1-x2).^2,2).^0.5;
centriodDis_mix = y;

x1 = centriod_memoryPrecision_pure;
x2 = centriod_meta_pure;
y = sum(abs(x1-x2).^2,2).^0.5;
centriodDis_pure = y;

x1 = centriod_memoryPrecision_all;
x2 = centriod_meta_baseline_all;
y = sum(abs(x1-x2).^2,2).^0.5;
centriodDis_AC_all = y;

x1 = centriod_meta_all;
x2 = centriod_meta_baseline_all;
y = sum(abs(x1-x2).^2,2).^0.5;
centriodDis_BC_all = y;

fprintf('centriodDis_all  = %.2f (shuffle=[%.2f,%.2f]).\n',centriodDis_all,centriodDis_all_shuffled_prctileA,centriodDis_all_shuffled_prctileB);
% fprintf('centriodDis_mix  = %.2f (shuffle=[%.2f,%.2f]).\n',centriodDis_mix,centriodDis_all_shuffled_prctileA,centriodDis_all_shuffled_prctileB);
% fprintf('centriodDis_pure = %.2f (shuffle=[%.2f,%.2f]).\n',centriodDis_pure,centriodDis_all_shuffled_prctileA,centriodDis_all_shuffled_prctileB);
fprintf('centriodDis_AC_all  = %.2f (shuffle=[%.2f,%.2f]).\n',centriodDis_AC_all,centriodDis_all_AC_shuffled_prctileA,centriodDis_all_AC_shuffled_prctileB);
fprintf('centriodDis_BC_all  = %.2f (shuffle=[%.2f,%.2f]).\n',centriodDis_BC_all,centriodDis_all_BC_shuffled_prctileA,centriodDis_all_BC_shuffled_prctileB);


%% Centriod comparison
temp_resampleNum;


r2_memoryPrecision_all;

temptempBoolIndex = ~isnan(r2_memoryPrecision_all);
temptempIndex = find(temptempBoolIndex==true);

roi_med_resampled = nan(size(roi_med,1),size(roi_med,2),temp_resampleNum);
r2_memoryPrecision_valid_resampled = nan(size(r2_memoryPrecision_all,1),temp_resampleNum);
r2_choiceMemory_valid_resampled = nan(size(r2_memoryPrecision_all,1),temp_resampleNum);

temptempIndex_resampled_multi = nan(length(temptempIndex),temp_resampleNum);

for tempi=1:temp_resampleNum
    tempResampleIndex = sort(randi(length(temptempIndex),1,length(temptempIndex)))';
    temptempIndex_resampled = temptempIndex(tempResampleIndex);
    
    temptempIndex_resampled_multi(:,tempi) = temptempIndex_resampled;
    
    temp_roi_med_resampled = roi_med;
    temp_roi_med_resampled(temptempIndex,:) = roi_med(temptempIndex_resampled,:);
    roi_med_resampled(:,:,tempi) = temp_roi_med_resampled;
    
    temp_r2_memoryPrecision_valid_resampled = r2_memoryPrecision_all;
    temp_r2_memoryPrecision_valid_resampled(temptempIndex) = r2_memoryPrecision_all(temptempIndex_resampled);
    r2_memoryPrecision_valid_resampled(:,tempi) = temp_r2_memoryPrecision_valid_resampled;
    
    temp_r2_choiceMemory_valid_resampled = r2_choiceMemory_all;
    temp_r2_choiceMemory_valid_resampled(temptempIndex) = r2_choiceMemory_all(temptempIndex_resampled);
    r2_choiceMemory_valid_resampled(:,tempi) = temp_r2_choiceMemory_valid_resampled;
end


%% Resampled centriod, offloading-tuning roi

temp_r2A_resampled = r2_memoryPrecision_valid_resampled;
temp_r2B_resampled = r2_choiceMemory_valid_resampled;

temp1 = squeeze(roi_med_resampled(:,2,:)) .* temp_r2A_resampled ./ sum(temp_r2A_resampled,1,'omitnan');
centriod_memoryPrecision_x = sum(temp1,1,'omitnan');
temp1 = squeeze(roi_med_resampled(:,1,:)) .* temp_r2A_resampled ./ sum(temp_r2A_resampled,1,'omitnan');
centriod_memoryPrecision_y = sum(temp1,1,'omitnan');
centriod_memoryPrecision_all_resampled = [centriod_memoryPrecision_x;centriod_memoryPrecision_y]';

temp1 = squeeze(roi_med_resampled(:,2,:)) .* temp_r2B_resampled ./ sum(temp_r2B_resampled,1,'omitnan');
centriod_choiceMemory_x = sum(temp1,1,'omitnan');
temp1 = squeeze(roi_med_resampled(:,1,:)) .* temp_r2B_resampled ./ sum(temp_r2B_resampled,1,'omitnan');
centriod_choiceMemory_y = sum(temp1,1,'omitnan');
centriod_meta_all_resampled = [centriod_choiceMemory_x;centriod_choiceMemory_y]';



%% Resampled centriod, offloading-tuning & memoryPrecisionHigh_choiceMemoryHigh
temptempBoolIndex_current = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
temptempIndex_current = find(temptempBoolIndex_current==true);


temp_r2A_resampled = r2_memoryPrecision_valid_resampled;
% temp_r2A_resampled(~temptempBoolIndex,:) = nan;
temp_r2B_resampled = r2_choiceMemory_valid_resampled;
% temp_r2B_resampled(~temptempBoolIndex,:) = nan;

for tempi=1:size(temptempIndex_resampled_multi,2)
    temp1 = ismember(temptempIndex_resampled_multi(:,tempi),temptempIndex_current);
    temp2 = temptempIndex(temp1);
    
    temptempBoolIndex_current_B = false(length(temptempBoolIndex_current),1);
    temptempBoolIndex_current_B(temp2) = true;
    
    temp_r2A_resampled(~temptempBoolIndex_current_B,tempi) = nan;
    temp_r2B_resampled(~temptempBoolIndex_current_B,tempi) = nan;
end

temp1 = squeeze(roi_med_resampled(:,2,:)) .* temp_r2A_resampled ./ sum(temp_r2A_resampled,1,'omitnan');
centriod_memoryPrecision_x = sum(temp1,1,'omitnan');
temp1 = squeeze(roi_med_resampled(:,1,:)) .* temp_r2A_resampled ./ sum(temp_r2A_resampled,1,'omitnan');
centriod_memoryPrecision_y = sum(temp1,1,'omitnan');
centriod_memoryPrecision_mix_resampled = [centriod_memoryPrecision_x;centriod_memoryPrecision_y]';

temp1 = squeeze(roi_med_resampled(:,2,:)) .* temp_r2B_resampled ./ sum(temp_r2B_resampled,1,'omitnan');
centriod_choiceMemory_x = sum(temp1,1,'omitnan');
temp1 = squeeze(roi_med_resampled(:,1,:)) .* temp_r2B_resampled ./ sum(temp_r2B_resampled,1,'omitnan');
centriod_choiceMemory_y = sum(temp1,1,'omitnan');
centriod_meta_mix_resampled = [centriod_choiceMemory_x;centriod_choiceMemory_y]';


%% Resampled centriod, offloading-tuning & memoryPrecisionHigh_choiceMemoryLow
temptempBoolIndex_current = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow);
temptempIndex_current = find(temptempBoolIndex_current==true);

temp_r2A_resampled = r2_memoryPrecision_valid_resampled;
% temp_r2A_resampled(~temptempBoolIndex_current,:) = nan;

for tempi=1:size(temptempIndex_resampled_multi,2)
    temp1 = ismember(temptempIndex_resampled_multi(:,tempi),temptempIndex_current);
    temp2 = temptempIndex(temp1);
    
    temptempBoolIndex_current_B = false(length(temptempBoolIndex_current),1);
    temptempBoolIndex_current_B(temp2) = true;
    
    temp_r2A_resampled(~temptempBoolIndex_current_B,tempi) = nan;
end

temp1 = squeeze(roi_med_resampled(:,2,:)) .* temp_r2A_resampled ./ sum(temp_r2A_resampled,1,'omitnan');
centriod_memoryPrecision_x = sum(temp1,1,'omitnan');
temp1 = squeeze(roi_med_resampled(:,1,:)) .* temp_r2A_resampled ./ sum(temp_r2A_resampled,1,'omitnan');
centriod_memoryPrecision_y = sum(temp1,1,'omitnan');
centriod_memoryPrecision_pure_resampled = [centriod_memoryPrecision_x;centriod_memoryPrecision_y]';



%% Resampled centriod, offloading-tuning & memoryPrecisionLow_choiceMemoryHigh
temptempBoolIndex_current = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh);
temptempIndex_current = find(temptempBoolIndex_current==true);

temp_r2B_resampled = r2_choiceMemory_valid_resampled;
% temp_r2B_resampled(~temptempBoolIndex_current,:) = nan;

for tempi=1:size(temptempIndex_resampled_multi,2)
    temp1 = ismember(temptempIndex_resampled_multi(:,tempi),temptempIndex_current);
    temp2 = temptempIndex(temp1);
    
    temptempBoolIndex_current_B = false(length(temptempBoolIndex_current),1);
    temptempBoolIndex_current_B(temp2) = true;
    
    temp_r2B_resampled(~temptempBoolIndex_current_B,tempi) = nan;
end

temp1 = squeeze(roi_med_resampled(:,2,:)) .* temp_r2B_resampled ./ sum(temp_r2B_resampled,1,'omitnan');
centriod_choiceMemory_x = sum(temp1,1,'omitnan');
temp1 = squeeze(roi_med_resampled(:,1,:)) .* temp_r2B_resampled ./ sum(temp_r2B_resampled,1,'omitnan');
centriod_choiceMemory_y = sum(temp1,1,'omitnan');
centriod_meta_pure_resampled = [centriod_choiceMemory_x;centriod_choiceMemory_y]';


%%
centriod_memoryPrecision_all_resampled;
centriod_meta_all_resampled;
centriod_memoryPrecision_mix_resampled;
centriod_meta_mix_resampled;
centriod_memoryPrecision_pure_resampled;
centriod_meta_pure_resampled;

centriod_memoryPrecision_all_resampled_xMean = mean(centriod_memoryPrecision_all_resampled(:,1));
centriod_memoryPrecision_all_resampled_yMean = mean(centriod_memoryPrecision_all_resampled(:,2));

centriod_meta_all_resampled_xMean = mean(centriod_meta_all_resampled(:,1));
centriod_meta_all_resampled_yMean = mean(centriod_meta_all_resampled(:,2));

centriod_memoryPrecision_mix_resampled_xMean = mean(centriod_memoryPrecision_mix_resampled(:,1));
centriod_memoryPrecision_mix_resampled_yMean = mean(centriod_memoryPrecision_mix_resampled(:,2));

centriod_meta_mix_resampled_xMean = mean(centriod_meta_mix_resampled(:,1));
centriod_meta_mix_resampled_yMean = mean(centriod_meta_mix_resampled(:,2));

centriod_memoryPrecision_pure_resampled_xMean = mean(centriod_memoryPrecision_pure_resampled(:,1));
centriod_memoryPrecision_pure_resampled_yMean = mean(centriod_memoryPrecision_pure_resampled(:,2));

centriod_meta_pure_resampled_xMean = mean(centriod_meta_pure_resampled(:,1));
centriod_meta_pure_resampled_yMean = mean(centriod_meta_pure_resampled(:,2));

[~,temp_p_x_all] = ttest2(centriod_memoryPrecision_all_resampled(:,1),centriod_meta_all_resampled(:,1));
[~,temp_p_y_all] = ttest2(centriod_memoryPrecision_all_resampled(:,2),centriod_meta_all_resampled(:,2));

[~,temp_p_x_mix] = ttest2(centriod_memoryPrecision_mix_resampled(:,1),centriod_meta_mix_resampled(:,1));
[~,temp_p_y_mix] = ttest2(centriod_memoryPrecision_mix_resampled(:,2),centriod_meta_mix_resampled(:,2));

[~,temp_p_x_pure] = ttest2(centriod_memoryPrecision_pure_resampled(:,1),centriod_meta_pure_resampled(:,1));
[~,temp_p_y_pure] = ttest2(centriod_memoryPrecision_pure_resampled(:,2),centriod_meta_pure_resampled(:,2));



%% roi anatomy info
if if_plot_twoFeatureMerge == 1
    fig = figure('Name',' ','NumberTitle','off');
    % set(gcf,'Position',[400 50 800 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[700 50 300 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    
    nexttile
    
    %% plot whole FOV
    % temp_template = template / max(template,[],'all');
    
    %temp_template2 = (temp_template.^0.7)*255;
    % temp_template2 = (temp_template.^0.55)*255;
    
    
    % inmin = 0;
    % inmax = 1500;%843
    
    image(temp_template2);
    hold on
    colormap(gray);
    
    
    %% plot roi
    if if_plot_contour0_r2Circle1 == 0
        temp_bound_offset = 1;%0-->1
        
        for tempi=1:3
            if tempi == 1
                temp_selectiveCellIndex_suite2p = cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow;
                temp_roiColor = [252,141,89]/255;
            elseif tempi == 2
                temp_selectiveCellIndex_suite2p = cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh;
                temp_roiColor = [145,191,219]/255;
            elseif tempi == 3
                temp_selectiveCellIndex_suite2p = cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh;
                temp_roiColor = [255,255,191]/255;
            end
            
            temp_selectiveCellIndex_suite2p;
            
            temp_roiNum = length(temp_selectiveCellIndex_suite2p);
            
            for tempj=1:temp_roiNum
                temp_cellIndex_suite2p = temp_selectiveCellIndex_suite2p(tempj);
                
                temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
                temp_roi_xpix = double(temp_roi_stat.xpix);
                temp_roi_ypix = double(temp_roi_stat.ypix);
                
                temp_I = false(512,512);
                for tempk=1:length(temp_roi_xpix)
                    temp_I(temp_roi_ypix(tempk),temp_roi_xpix(tempk)) = true;
                end
                
                %se = strel('disk',3);%1,2,3
                se = strel('disk',temp_bound_offset);%1,2,3
                temp_I = imdilate(temp_I,se);
                
                B = bwboundaries(temp_I,'noholes');
                for k = 1:length(B)
                    boundary = B{k};
                    %plot(boundary(:,2)+1,boundary(:,1)+1,'color',temp_roiColor, 'LineWidth', 0.5);%1.5-->0.5
                    patch(boundary(:,2)+1,boundary(:,1)+1,temp_roiColor);
                    hold on
                end
                
            end
            
        end
        
        axis equal;
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'Visible','off');
        
    elseif if_plot_contour0_r2Circle1 == 1
        
        tempMappingCellIndex_suite2p;
        r2n11n_memoryPrecision_valid;
        r2n11n_choiceMemory_valid;
        cellType;
        
        
        r2n11n_memoryPrecision_valid_size = r2n11n_memoryPrecision_valid.^0.00001;
        %r2n11n_memoryPrecision_valid_size = rescale(r2n11n_memoryPrecision_valid_size,7,40);
        r2n11n_memoryPrecision_valid_size = rescale(r2n11n_memoryPrecision_valid_size,0.2,1);
        
        
        r2n11n_choiceMemory_valid_size = r2n11n_choiceMemory_valid.^0.00001;%0.00001
        %r2n11n_choiceMemory_valid_size = rescale(r2n11n_choiceMemory_valid_size,0.2,1);
        r2n11n_choiceMemory_valid_size = rescale(r2n11n_choiceMemory_valid_size,7,40);
        
        colormap parula %parula, cool, jet
        cmap = colormap;
        
        colormap gray
        
        %r2n11n_choiceMemory_valid_sizeB = round(r2n11n_choiceMemory_valid_size*255)+1;
        %[tempB_min,tempB_max] = bounds(r2n11n_choiceMemory_valid_sizeB);
        r2n11n_memoryPrecision_valid_sizeB = round(r2n11n_memoryPrecision_valid_size*255)+1;
        [tempB_min,tempB_max] = bounds(r2n11n_memoryPrecision_valid_sizeB);
        
        
        
        
        temp_roiNum = length(tempMappingCellIndex_suite2p);
        
        for tempj=1:temp_roiNum
            if isnan(r2n11n_memoryPrecision_valid(tempj))
                continue
            end
            
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            
            %temp_size = r2n11n_memoryPrecision_valid_size(tempj);
            %temp_roiColor = cmap(r2n11n_choiceMemory_valid_sizeB(tempj),:);
            
            temp_size = r2n11n_choiceMemory_valid_size(tempj);
            temp_roiColor = cmap(r2n11n_memoryPrecision_valid_sizeB(tempj),:);
            
            
            %if cellType(tempj) == 3
            %    scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temp_roiColor,'filled','MarkerEdgeColor',[1 1 1]);
            %else
            %    scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %end
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            hold on
            
        end
        
        
        axis equal;
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'Visible','off');
        
    end
    
end

%% Plot
temp_roiColor_purePrecision = [252,141,89]/255;
temp_roiColor_pureMeta = [145,191,219]/255;
% temp_roiColor_mix = [255,255,191]/255;
temp_roiColor_mix = 0.85*[255,255,191]/255;

if if_plot_twoFeatureSplit == 1
    fig = figure('Name',' ','NumberTitle','off');
    set(gcf,'Position',[700 500 480 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[700 500 720 300*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
    
    temp_roiNum = length(tempMappingCellIndex_suite2p);
    
    % temp_roiColor = [0.9290 0.6940 0.1250];
    
    temp_size_min = 10;%7
    temp_size_max = 40;%40
    
    temp_lineWidth = 1;
    
    %% Plot A, memoryPrecision
    nexttile
    
    image(temp_template2);
    hold on
    colormap(gray);
    
    % r2n11n_memoryPrecision_valid_size = r2n11n_memoryPrecision_valid.^0.00001;
    r2n11n_memoryPrecision_valid_size = r2n11n_memoryPrecision_valid.^1;
    r2n11n_memoryPrecision_valid_size = rescale(r2n11n_memoryPrecision_valid_size,temp_size_min,temp_size_max);
    
    % roi_med = nan(length(r2n11n_memoryPrecision_valid),2); % (y,x)
    
    temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
    temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow);
    temptempBoolIndex12 = temptempBoolIndex1 | temptempBoolIndex2;
    
    for tempj=1:temp_roiNum
        %if isnan(r2n11n_memoryPrecision_valid(tempj))
        if temptempBoolIndex12(tempj) == false
            continue
        end
        
        
        temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
        
        temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
        temp_roi_xpix = double(temp_roi_stat.xpix);
        temp_roi_ypix = double(temp_roi_stat.ypix);
        temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
        
        %roi_med(tempj,1) = temp_roi_med(1);
        %roi_med(tempj,2) = temp_roi_med(2);
        
        temp_size = r2n11n_memoryPrecision_valid_size(tempj);
        
        %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
        
        if cellType(tempj) == 3
            temptemp_roiColor = temp_roiColor_mix;
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temp_roiColor_mix,'filled','MarkerEdgeColor',[1 1 1]*0);
        else
            temptemp_roiColor = temp_roiColor_purePrecision;
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temp_roiColor_purePrecision,'filled','MarkerEdgeColor',[1 1 1]*0);
        end
        scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
        %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
        hold on
        
    end
    
    y_min = 0;
    y_max = 511;
    
    x_min = 0;
    x_max = 511;
    
    %
    %     temp1 = centriod_memoryPrecision_all;
    %     plot([1 1].*temp1(1),[y_min y_max],...
    %         'LineWidth',temp_lineWidth,'color',[0.75 0.25 0.25]);%[0 0.4470 0.7410]
    %     hold on
    %     plot([x_min x_max],[1 1].*temp1(2),...
    %         'LineWidth',temp_lineWidth,'color',[0.75 0.25 0.25]);%[0 0.4470 0.7410]
    %     hold on
    
    temp1 = centriod_memoryPrecision_mix;
    plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*0.08 temp1(2)+(y_max-y_min)*0.08],...
        'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
    hold on
    plot([temp1(1)-(y_max-y_min)*0.08 temp1(1)+(y_max-y_min)*0.08],[1 1].*temp1(2),...
        'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
    hold on
    
    
    temp1 = centriod_memoryPrecision_pure;
    plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*0.08 temp1(2)+(y_max-y_min)*0.08],...
        'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
    hold on
    plot([temp1(1)-(y_max-y_min)*0.08 temp1(1)+(y_max-y_min)*0.08],[1 1].*temp1(2),...
        'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
    hold on
    
    
    axis equal;
    xlim([x_min x_max])
    ylim([y_min y_max])
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    % set(gca, 'YDir','reverse');
    set(gca,'Visible','off');
    temp_title = title(sprintf('Memory precision tuning'),'FontSize',11);
    temp_title.Visible = 'on';
    
    %% Plot B
    nexttile
    
    image(temp_template2);
    hold on
    colormap(gray);
    
    % r2n11n_choiceMemory_valid_size = r2n11n_choiceMemory_valid.^0.00001;%0.00001
    r2n11n_choiceMemory_valid_size = r2n11n_choiceMemory_valid.^1;%0.00001
    % r2n11n_choiceMemory_valid_size = rescale(r2n11n_choiceMemory_valid_size,1,40);
    r2n11n_choiceMemory_valid_size = rescale(r2n11n_choiceMemory_valid_size,temp_size_min,temp_size_max*(max_r2_choiceMemory_valid/max_r2_memoryPrecision_valid));
    
    temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
    temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh);
    temptempBoolIndex12 = temptempBoolIndex1 | temptempBoolIndex2;
    
    for tempj=1:temp_roiNum
        %if isnan(r2n11n_memoryPrecision_valid(tempj))
        if temptempBoolIndex12(tempj) == false
            continue
        end
        
        
        
        temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
        
        temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
        temp_roi_xpix = double(temp_roi_stat.xpix);
        temp_roi_ypix = double(temp_roi_stat.ypix);
        temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
        
        temp_size = r2n11n_choiceMemory_valid_size(tempj);
        
        %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
        
        if cellType(tempj) == 3
            temptemp_roiColor = temp_roiColor_mix;
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temp_roiColor_mix,'filled','MarkerEdgeColor',[1 1 1]*0);
        else
            temptemp_roiColor = temp_roiColor_pureMeta;
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temp_roiColor_pureMeta,'filled','MarkerEdgeColor',[1 1 1]*0);
        end
        scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
        %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
        hold on
        
    end
    
    y_min = 0;
    y_max = 511;
    
    x_min = 0;
    x_max = 511;
    
    
    %     temp1 = centriod_meta_all;
    %     plot([1 1].*temp1(1),[y_min y_max],...
    %         'LineWidth',temp_lineWidth,'color',[0.75 0.25 0.25]);%[0 0.4470 0.7410]
    %     hold on
    %     plot([x_min x_max],[1 1].*temp1(2),...
    %         'LineWidth',temp_lineWidth,'color',[0.75 0.25 0.25]);%[0 0.4470 0.7410]
    %     hold on
    
    temp1 = centriod_meta_mix;
    plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*0.08 temp1(2)+(y_max-y_min)*0.08],...
        'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
    hold on
    plot([temp1(1)-(y_max-y_min)*0.08 temp1(1)+(y_max-y_min)*0.08],[1 1].*temp1(2),...
        'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
    hold on
    
    
    temp1 = centriod_meta_pure;
    plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*0.08 temp1(2)+(y_max-y_min)*0.08],...
        'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
    hold on
    plot([temp1(1)-(y_max-y_min)*0.08 temp1(1)+(y_max-y_min)*0.08],[1 1].*temp1(2),...
        'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
    hold on
    
    
    axis equal;
    xlim([x_min x_max])
    ylim([y_min y_max])
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    % set(gca, 'YDir','reverse');
    set(gca,'Visible','off');
    temp_title = title(sprintf('Meta-memory tuning'),'FontSize',11);
    temp_title.Visible = 'on';
    
    
    %% Plot centriod distance distribution
    fig = figure('Name',' ','NumberTitle','off');
    set(gcf,'Position',[1200 520 240 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    h = histogram(centriodDis_all_shuffled,'Normalization','pdf');
    h.FaceColor = [1 1 1]*0.75;
    hold on
    
    [y_min,y_max] = bounds(h.Values);
    
    plot([1 1].*centriodDis_all_shuffled_prctileA,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
        'LineWidth',2,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
    hold on
    plot([1 1].*centriodDis_all_shuffled_prctileB,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
        'LineWidth',2,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
    hold on
    
    %     plot([1 1].*centriodDis_all,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
    %         'LineWidth',2,'color',[0.75 0.25 0.25]);%[0 0.4470 0.7410]
    %     hold on
    plot([1 1].*centriodDis_mix,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
        'LineWidth',5,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
    hold on
    plot([1 1].*centriodDis_pure,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
        'LineWidth',5,'color',[0.25 0.25 0.75]);%[0 0.4470 0.7410]
    hold on
    
    set(gca,'linewidth',1.5)
    %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
    ylim([y_min-(y_max-y_min)*0.00 y_max+(y_max-y_min)*0.00]);
    set(gca, 'FontSize', 12)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Centriod distance', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Pdf', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('Centriod distance'),'FontSize',11);
    temp_title.Interpreter = 'none';
    
end


%% Plot centriod comparison
if if_plot_centriodComparison == 1
    %close all
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[50+0 50+0 240 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    
    centriod_memoryPrecision_mix_resampled;
    centriod_meta_mix_resampled;
    centriod_memoryPrecision_pure_resampled;
    centriod_meta_pure_resampled;
    
    temp_p_x_mix;
    temp_p_y_mix;
    temp_p_x_pure;
    temp_p_y_pure;
    
    %% Mix
    nexttile
    
    temp_centriod1 = centriod_memoryPrecision_mix_resampled;
    temp_centriod2 = centriod_meta_mix_resampled;
    
    temp_p_x = temp_p_x_mix;
    temp_p_y = temp_p_y_mix;
    
    temp_centriod1_mean = [mean(temp_centriod1(:,1)) mean(temp_centriod1(:,2))];
    temp_centriod2_mean = [mean(temp_centriod2(:,1)) mean(temp_centriod2(:,2))];
    
    temp_centriod1_sem = [std(temp_centriod1(:,1))./sqrt(size(temp_centriod1,1)) ...
        std(temp_centriod1(:,2))./sqrt(size(temp_centriod1,1))];
    temp_centriod2_sem = [std(temp_centriod2(:,1))./sqrt(size(temp_centriod2,1)) ...
        std(temp_centriod2(:,2))./sqrt(size(temp_centriod2,1))];
    
    y_bar = [temp_centriod1_mean; temp_centriod2_mean]';
    
    y_min = min(y_bar,[],'all') - max([temp_centriod1_sem temp_centriod1_sem]);
    y_max = max(y_bar,[],'all') + max([temp_centriod1_sem temp_centriod1_sem]);
    
    bar_width = 0.8;
    b = bar([0 1],y_bar,bar_width,'grouped',...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    b(1).FaceColor = temp_roiColor_purePrecision;
    b(2).FaceColor = temp_roiColor_pureMeta;
    
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(y_bar);
    % Get the x coordinate of the bars
    x_bar = nan(ngroups, nbars);
    for i = 1:nbars
        x_bar(:,i) = b(i).XEndPoints;
    end
    
    errorbar(x_bar, y_bar,[temp_centriod1_sem; temp_centriod2_sem]', '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 7);
    hold on
    
    %     le = legend('Precision','Meta',...
    %         'Location','northeast','fontsize',9);
    %     le.ItemTokenSize = ones(1,2)*10;
    
    tempTxt_x = sprintf('');
    if temp_p_x < 0.001
        tempTxt_x = sprintf('***');
    elseif temp_p_x < 0.01
        tempTxt_x = sprintf('**');
    elseif temp_p_x < 0.05
        tempTxt_x = sprintf('*');
    end
    
    tempTxt_y = sprintf('');
    if temp_p_y < 0.001
        tempTxt_y = sprintf('***');
    elseif temp_p_y < 0.01
        tempTxt_y = sprintf('**');
    elseif temp_p_y < 0.05
        tempTxt_y = sprintf('*');
    end
    
    text(0,max(y_bar(1,:))+(y_max-y_min)*0.1,tempTxt_x, 'fontsize',12,'FontWeight','bold','HorizontalAlignment','center');
    text(1,max(y_bar(2,:))+(y_max-y_min)*0.1,tempTxt_y, 'fontsize',12,'FontWeight','bold','HorizontalAlignment','center');
    
    
    set(gca,'linewidth',1.5)
    ylim([y_min-(y_max-y_min)*0.15 y_max+(y_max-y_min)*0.15]);
    set(gca, 'FontSize', 12) %14
    set(gca,'box','off');% 取消右、上边框
    
    
    %temp_category1 = 'X';
    %temp_category2 = 'Y';
    set(gca,'XTick',0:1);
    %set(gca,'XTickLabel', [temp_category1; temp_category2]);%给坐标加标签
    set(gca,'XTickLabel','');
    set(gca,'YTick',round([y_min y_max]));
    
    ylabel('Centriod', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('Mix tuning'),'fontsize',12);
    temp_title.Interpreter = 'none';
    
    
    %% Pure
    nexttile
    
    temp_centriod1 = centriod_memoryPrecision_pure_resampled;
    temp_centriod2 = centriod_meta_pure_resampled;
    
    temp_p_x = temp_p_x_pure;
    temp_p_y = temp_p_y_pure;
    
    temp_centriod1_mean = [mean(temp_centriod1(:,1)) mean(temp_centriod1(:,2))];
    temp_centriod2_mean = [mean(temp_centriod2(:,1)) mean(temp_centriod2(:,2))];
    
    temp_centriod1_sem = [std(temp_centriod1(:,1))./sqrt(size(temp_centriod1,1)) ...
        std(temp_centriod1(:,2))./sqrt(size(temp_centriod1,1))];
    temp_centriod2_sem = [std(temp_centriod2(:,1))./sqrt(size(temp_centriod2,1)) ...
        std(temp_centriod2(:,2))./sqrt(size(temp_centriod2,1))];
    
    y_bar = [temp_centriod1_mean; temp_centriod2_mean]';
    
    y_min = min(y_bar,[],'all') - max([temp_centriod1_sem temp_centriod1_sem]);
    y_max = max(y_bar,[],'all') + max([temp_centriod1_sem temp_centriod1_sem]);
    
    bar_width = 0.8;
    b = bar([0 1],y_bar,bar_width,'grouped',...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    b(1).FaceColor = temp_roiColor_purePrecision;
    b(2).FaceColor = temp_roiColor_pureMeta;
    
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(y_bar);
    % Get the x coordinate of the bars
    x_bar = nan(ngroups, nbars);
    for i = 1:nbars
        x_bar(:,i) = b(i).XEndPoints;
    end
    
    errorbar(x_bar, y_bar,[temp_centriod1_sem; temp_centriod2_sem]', '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 7);
    hold on
    
    le = legend('Precision','Meta',...
        'Location','northeast','fontsize',9);
    le.ItemTokenSize = ones(1,2)*10;
    
    tempTxt_x = sprintf('');
    if temp_p_x < 0.001
        tempTxt_x = sprintf('***');
    elseif temp_p_x < 0.01
        tempTxt_x = sprintf('**');
    elseif temp_p_x < 0.05
        tempTxt_x = sprintf('*');
    end
    
    tempTxt_y = sprintf('');
    if temp_p_y < 0.001
        tempTxt_y = sprintf('***');
    elseif temp_p_y < 0.01
        tempTxt_y = sprintf('**');
    elseif temp_p_y < 0.05
        tempTxt_y = sprintf('*');
    end
    
    text(0,max(y_bar(1,:))+(y_max-y_min)*0.1,tempTxt_x, 'fontsize',12,'FontWeight','bold','HorizontalAlignment','center');
    text(1,max(y_bar(2,:))+(y_max-y_min)*0.1,tempTxt_y, 'fontsize',12,'FontWeight','bold','HorizontalAlignment','center');
    
    
    set(gca,'linewidth',1.5)
    ylim([y_min-(y_max-y_min)*0.15 y_max+(y_max-y_min)*0.15]);
    set(gca, 'FontSize', 12) %14
    set(gca,'box','off');% 取消右、上边框
    
    
    temp_category1 = 'X';
    temp_category2 = 'Y';
    set(gca,'XTick',0:1);
    set(gca,'XTickLabel', [temp_category1; temp_category2]);%给坐标加标签
    set(gca,'YTick',round([y_min y_max]));
    
    ylabel('Centriod', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('Pure tuning'),'fontsize',12);
    temp_title.Interpreter = 'none';
    
end





if if_plot_twoFeatureSplitB == 1
    
    temp_size_min = 5;%2-->5
    temp_size_max = 40;%40-->70-->50-->40
    
    temp_size_threshold_prctile = 87.5;%50,75,85
    
    if if_r2_binary == 1
        temp_size_threshold_prctile = 0;
        temp_size_min = 8;
        temp_size_max = 9;
    end
    
    temp_lineWidth = 1.5;%1.5-->3-->1.5
    temp_lineWidthB = 3;
    
    
    %     y_min = 0;
    %     y_max = 511;
    %
    %     x_min = 0;
    %     x_max = 511;
    
    y_min = 1;
    y_max = 512;
    
    x_min = 1;
    x_max = 512;
    
    r2n11n_memoryPrecision_valid_size = r2n11n_memoryPrecision_valid.^1;
    r2n11n_memoryPrecision_valid_size = rescale(r2n11n_memoryPrecision_valid_size,temp_size_min,temp_size_max);
    
    temp_size_threshold = prctile(r2n11n_memoryPrecision_valid_size,temp_size_threshold_prctile);
    temp1 = r2n11n_memoryPrecision_valid_size>=temp_size_threshold;
    temp2 = r2n11n_memoryPrecision_valid_size(temp1);
    temp2_n11n = rescale(temp2,temp_size_min,temp_size_max);
    r2n11n_memoryPrecision_valid_size(temp1) = temp2_n11n;
    r2n11n_memoryPrecision_valid_size(~temp1) = temp_size_min;
    
    temp1_mean = mean(r2n11n_memoryPrecision_valid_size(r2n11n_memoryPrecision_valid_size>temp_size_min));
    
    
    
    r2n11n_choiceMemory_valid_size = r2n11n_choiceMemory_valid.^1;%0.00001
    %r2n11n_choiceMemory_valid_size = rescale(r2n11n_choiceMemory_valid_size,temp_size_min,temp_size_max*(max_r2_choiceMemory_valid/max_r2_memoryPrecision_valid));
    r2n11n_choiceMemory_valid_size = rescale(r2n11n_choiceMemory_valid_size,temp_size_min,temp_size_max);
    
    temp_size_threshold = prctile(r2n11n_choiceMemory_valid_size,temp_size_threshold_prctile);
    temp1 = r2n11n_choiceMemory_valid_size>=temp_size_threshold;
    temp2 = r2n11n_choiceMemory_valid_size(temp1);
    temp2_n11n = rescale(temp2,temp_size_min,temp_size_max);
    r2n11n_choiceMemory_valid_size(temp1) = temp2_n11n;
    r2n11n_choiceMemory_valid_size(~temp1) = temp_size_min;
    
    temp2_mean = mean(r2n11n_choiceMemory_valid_size(r2n11n_choiceMemory_valid_size>temp_size_min));
    
    r2n11n_choiceMemory_valid_size(r2n11n_choiceMemory_valid_size<=temp_size_min) = 0;
    r2n11n_choiceMemory_valid_size = r2n11n_choiceMemory_valid_size * ...
        (temp1_mean/temp2_mean);
    r2n11n_choiceMemory_valid_size(r2n11n_choiceMemory_valid_size==0) = temp_size_min;
    
    
    
    r2n11n_choiceMemory_baseline_valid_size = r2n11n_choiceMemory_baseline_valid.^1;%0.00001
    %r2n11n_choiceMemory_baseline_valid_size = rescale(r2n11n_choiceMemory_baseline_valid_size,temp_size_min,temp_size_max*(max_r2_choiceMemory_baseline_valid/max_r2_memoryPrecision_valid));
    r2n11n_choiceMemory_baseline_valid_size = rescale(r2n11n_choiceMemory_baseline_valid_size,temp_size_min,temp_size_max);
    
    temp_size_threshold = prctile(r2n11n_choiceMemory_baseline_valid_size,temp_size_threshold_prctile);
    temp1 = r2n11n_choiceMemory_baseline_valid_size>=temp_size_threshold;
    temp2 = r2n11n_choiceMemory_baseline_valid_size(temp1);
    temp2_n11n = rescale(temp2,temp_size_min,temp_size_max);
    r2n11n_choiceMemory_baseline_valid_size(temp1) = temp2_n11n;
    r2n11n_choiceMemory_baseline_valid_size(~temp1) = temp_size_min;
    
    temp3_mean = mean(r2n11n_choiceMemory_baseline_valid_size(r2n11n_choiceMemory_baseline_valid_size>temp_size_min));
    
    r2n11n_choiceMemory_baseline_valid_size(r2n11n_choiceMemory_baseline_valid_size<=temp_size_min) = 0;
    r2n11n_choiceMemory_baseline_valid_size = r2n11n_choiceMemory_baseline_valid_size * ...
        (temp1_mean/temp3_mean);
    r2n11n_choiceMemory_baseline_valid_size(r2n11n_choiceMemory_baseline_valid_size==0) = temp_size_min;
    
    if false
        %% Plot spatial organization in one panel
        fig = figure('Name',' ','NumberTitle','off');
        %set(gcf,'Position',[700 100 480 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[700 100 720 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 100 240 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 100 240*0.7 300*0.7]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[400 100 240*0.7*0.8*0.9*1.01 300*0.7*0.8*0.9*1.01]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        temp_roiNum = length(tempMappingCellIndex_suite2p);
        
        
        
        nexttile
        
        
        temp_template_onlyBorder = temp_template2;
        
        temp_template_onlyBorder(:,:) = 255;
        temp1 = 6+2;%10
        temp_template_onlyBorder(1:temp1,:) = 0;
        temp_template_onlyBorder(end+1-temp1:end,:) = 0;
        temp_template_onlyBorder(:,1:temp1) = 0;
        temp_template_onlyBorder(:,end+1-temp1:end) = 0;
        
        
        %image(temp_template2);
        image(temp_template_onlyBorder);
        hold on
        colormap(gray);
        
        
        
        %%
        r2n11n_memoryPrecision_valid_size_B = r2n11n_memoryPrecision_valid_size;
        r2n11n_choiceMemory_valid_size_B = r2n11n_choiceMemory_valid_size;
        r2n11n_choiceMemory_baseline_valid_size_B = r2n11n_choiceMemory_baseline_valid_size;
        
        for tempj=1:temp_roiNum
            %temp1 = r2n11n_memoryPrecision_valid_size_B(tempj);
            %temp2 = r2n11n_choiceMemory_valid_size_B(tempj);
            %temp3 = r2n11n_choiceMemory_baseline_valid_size_B(tempj);
            
            temp1 = r2n11n_memoryPrecision_raw(tempj);
            temp2 = r2n11n_choiceMemory_raw(tempj);
            temp3 = r2n11n_choiceMemory_baseline_raw(tempj);
            
            temp123 = [temp1,temp2,temp3];
            [M,I] = max(temp123);
            
            if I == 1
                r2n11n_choiceMemory_valid_size_B(tempj) = temp_size_min;
                r2n11n_choiceMemory_baseline_valid_size_B(tempj) = temp_size_min;
            elseif I == 2
                r2n11n_memoryPrecision_valid_size_B(tempj) = temp_size_min;
                r2n11n_choiceMemory_baseline_valid_size_B(tempj) = temp_size_min;
            elseif I == 3
                r2n11n_memoryPrecision_valid_size_B(tempj) = temp_size_min;
                r2n11n_choiceMemory_valid_size_B(tempj) = temp_size_min;
            end
        end
        
        %% Memory precision
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_memoryPrecision_valid_size_B(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_purePrecision;
            
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
            
        end
        
        
        %% Meta-memory
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_valid_size_B(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_pureMeta;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
            
        end
        
        
        
        %% Baseline meta
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_baseline_valid_size_B(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_mix;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
        end
        
        
        
        %%
        
        temptemp_coeff = 0.08;
        temptemp_coeff_B = 0.1;
        
        temp1 = centriod_memoryPrecision_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff_B temp1(2)+(y_max-y_min)*temptemp_coeff_B],...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff_B temp1(1)+(y_max-y_min)*temptemp_coeff_B],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        temp1 = centriod_meta_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff_B temp1(1)+(y_max-y_min)*temptemp_coeff_B],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff_B temp1(2)+(y_max-y_min)*temptemp_coeff_B],...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        temp1 = centriod_meta_baseline_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff_B temp1(1)+(y_max-y_min)*temptemp_coeff_B],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff_B temp1(2)+(y_max-y_min)*temptemp_coeff_B],...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        
        
        
        temp1 = centriod_memoryPrecision_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff temp1(2)+(y_max-y_min)*temptemp_coeff],...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_meta_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff temp1(2)+(y_max-y_min)*temptemp_coeff],...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_meta_baseline_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff temp1(2)+(y_max-y_min)*temptemp_coeff],...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
        hold on
        
        
        
        
        temp1 = centriod_meta_baseline_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff temp1(1)+(y_max-y_min)*temptemp_coeff],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_meta_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff temp1(1)+(y_max-y_min)*temptemp_coeff],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_memoryPrecision_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff temp1(1)+(y_max-y_min)*temptemp_coeff],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
        hold on
        
        
        
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'Visible','off');
        
    end
    
    if false
        fig = figure('Name',' ','NumberTitle','off');
        %set(gcf,'Position',[700 100 480 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[700 100 720 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[700 100 720*0.7 300*0.7]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[700 100 720*0.7*0.93 300*0.7*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[700 100 720*0.7*0.93*0.9*0.95 300*0.7*0.93*0.9*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
        
        temp_roiNum = length(tempMappingCellIndex_suite2p);
        
        %         %temp_size_min = 5;%7,10
        %         %temp_size_max = 80;%40
        %
        %         temp_size_min = 2;%7,10,5
        %         temp_size_max = 40;%40,60
        %
        %         %temp_size_min = 1;%7,10
        %         %temp_size_max = 40;%40
        %
        %         temp_size_threshold_prctile = 75;%50,75
        %
        %         temp_lineWidth = 3;%1.5
        
        %% Plot A, memoryPrecision
        nexttile
        
        image(temp_template2);
        hold on
        colormap(gray);
        
        %         r2n11n_memoryPrecision_valid_size = r2n11n_memoryPrecision_valid.^1;
        %         r2n11n_memoryPrecision_valid_size = rescale(r2n11n_memoryPrecision_valid_size,temp_size_min,temp_size_max);
        %
        %         temp_size_threshold = prctile(r2n11n_memoryPrecision_valid_size,temp_size_threshold_prctile);
        %         temp1 = r2n11n_memoryPrecision_valid_size>temp_size_threshold;
        %         temp2 = r2n11n_memoryPrecision_valid_size(temp1);
        %         temp2_n11n = rescale(temp2,temp_size_min,temp_size_max);
        %         r2n11n_memoryPrecision_valid_size(temp1) = temp2_n11n;
        %         r2n11n_memoryPrecision_valid_size(~temp1) = temp_size_min;
        %
        %         temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
        %         temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow);
        %         temptempBoolIndex12 = temptempBoolIndex1 | temptempBoolIndex2;
        
        for tempj=1:temp_roiNum
            %if temptempBoolIndex12(tempj) == false
            %    continue
            %end
            
            
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_memoryPrecision_valid_size(tempj);
            
            %if temp_size < temp_size_threshold
            %    %continue
            %    temp_size = temp_size_min;
            %end
            if temp_size <= temp_size_min
                continue
            end
            
            %if cellType(tempj) == 3
            %    temptemp_roiColor = temp_roiColor_mix;
            %else
            %    temptemp_roiColor = temp_roiColor_purePrecision;
            %end
            temptemp_roiColor = temp_roiColor_purePrecision;
            
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %         scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*1);
            %         scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            
            %         if temptempBoolIndex12(tempj) == true
            %             scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %         else
            %             scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            %             %if temp_size > 30
            %                 a = 1;
            %             %end
            %         end
            hold on
            
        end
        
        y_min = 0;
        y_max = 511;
        
        x_min = 0;
        x_max = 511;
        
        
        %         temp1 = centriod_memoryPrecision_all;
        %         plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*0.08 temp1(2)+(y_max-y_min)*0.08],...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
        %         hold on
        %         plot([temp1(1)-(y_max-y_min)*0.08 temp1(1)+(y_max-y_min)*0.08],[1 1].*temp1(2),...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
        %         hold on
        
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'Visible','off');
        if if_tuning_location0_precision1 == 1
            temp_title = title(sprintf('Memory precision'),'FontSize',9);
        elseif if_tuning_location0_precision1 == 0
            temp_title = title(sprintf('WM'),'FontSize',9);
        end
        temp_title.Visible = 'on';
        
        %% Plot B, meta-memory
        nexttile
        
        image(temp_template2);
        hold on
        colormap(gray);
        
        %         r2n11n_choiceMemory_valid_size = r2n11n_choiceMemory_valid.^1;%0.00001
        %         %r2n11n_choiceMemory_valid_size = rescale(r2n11n_choiceMemory_valid_size,temp_size_min,temp_size_max*(max_r2_choiceMemory_valid/max_r2_memoryPrecision_valid));
        %         r2n11n_choiceMemory_valid_size = rescale(r2n11n_choiceMemory_valid_size,temp_size_min,temp_size_max);
        %
        %
        %         temp_size_threshold = prctile(r2n11n_choiceMemory_valid_size,temp_size_threshold_prctile);
        %         temp1 = r2n11n_choiceMemory_valid_size>temp_size_threshold;
        %         temp2 = r2n11n_choiceMemory_valid_size(temp1);
        %         temp2_n11n = rescale(temp2,temp_size_min,temp_size_max);
        %         r2n11n_choiceMemory_valid_size(temp1) = temp2_n11n;
        %         r2n11n_choiceMemory_valid_size(~temp1) = temp_size_min;
        %
        %
        %         temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
        %         temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh);
        %         temptempBoolIndex12 = temptempBoolIndex1 | temptempBoolIndex2;
        
        for tempj=1:temp_roiNum
            %if temptempBoolIndex12(tempj) == false
            %    continue
            %end
            
            
            
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_valid_size(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            
            %if cellType(tempj) == 3
            %    temptemp_roiColor = temp_roiColor_mix;
            %else
            %    temptemp_roiColor = temp_roiColor_pureMeta;
            %end
            temptemp_roiColor = temp_roiColor_pureMeta;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %         scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*1);
            %         scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
            
        end
        
        y_min = 0;
        y_max = 511;
        
        x_min = 0;
        x_max = 511;
        
        
        %         temp1 = centriod_meta_all;
        %         plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*0.08 temp1(2)+(y_max-y_min)*0.08],...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
        %         hold on
        %         plot([temp1(1)-(y_max-y_min)*0.08 temp1(1)+(y_max-y_min)*0.08],[1 1].*temp1(2),...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
        %         hold on
        
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        % set(gca, 'YDir','reverse');
        set(gca,'Visible','off');
        temp_title = title(sprintf('Meta-WM'),'FontSize',9);
        temp_title.Visible = 'on';
        
        
        %% Plot C, baseline meta-memory
        nexttile
        
        image(temp_template2);
        hold on
        colormap(gray);
        
        %         r2n11n_choiceMemory_baseline_valid_size = r2n11n_choiceMemory_baseline_valid.^1;%0.00001
        %         %r2n11n_choiceMemory_baseline_valid_size = rescale(r2n11n_choiceMemory_baseline_valid_size,temp_size_min,temp_size_max*(max_r2_choiceMemory_baseline_valid/max_r2_memoryPrecision_valid));
        %         r2n11n_choiceMemory_baseline_valid_size = rescale(r2n11n_choiceMemory_baseline_valid_size,temp_size_min,temp_size_max);
        %
        %
        %         temp_size_threshold = prctile(r2n11n_choiceMemory_baseline_valid_size,temp_size_threshold_prctile);
        %         temp1 = r2n11n_choiceMemory_baseline_valid_size>temp_size_threshold;
        %         temp2 = r2n11n_choiceMemory_baseline_valid_size(temp1);
        %         temp2_n11n = rescale(temp2,temp_size_min,temp_size_max);
        %         r2n11n_choiceMemory_baseline_valid_size(temp1) = temp2_n11n;
        %         r2n11n_choiceMemory_baseline_valid_size(~temp1) = temp_size_min;
        
        
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_baseline_valid_size(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            
            temptemp_roiColor = temp_roiColor_mix;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            hold on
            
        end
        
        y_min = 0;
        y_max = 511;
        
        x_min = 0;
        x_max = 511;
        
        
        %         temp1 = centriod_meta_baseline_all;
        %         plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*0.08 temp1(2)+(y_max-y_min)*0.08],...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
        %         hold on
        %         plot([temp1(1)-(y_max-y_min)*0.08 temp1(1)+(y_max-y_min)*0.08],[1 1].*temp1(2),...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
        %         hold on
        
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        % set(gca, 'YDir','reverse');
        set(gca,'Visible','off');
        %temp_title = title(sprintf('Baseline Meta'),'FontSize',9);
        temp_title = title(sprintf('Baseline'),'FontSize',9);
        temp_title.Visible = 'on';
        
    end
    
    
    if true
        fig = figure('Name',' ','NumberTitle','off');
        %set(gcf,'Position',[700 100 480 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[700 100 720 300]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[700 100 720*0.7 300*0.7]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[700 100 720*0.7*0.93 300*0.7*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[700 100 720*0.7*0.93*0.9*0.95*1.2*1.05 300*0.7*0.93*0.9*0.95*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
        
        temp_roiNum = length(tempMappingCellIndex_suite2p);
        
        %         %temp_size_min = 5;%7,10
        %         %temp_size_max = 80;%40
        %
        %         temp_size_min = 2;%7,10,5
        %         temp_size_max = 40;%40,60
        %
        %         %temp_size_min = 1;%7,10
        %         %temp_size_max = 40;%40
        %
        %         temp_size_threshold_prctile = 75;%50,75
        %
        %         temp_lineWidth = 3;%1.5
        
        %% Plot A, memoryPrecision
        nexttile
        
        image(temp_template2);
        hold on
        colormap(gray);
        
        %         r2n11n_memoryPrecision_valid_size = r2n11n_memoryPrecision_valid.^1;
        %         r2n11n_memoryPrecision_valid_size = rescale(r2n11n_memoryPrecision_valid_size,temp_size_min,temp_size_max);
        %
        %         temp_size_threshold = prctile(r2n11n_memoryPrecision_valid_size,temp_size_threshold_prctile);
        %         temp1 = r2n11n_memoryPrecision_valid_size>temp_size_threshold;
        %         temp2 = r2n11n_memoryPrecision_valid_size(temp1);
        %         temp2_n11n = rescale(temp2,temp_size_min,temp_size_max);
        %         r2n11n_memoryPrecision_valid_size(temp1) = temp2_n11n;
        %         r2n11n_memoryPrecision_valid_size(~temp1) = temp_size_min;
        %
        %         temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
        %         temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow);
        %         temptempBoolIndex12 = temptempBoolIndex1 | temptempBoolIndex2;
        
        for tempj=1:temp_roiNum
            %if temptempBoolIndex12(tempj) == false
            %    continue
            %end
            
            
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_memoryPrecision_valid_size(tempj);
            
            %if temp_size < temp_size_threshold
            %    %continue
            %    temp_size = temp_size_min;
            %end
            if temp_size <= temp_size_min
                continue
            end
            
            %if cellType(tempj) == 3
            %    temptemp_roiColor = temp_roiColor_mix;
            %else
            %    temptemp_roiColor = temp_roiColor_purePrecision;
            %end
            temptemp_roiColor = temp_roiColor_purePrecision;
            
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %         scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*1);
            %         scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            
            %         if temptempBoolIndex12(tempj) == true
            %             scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %         else
            %             scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            %             %if temp_size > 30
            %                 a = 1;
            %             %end
            %         end
            hold on
            
        end
        
        y_min = 0;
        y_max = 511;
        
        x_min = 0;
        x_max = 511;
        
        
        %         temp1 = centriod_memoryPrecision_all;
        %         plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*0.08 temp1(2)+(y_max-y_min)*0.08],...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
        %         hold on
        %         plot([temp1(1)-(y_max-y_min)*0.08 temp1(1)+(y_max-y_min)*0.08],[1 1].*temp1(2),...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
        %         hold on
        
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'Visible','off');
        if if_tuning_location0_precision1 == 1
            temp_title = title(sprintf('Memory precision'),'FontSize',9);
        elseif if_tuning_location0_precision1 == 0
            if if_r2Binary_unitView0_subspaceView1 == 0
                temp_title = title(sprintf('WM'),'FontSize',9);
            elseif if_r2Binary_unitView0_subspaceView1 == 1
                temp_title = title(sprintf('WM strength'),'FontSize',9);
            end
        end
        temp_title.Visible = 'on';
        
        %% Plot B, meta-memory
        nexttile
        
        image(temp_template2);
        hold on
        colormap(gray);
        
        %         r2n11n_choiceMemory_valid_size = r2n11n_choiceMemory_valid.^1;%0.00001
        %         %r2n11n_choiceMemory_valid_size = rescale(r2n11n_choiceMemory_valid_size,temp_size_min,temp_size_max*(max_r2_choiceMemory_valid/max_r2_memoryPrecision_valid));
        %         r2n11n_choiceMemory_valid_size = rescale(r2n11n_choiceMemory_valid_size,temp_size_min,temp_size_max);
        %
        %
        %         temp_size_threshold = prctile(r2n11n_choiceMemory_valid_size,temp_size_threshold_prctile);
        %         temp1 = r2n11n_choiceMemory_valid_size>temp_size_threshold;
        %         temp2 = r2n11n_choiceMemory_valid_size(temp1);
        %         temp2_n11n = rescale(temp2,temp_size_min,temp_size_max);
        %         r2n11n_choiceMemory_valid_size(temp1) = temp2_n11n;
        %         r2n11n_choiceMemory_valid_size(~temp1) = temp_size_min;
        %
        %
        %         temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
        %         temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh);
        %         temptempBoolIndex12 = temptempBoolIndex1 | temptempBoolIndex2;
        
        for tempj=1:temp_roiNum
            %if temptempBoolIndex12(tempj) == false
            %    continue
            %end
            
            
            
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_valid_size(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            
            %if cellType(tempj) == 3
            %    temptemp_roiColor = temp_roiColor_mix;
            %else
            %    temptemp_roiColor = temp_roiColor_pureMeta;
            %end
            temptemp_roiColor = temp_roiColor_pureMeta;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %         scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*1);
            %         scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
            
        end
        
        y_min = 0;
        y_max = 511;
        
        x_min = 0;
        x_max = 511;
        
        
        %         temp1 = centriod_meta_all;
        %         plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*0.08 temp1(2)+(y_max-y_min)*0.08],...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
        %         hold on
        %         plot([temp1(1)-(y_max-y_min)*0.08 temp1(1)+(y_max-y_min)*0.08],[1 1].*temp1(2),...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
        %         hold on
        
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        % set(gca, 'YDir','reverse');
        set(gca,'Visible','off');
        temp_title = title(sprintf('Meta-WM'),'FontSize',9);
        temp_title.Visible = 'on';
        
        
        %% Plot C, baseline meta-memory
        nexttile
        
        image(temp_template2);
        hold on
        colormap(gray);
        
        %         r2n11n_choiceMemory_baseline_valid_size = r2n11n_choiceMemory_baseline_valid.^1;%0.00001
        %         %r2n11n_choiceMemory_baseline_valid_size = rescale(r2n11n_choiceMemory_baseline_valid_size,temp_size_min,temp_size_max*(max_r2_choiceMemory_baseline_valid/max_r2_memoryPrecision_valid));
        %         r2n11n_choiceMemory_baseline_valid_size = rescale(r2n11n_choiceMemory_baseline_valid_size,temp_size_min,temp_size_max);
        %
        %
        %         temp_size_threshold = prctile(r2n11n_choiceMemory_baseline_valid_size,temp_size_threshold_prctile);
        %         temp1 = r2n11n_choiceMemory_baseline_valid_size>temp_size_threshold;
        %         temp2 = r2n11n_choiceMemory_baseline_valid_size(temp1);
        %         temp2_n11n = rescale(temp2,temp_size_min,temp_size_max);
        %         r2n11n_choiceMemory_baseline_valid_size(temp1) = temp2_n11n;
        %         r2n11n_choiceMemory_baseline_valid_size(~temp1) = temp_size_min;
        
        
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_baseline_valid_size(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            
            temptemp_roiColor = temp_roiColor_mix;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            hold on
            
        end
        
        y_min = 0;
        y_max = 511;
        
        x_min = 0;
        x_max = 511;
        
        
        %         temp1 = centriod_meta_baseline_all;
        %         plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*0.08 temp1(2)+(y_max-y_min)*0.08],...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
        %         hold on
        %         plot([temp1(1)-(y_max-y_min)*0.08 temp1(1)+(y_max-y_min)*0.08],[1 1].*temp1(2),...
        %             'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
        %         hold on
        
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        % set(gca, 'YDir','reverse');
        set(gca,'Visible','off');
        %temp_title = title(sprintf('Baseline Meta'),'FontSize',9);
        temp_title = title(sprintf('Baseline'),'FontSize',9);
        temp_title.Visible = 'on';
        
        
        
        %% Plot D, merged
        nexttile
        
        temp_template_onlyBorder = temp_template2;
        
        temp_template_onlyBorder(:,:) = 255;
        temp1 = 6+2;%10
        temp_template_onlyBorder(1:temp1,:) = 0;
        temp_template_onlyBorder(end+1-temp1:end,:) = 0;
        temp_template_onlyBorder(:,1:temp1) = 0;
        temp_template_onlyBorder(:,end+1-temp1:end) = 0;
        
        
        %image(temp_template2);
        image(temp_template_onlyBorder);
        hold on
        colormap(gray);
        
        
        
        %%
        r2n11n_memoryPrecision_valid_size_B = r2n11n_memoryPrecision_valid_size;
        r2n11n_choiceMemory_valid_size_B = r2n11n_choiceMemory_valid_size;
        r2n11n_choiceMemory_baseline_valid_size_B = r2n11n_choiceMemory_baseline_valid_size;
        
        for tempj=1:temp_roiNum
            %temp1 = r2n11n_memoryPrecision_valid_size_B(tempj);
            %temp2 = r2n11n_choiceMemory_valid_size_B(tempj);
            %temp3 = r2n11n_choiceMemory_baseline_valid_size_B(tempj);
            
            temp1 = r2n11n_memoryPrecision_raw(tempj);
            temp2 = r2n11n_choiceMemory_raw(tempj);
            temp3 = r2n11n_choiceMemory_baseline_raw(tempj);
            
            temp123 = [temp1,temp2,temp3];
            [M,I] = max(temp123);
            
            if I == 1
                r2n11n_choiceMemory_valid_size_B(tempj) = temp_size_min;
                r2n11n_choiceMemory_baseline_valid_size_B(tempj) = temp_size_min;
            elseif I == 2
                r2n11n_memoryPrecision_valid_size_B(tempj) = temp_size_min;
                r2n11n_choiceMemory_baseline_valid_size_B(tempj) = temp_size_min;
            elseif I == 3
                r2n11n_memoryPrecision_valid_size_B(tempj) = temp_size_min;
                r2n11n_choiceMemory_valid_size_B(tempj) = temp_size_min;
            end
        end
        
        %% Memory precision
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_memoryPrecision_valid_size_B(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_purePrecision;
            
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
            
        end
        
        
        %% Meta-memory
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_valid_size_B(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_pureMeta;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
            
        end
        
        
        
        %% Baseline meta
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_baseline_valid_size_B(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_mix;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
        end
        
        
        
        %%
        
        temptemp_coeff = 0.08;
        temptemp_coeff_B = 0.1;
        
        temp1 = centriod_memoryPrecision_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff_B temp1(2)+(y_max-y_min)*temptemp_coeff_B],...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff_B temp1(1)+(y_max-y_min)*temptemp_coeff_B],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        temp1 = centriod_meta_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff_B temp1(1)+(y_max-y_min)*temptemp_coeff_B],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff_B temp1(2)+(y_max-y_min)*temptemp_coeff_B],...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        temp1 = centriod_meta_baseline_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff_B temp1(1)+(y_max-y_min)*temptemp_coeff_B],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff_B temp1(2)+(y_max-y_min)*temptemp_coeff_B],...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        
        
        
        temp1 = centriod_memoryPrecision_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff temp1(2)+(y_max-y_min)*temptemp_coeff],...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_meta_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff temp1(2)+(y_max-y_min)*temptemp_coeff],...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_meta_baseline_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff temp1(2)+(y_max-y_min)*temptemp_coeff],...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
        hold on
        
        
        
        
        temp1 = centriod_meta_baseline_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff temp1(1)+(y_max-y_min)*temptemp_coeff],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_meta_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff temp1(1)+(y_max-y_min)*temptemp_coeff],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_memoryPrecision_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff temp1(1)+(y_max-y_min)*temptemp_coeff],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
        hold on
        
        
        
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'Visible','off');        
        
        temp_title = title(sprintf('Merged'),'FontSize',9);
        temp_title.Visible = 'on';
        
    end
    
    
    %% Plot 3 feature and merged map.
    if true
        fig = figure('Name',' ','NumberTitle','off');
        set(gcf,'Position',[700 100 720*0.7*0.93*0.9*0.95*1.2*1.05 300*0.7*0.93*0.9*0.95*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
        
        temp_roiNum = length(tempMappingCellIndex_suite2p);
        
        
        %% Plot A, baseline meta-memory
        nexttile
        
        image(temp_template2);
        hold on
        colormap(gray);
                
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_baseline_valid_size(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            
            temptemp_roiColor = temp_roiColor_mix;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            hold on
            
        end
        
        y_min = 0;
        y_max = 511;
        
        x_min = 0;
        x_max = 511;
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        % set(gca, 'YDir','reverse');
        set(gca,'Visible','off');
        %temp_title = title(sprintf('Baseline Meta'),'FontSize',9);
        temp_title = title(sprintf('Baseline'),'FontSize',9);
        temp_title.Visible = 'on';
        
        
        %% Plot B, memoryPrecision
        nexttile
        
        image(temp_template2);
        hold on
        colormap(gray);
                
        for tempj=1:temp_roiNum            
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_memoryPrecision_valid_size(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_purePrecision;
            
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            hold on
            
        end
        
        y_min = 0;
        y_max = 511;
        
        x_min = 0;
        x_max = 511;
                       
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'Visible','off');
        if if_tuning_location0_precision1 == 1
            temp_title = title(sprintf('Memory precision'),'FontSize',9);
        elseif if_tuning_location0_precision1 == 0
            %temp_title = title(sprintf('WM'),'FontSize',9);
            if if_r2Binary_unitView0_subspaceView1 == 0
                temp_title = title(sprintf('WM'),'FontSize',9);
            elseif if_r2Binary_unitView0_subspaceView1 == 1
                temp_title = title(sprintf('WM strength'),'FontSize',9);
            end            
        end
        temp_title.Visible = 'on';
        
        %% Plot C, meta-memory
        nexttile
        
        image(temp_template2);
        hold on
        colormap(gray);
                
        for tempj=1:temp_roiNum                    
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_valid_size(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_pureMeta;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            hold on
            
        end
        
        y_min = 0;
        y_max = 511;
        
        x_min = 0;
        x_max = 511;
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        % set(gca, 'YDir','reverse');
        set(gca,'Visible','off');
        temp_title = title(sprintf('Meta-WM'),'FontSize',9);
        temp_title.Visible = 'on';
        
        
        
        
        
        %% Plot D, merged
        nexttile
        
        temp_template_onlyBorder = temp_template2;
        
        temp_template_onlyBorder(:,:) = 255;
        temp1 = 6+2;%10
        temp_template_onlyBorder(1:temp1,:) = 0;
        temp_template_onlyBorder(end+1-temp1:end,:) = 0;
        temp_template_onlyBorder(:,1:temp1) = 0;
        temp_template_onlyBorder(:,end+1-temp1:end) = 0;
        
        
        %image(temp_template2);
        image(temp_template_onlyBorder);
        hold on
        colormap(gray);
        
        
        
        %%
        r2n11n_memoryPrecision_valid_size_B = r2n11n_memoryPrecision_valid_size;
        r2n11n_choiceMemory_valid_size_B = r2n11n_choiceMemory_valid_size;
        r2n11n_choiceMemory_baseline_valid_size_B = r2n11n_choiceMemory_baseline_valid_size;
        
        for tempj=1:temp_roiNum            
            temp1 = r2n11n_memoryPrecision_raw(tempj);
            temp2 = r2n11n_choiceMemory_raw(tempj);
            temp3 = r2n11n_choiceMemory_baseline_raw(tempj);
            
            temp123 = [temp1,temp2,temp3];
            [M,I] = max(temp123);
            
            if I == 1
                r2n11n_choiceMemory_valid_size_B(tempj) = temp_size_min;
                r2n11n_choiceMemory_baseline_valid_size_B(tempj) = temp_size_min;
            elseif I == 2
                r2n11n_memoryPrecision_valid_size_B(tempj) = temp_size_min;
                r2n11n_choiceMemory_baseline_valid_size_B(tempj) = temp_size_min;
            elseif I == 3
                r2n11n_memoryPrecision_valid_size_B(tempj) = temp_size_min;
                r2n11n_choiceMemory_valid_size_B(tempj) = temp_size_min;
            end
        end
        
        %% Memory precision
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_memoryPrecision_valid_size_B(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_purePrecision;
            
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
            
        end
        
        
        %% Meta-memory
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_valid_size_B(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_pureMeta;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
            
        end
        
        
        
        %% Baseline meta
        for tempj=1:temp_roiNum
            temp_cellIndex_suite2p = tempMappingCellIndex_suite2p(tempj);
            
            temp_roi_stat = roi_stats{cellIndex_suite2p==temp_cellIndex_suite2p};
            temp_roi_xpix = double(temp_roi_stat.xpix);
            temp_roi_ypix = double(temp_roi_stat.ypix);
            temp_roi_med = round([median(temp_roi_ypix),median(temp_roi_xpix)]);
            
            temp_size = r2n11n_choiceMemory_baseline_valid_size_B(tempj);
            
            if temp_size <= temp_size_min
                continue
            end
            
            temptemp_roiColor = temp_roiColor_mix;
            
            scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled','MarkerEdgeColor',[1 1 1]*0);
            %scatter(temp_roi_med(2),temp_roi_med(1),temp_size,temptemp_roiColor,'filled');
            hold on
        end
        
        
        
        %%
        
        temptemp_coeff = 0.08;
        temptemp_coeff_B = 0.1;
        
        temp1 = centriod_memoryPrecision_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff_B temp1(2)+(y_max-y_min)*temptemp_coeff_B],...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff_B temp1(1)+(y_max-y_min)*temptemp_coeff_B],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        temp1 = centriod_meta_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff_B temp1(1)+(y_max-y_min)*temptemp_coeff_B],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff_B temp1(2)+(y_max-y_min)*temptemp_coeff_B],...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        temp1 = centriod_meta_baseline_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff_B temp1(1)+(y_max-y_min)*temptemp_coeff_B],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff_B temp1(2)+(y_max-y_min)*temptemp_coeff_B],...
            'LineWidth',temp_lineWidthB,'color',[0 0 0]);%[0 0.4470 0.7410]
        hold on
        
        
        
        temp1 = centriod_memoryPrecision_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff temp1(2)+(y_max-y_min)*temptemp_coeff],...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_meta_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff temp1(2)+(y_max-y_min)*temptemp_coeff],...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_meta_baseline_all;
        plot([1 1].*temp1(1),[temp1(2)-(y_max-y_min)*temptemp_coeff temp1(2)+(y_max-y_min)*temptemp_coeff],...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
        hold on
        
        
        
        
        temp1 = centriod_meta_baseline_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff temp1(1)+(y_max-y_min)*temptemp_coeff],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_mix);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_meta_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff temp1(1)+(y_max-y_min)*temptemp_coeff],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_pureMeta);%[0 0.4470 0.7410]
        hold on
        
        temp1 = centriod_memoryPrecision_all;
        plot([temp1(1)-(y_max-y_min)*temptemp_coeff temp1(1)+(y_max-y_min)*temptemp_coeff],[1 1].*temp1(2),...
            'LineWidth',temp_lineWidth,'color',temp_roiColor_purePrecision);%[0 0.4470 0.7410]
        hold on
        
        
        
        
        axis equal;
        xlim([x_min x_max])
        ylim([y_min y_max])
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'Visible','off');        
        
        temp_title = title(sprintf('Merged'),'FontSize',9);
        temp_title.Visible = 'on';
        
    end
    
    
    
    %% Plot centriod distance distribution
    if false
        fig = figure('Name',' ','NumberTitle','off');
        %set(gcf,'Position',[1200 120 240 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[1500 120 240*1.1 260*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[700 500 720 240]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[700 500 720 240*0.85*0.9*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        h = histogram(centriodDis_all_shuffled,'Normalization','pdf');
        h.FaceColor = [1 1 1]*0.75;
        hold on
        
        [y_min,y_max] = bounds(h.Values);
        
        plot([1 1].*centriodDis_all_shuffled_prctileA,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
            'LineWidth',2,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        plot([1 1].*centriodDis_all_shuffled_prctileB,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
            'LineWidth',2,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        
        plot([1 1].*centriodDis_all,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
            'LineWidth',2,'color',[0.75 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        
        x_min = min([centriodDis_all_shuffled_prctileA,centriodDis_all_shuffled_prctileB,centriodDis_all]);
        x_max = max([centriodDis_all_shuffled_prctileA,centriodDis_all_shuffled_prctileB,centriodDis_all]);
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.15 x_max+(x_max-x_min)*0.15]);
        ylim([y_min-(y_max-y_min)*0.00 y_max+(y_max-y_min)*0.00]);
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Centriod distance', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Pdf', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Precision and Meta'),'FontSize',11);
        temp_title.Interpreter = 'none';
        
        
        
        nexttile
        
        h = histogram(centriodDis_all_AC_shuffled,'Normalization','pdf');
        h.FaceColor = [1 1 1]*0.75;
        hold on
        
        [y_min,y_max] = bounds(h.Values);
        
        plot([1 1].*centriodDis_all_AC_shuffled_prctileA,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
            'LineWidth',2,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        plot([1 1].*centriodDis_all_AC_shuffled_prctileB,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
            'LineWidth',2,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        
        plot([1 1].*centriodDis_AC_all,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
            'LineWidth',2,'color',[0.75 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        
        x_min = min([centriodDis_all_AC_shuffled_prctileA,centriodDis_all_AC_shuffled_prctileB,centriodDis_AC_all]);
        x_max = max([centriodDis_all_AC_shuffled_prctileA,centriodDis_all_AC_shuffled_prctileB,centriodDis_AC_all]);
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.15 x_max+(x_max-x_min)*0.15]);
        ylim([y_min-(y_max-y_min)*0.00 y_max+(y_max-y_min)*0.00]);
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Centriod distance', 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel('Pdf', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Precision and Baseline'),'FontSize',11);
        temp_title.Interpreter = 'none';
        
        
        
        nexttile
        
        h = histogram(centriodDis_all_BC_shuffled,'Normalization','pdf');
        h.FaceColor = [1 1 1]*0.75;
        hold on
        
        [y_min,y_max] = bounds(h.Values);
        
        plot([1 1].*centriodDis_all_BC_shuffled_prctileA,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
            'LineWidth',2,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        plot([1 1].*centriodDis_all_BC_shuffled_prctileB,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
            'LineWidth',2,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        
        plot([1 1].*centriodDis_BC_all,[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
            'LineWidth',2,'color',[0.75 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        
        x_min = min([centriodDis_all_BC_shuffled_prctileA,centriodDis_all_BC_shuffled_prctileB,centriodDis_BC_all]);
        x_max = max([centriodDis_all_BC_shuffled_prctileA,centriodDis_all_BC_shuffled_prctileB,centriodDis_BC_all]);
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.15 x_max+(x_max-x_min)*0.15]);
        ylim([y_min-(y_max-y_min)*0.00 y_max+(y_max-y_min)*0.00]);
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Centriod distance', 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel('Pdf', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Meta and Baseline'),'FontSize',11);
        temp_title.Interpreter = 'none';
        
    end
    
    
    %% Plot centriod distance distribution in one panel
    if true
        fig = figure('Name','Centroid distance','NumberTitle','off');
        %set(gcf,'Position',[50+0 50+0 240*0.80 336*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[50+0 50+0 240*0.80 336*1.11*0.67*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_1 = centriodDis_all_shuffled';
        temp_2 = centriodDis_all_AC_shuffled';
        temp_3 = centriodDis_all_BC_shuffled';
        
        temp_x1_mean = centriodDis_all;
        temp_x2_mean = centriodDis_AC_all;
        temp_x3_mean = centriodDis_BC_all;
        
        %[~,temp_p_1] = ttest(temp_1,temp_x1_mean);
        %[~,temp_p_2] = ttest(temp_2,temp_x2_mean);
        %[~,temp_p_3] = ttest(temp_3,temp_x3_mean);
        
        temp_1_sorted = sort(temp_1);
        temp_p_1 = sum(temp_x1_mean<temp_1_sorted)/length(temp_1_sorted);
        
        temp_2_sorted = sort(temp_2);
        temp_p_2 = sum(temp_x2_mean<temp_2_sorted)/length(temp_2_sorted);
        
        temp_3_sorted = sort(temp_3);
        temp_p_3 = sum(temp_x3_mean<temp_3_sorted)/length(temp_3_sorted);
        
        
        temp_y_min = min([temp_1 temp_2 temp_3 temp_x1_mean temp_x2_mean temp_x3_mean]);
        temp_y_max = max([temp_1 temp_2 temp_3 temp_x1_mean temp_x2_mean temp_x3_mean]);
        
        temp_data = [temp_1';temp_2';temp_3'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        
        temp_label = [g1;g2;g3];
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 3, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        h(3).ViolinPlot.FaceAlpha = 0.1;
        
        
        plot([1-0.45 1+0.45],temp_x1_mean*[1 1],'color',[0.6350 0.0780 0.1840],'linewidth',4);
        hold on
        
        plot([2-0.45 2+0.45],temp_x2_mean*[1 1],'color',[0.6350 0.0780 0.1840],'linewidth',4);
        hold on
        
        plot([3-0.45 3+0.45],temp_x3_mean*[1 1],'color',[0.6350 0.0780 0.1840],'linewidth',4);
        hold on
        
        tempTxt = sprintf('');
        if temp_p_1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p_1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p_1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y_min-(temp_y_max-temp_y_min)*0.15,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        tempTxt = sprintf('');
        if temp_p_2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p_2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p_2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y_min-(temp_y_max-temp_y_min)*0.15,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        %plot([2-0.2 2+0.2],temp_y_max+(temp_y_max-temp_y_min)*0.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        %hold on
        
        tempTxt = sprintf('');
        if temp_p_3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p_3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p_3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y_min-(temp_y_max-temp_y_min)*0.15,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        %plot([3-0.2 3+0.2],temp_y_max+(temp_y_max-temp_y_min)*0.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
        %hold on
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0 4]);
        %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.1]);%0.3
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.20 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
        set(gca, 'FontSize', 8) %14
        %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
        set(gca,'box','off');% 取消右、上边框
        
        
        %xtl = ["Precision", "Meta", "Baseline"];
        if if_tuning_location0_precision1 == 1
            xtl = ["Q & Meta", "Precision & Baseline", "Meta & Baseline"];
        elseif if_tuning_location0_precision1 == 0
            xtl = ["WM & Meta", "WM & Baseline", "Meta & Baseline"];
        end
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.25;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
        
        set(gca,'xticklabel','');
        
        yticks([0 40 80]);
        
        ylabel('Centroid distance', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('Centroid distance'),'fontsize',9);
        %temp_title.Interpreter = 'none';
        
    end
    
end

%% Clustering analysis A
if if_plot_clusterA == 1
    KList = 1:30;
    % temp_shuffleNum_B = 100;%1000
    
    criterion = 'CalinskiHarabasz';
    % criterion = 'DaviesBouldin';
    % criterion = 'gap';
    % criterion = 'silhouette';
    
    
    %threshold_prctile_clusterA = 95;%95
    
    temp_t0 = tic;
    
    % memoryPrecision
    temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow);
    temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
    temptempBoolIndex12 = temptempBoolIndex1 | temptempBoolIndex2;
    eva_memoryPrecision = evalclusters(roi_med(temptempBoolIndex12,:),'kmeans',criterion,'KList',KList);
    CriterionValues_memoryPrecision = eva_memoryPrecision.CriterionValues;
    
    CriterionValues_memoryPrecision_shuffled = nan(length(KList),temp_shuffleNum_A);
    % for tempi=1:temp_shuffleNum_A
    parfor tempi=1:temp_shuffleNum_A
        temp_roi_med_shuffled = roi_med_shuffled(temptempBoolIndex12,:,tempi);
        temp_eva_memoryPrecision_shuffled = evalclusters(temp_roi_med_shuffled,'kmeans',criterion,'KList',KList);
        CriterionValues_memoryPrecision_shuffled(:,tempi) = temp_eva_memoryPrecision_shuffled.CriterionValues;
    end
    % fprintf('t_cluster_shuffled=%.3f\n',toc(temp_t0));
    
    temp_prctile = prctile(CriterionValues_memoryPrecision_shuffled,threshold_prctile_clusterA,2)';
    CriterionValues_memoryPrecision_sigBoolIndex = CriterionValues_memoryPrecision>temp_prctile;
    
    
    
    
    % choiceMemory
    temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
    temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh);
    temptempBoolIndex12 = temptempBoolIndex1 | temptempBoolIndex2;
    eva_choiceMemory = evalclusters(roi_med(temptempBoolIndex12,:),'kmeans',criterion,'KList',KList);
    CriterionValues_choiceMemory = eva_choiceMemory.CriterionValues;
    
    
    % temp_t0 = tic;
    CriterionValues_choiceMemory_shuffled = nan(length(KList),temp_shuffleNum_A);
    % for tempi=1:temp_shuffleNum_A
    parfor tempi=1:temp_shuffleNum_A
        temp_roi_med_shuffled = roi_med_shuffled(temptempBoolIndex12,:,tempi);
        temp_eva_choiceMemory_shuffled = evalclusters(temp_roi_med_shuffled,'kmeans',criterion,'KList',KList);
        CriterionValues_choiceMemory_shuffled(:,tempi) = temp_eva_choiceMemory_shuffled.CriterionValues;
    end
    
    temp_prctile = prctile(CriterionValues_choiceMemory_shuffled,threshold_prctile_clusterA,2)';
    CriterionValues_choiceMemory_sigBoolIndex = CriterionValues_choiceMemory>temp_prctile;
    
    
    fprintf('t_clusterA_shuffled=%.3f s\n',toc(temp_t0));
    
end


%% Clustering analysis B
if if_plot_clusterB == 1 && if_computeShuffle == 1
    % dis
    dis  = nan(temp_roiNum,temp_roiNum);
    for tempi=1:temp_roiNum
        temp_medA = roi_med(tempi,:);
        for tempj=1:temp_roiNum
            temp_medB = roi_med(tempj,:);
            temp_dis = norm(temp_medA-temp_medB);
            %dis(tempi,tempj) = temp_dis;
            dis(tempi,tempj) = temp_dis*700/512;
        end
    end
    
    if_smallToZero = 0;
    smallToZero_threshold = 0.002;%0.01,0.005,0.001,eps
    smallToZero_target = 10;%0.0001
    
    % r2_diff_memoryPrecision
    r2_diff_memoryPrecision = nan(temp_roiNum,temp_roiNum);
    for tempi=1:temp_roiNum
        temp_r2A = r2n11n_memoryPrecision_raw2(tempi);
        for tempj=1:temp_roiNum
            temp_r2B = r2n11n_memoryPrecision_raw2(tempj);
            %             if temp_r2A < 0
            %                 temp_r2A = 0;
            %             end
            %             if temp_r2B < 0
            %                 temp_r2B = 0;
            %             end
            temp_dff = abs((temp_r2A-temp_r2B)/(temp_r2A+temp_r2B));
            
            if if_smallToZero == 1
                if (abs(temp_r2A)+abs(temp_r2B)) < smallToZero_threshold
                    temp_dff = smallToZero_target;
                end
                %if temp_r2A < smallToZero_threshold || temp_r2B < smallToZero_threshold
                %   temp_dff = smallToZero_target;
                %end
            end
            
            r2_diff_memoryPrecision(tempi,tempj) = temp_dff;
        end
    end
    r2_diff_memoryPrecision(r2_diff_memoryPrecision>=1) = nan;
    
    % r2_diff_choiceMemory
    r2_diff_choiceMemory = nan(temp_roiNum,temp_roiNum);
    for tempi=1:temp_roiNum
        temp_r2A = r2n11n_choiceMemory_raw2(tempi);
        for tempj=1:temp_roiNum
            temp_r2B = r2n11n_choiceMemory_raw2(tempj);
            %             if temp_r2A < 0
            %                 temp_r2A = 0;
            %             end
            %             if temp_r2B < 0
            %                 temp_r2B = 0;
            %             end
            
            temp_dff = abs((temp_r2A-temp_r2B)/(temp_r2A+temp_r2B));
            
            if if_smallToZero == 1
                if (abs(temp_r2A)+abs(temp_r2B)) < smallToZero_threshold
                    temp_dff = smallToZero_target;
                end
                %if temp_r2A < smallToZero_threshold || temp_r2B < smallToZero_threshold
                %   temp_dff = smallToZero_target;
                %end
            end
            
            r2_diff_choiceMemory(tempi,tempj) = temp_dff;
        end
    end
    r2_diff_choiceMemory(r2_diff_choiceMemory>=1) = nan;
    
    % r2_diff_choiceMemory_baseline
    r2_diff_choiceMemory_baseline = nan(temp_roiNum,temp_roiNum);
    for tempi=1:temp_roiNum
        temp_r2A = r2n11n_choiceMemory_baseline_raw2(tempi);
        for tempj=1:temp_roiNum
            temp_r2B = r2n11n_choiceMemory_baseline_raw2(tempj);
            %             if temp_r2A < 0
            %                 temp_r2A = 0;
            %             end
            %             if temp_r2B < 0
            %                 temp_r2B = 0;
            %             end
            
            temp_dff = abs((temp_r2A-temp_r2B)/(temp_r2A+temp_r2B));
            
            if if_smallToZero == 1
                if (abs(temp_r2A)+abs(temp_r2B)) < smallToZero_threshold
                    temp_dff = smallToZero_target;
                end
                %if temp_r2A < smallToZero_threshold || temp_r2B < smallToZero_threshold
                %   temp_dff = smallToZero_target;
                %end
            end
            
            r2_diff_choiceMemory_baseline(tempi,tempj) = temp_dff;
        end
    end
    r2_diff_choiceMemory_baseline(r2_diff_choiceMemory_baseline>=1) = nan;
    
    
    %     dis_range = ...
    %         [0 50;
    %         %51 100;
    %         100 150;
    %         %151 200;
    %         200 250;
    %         %251 300;
    %         300 350;
    %         %351 400;
    %         400 450;
    %         %451 500;
    %         500 550];
    
    %     dis_range = ...
    %         [0 50;
    %         51 100;
    %         100 150;
    %         151 200;
    %         200 250;
    %         251 300;
    %         300 350;
    %         351 400;
    %         400 450;
    %         451 500;
    %         500 550];
    
    %     dis_range = ...
    %         [0 100;
    %         100 200;
    %         200 300;
    %         300 400;
    %         400 500;
    %         500 600;
    %         600 700;
    %         700 800];
    
    %     dis_range = ...
    %         [0 150;
    %         150 300;
    %         300 450;
    %         450 600;
    %         600 750];
    
    %     dis_range = ...
    %         [0 200;
    %         200 400;
    %         400 600;
    %         600 800];
    
    %     dis_range = ...
    %         [0 200;
    %         100 300;
    %         300 500;
    %         400 600;
    %         500 700];
    %     dis_range = dis_range + 50;
    
    %     dis_range = ...
    %         [0 100;
    %         100 200;
    %         300 400;
    %         400 500;
    %         500 600];
    
    dis_range = ...
        [0 200;
        150 350;
        300 500;
        450 650;
        600 800];
    
    
    
    r2_diff_memoryPrecision_disBin = cell(size(dis_range,1),1);
    r2_diff_choiceMemory_disBin = cell(size(dis_range,1),1);
    r2_diff_choiceMemory_baseline_disBin = cell(size(dis_range,1),1);
    clusterIndex_memoryPrecision_disBin_mean = nan(size(dis_range,1),1);
    clusterIndex_choiceMemory_disBin_mean = nan(size(dis_range,1),1);
    clusterIndex_choiceMemory_baseline_disBin_mean = nan(size(dis_range,1),1);
    clusterIndex_memoryPrecision_disBin_sem = nan(size(dis_range,1),1);
    clusterIndex_choiceMemory_disBin_sem = nan(size(dis_range,1),1);
    clusterIndex_choiceMemory_baseline_disBin_sem = nan(size(dis_range,1),1);
    
    for tempi=1:size(dis_range,1)
        temptempBoolIndex_disRange = (dis>(dis_range(tempi,1))) & (dis<(dis_range(tempi,2)));
        temp_pairNum = sum(temptempBoolIndex_disRange,'all');
        
        temp_r2_diff_memoryPrecision = r2_diff_memoryPrecision(temptempBoolIndex_disRange);
        temp_r2_diff_choiceMemory = r2_diff_choiceMemory(temptempBoolIndex_disRange);
        temp_r2_diff_choiceMemory_baseline = r2_diff_choiceMemory_baseline(temptempBoolIndex_disRange);
        
        temp_r2_diff_memoryPrecision = temp_r2_diff_memoryPrecision(~isnan(temp_r2_diff_memoryPrecision));
        temp_r2_diff_choiceMemory = temp_r2_diff_choiceMemory(~isnan(temp_r2_diff_choiceMemory));
        temp_r2_diff_choiceMemory_baseline = temp_r2_diff_choiceMemory_baseline(~isnan(temp_r2_diff_choiceMemory_baseline));
        
        r2_diff_memoryPrecision_disBin{tempi} = temp_r2_diff_memoryPrecision;
        r2_diff_choiceMemory_disBin{tempi} = temp_r2_diff_choiceMemory;
        r2_diff_choiceMemory_baseline_disBin{tempi} = temp_r2_diff_choiceMemory_baseline;
        
        %         clusterIndex_memoryPrecision_disBin_mean(tempi) = 1./mean(temp_r2_diff_memoryPrecision);
        %         clusterIndex_choiceMemory_disBin_mean(tempi) = 1./mean(temp_r2_diff_choiceMemory);
        %         clusterIndex_memoryPrecision_disBin_sem(tempi) = 1./(std(temp_r2_diff_memoryPrecision)/sqrt(temp_pairNum));
        %         clusterIndex_choiceMemory_disBin_sem(tempi) = 1./(std(temp_r2_diff_choiceMemory)/sqrt(temp_pairNum));
        
        %         clusterIndex_memoryPrecision_disBin_mean(tempi) = 1./mean(temp_r2_diff_memoryPrecision);
        %         clusterIndex_choiceMemory_disBin_mean(tempi) = 1./mean(temp_r2_diff_choiceMemory);
        %         clusterIndex_choiceMemory_baseline_disBin_mean(tempi) = 1./mean(temp_r2_diff_choiceMemory_baseline);
        
        clusterIndex_memoryPrecision_disBin_mean(tempi) = 1./median(temp_r2_diff_memoryPrecision);
        clusterIndex_choiceMemory_disBin_mean(tempi) = 1./median(temp_r2_diff_choiceMemory);
        clusterIndex_choiceMemory_baseline_disBin_mean(tempi) = 1./median(temp_r2_diff_choiceMemory_baseline);
        clusterIndex_memoryPrecision_disBin_sem(tempi) = 1./(std(temp_r2_diff_memoryPrecision)/sqrt(temp_pairNum));
        clusterIndex_choiceMemory_disBin_sem(tempi) = 1./(std(temp_r2_diff_choiceMemory)/sqrt(temp_pairNum));
        clusterIndex_choiceMemory_baseline_disBin_sem(tempi) = 1./(std(temp_r2_diff_choiceMemory_baseline)/sqrt(temp_pairNum));
    end
    
    
    
    
    %     for tempi=1:size(dis_range,1)
    %         temptempBoolIndex_disRange = (dis>(dis_range(tempi,1))) & (dis<(dis_range(tempi,2)));
    %
    %         temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow);
    %         temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
    %         temptempBoolIndex12A = temptempBoolIndex1 | temptempBoolIndex2;
    %
    %         temptempBoolIndex_disRangeA = temptempBoolIndex_disRange;
    %         temptempBoolIndex_disRangeA(~temptempBoolIndex12A,:) = false;
    %         temptempBoolIndex_disRangeA(:,~temptempBoolIndex12A) = false;
    %
    %
    %         temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
    %         temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh);
    %         temptempBoolIndex12B = temptempBoolIndex1 | temptempBoolIndex2;
    %
    %         temptempBoolIndex_disRangeB = temptempBoolIndex_disRange;
    %         temptempBoolIndex_disRangeB(~temptempBoolIndex12B,:) = false;
    %         temptempBoolIndex_disRangeB(:,~temptempBoolIndex12B) = false;
    %
    %
    %         temptempBoolIndexC = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_choiceMemoryBaselineHigh);
    %
    %         temptempBoolIndex_disRangeC = temptempBoolIndex_disRange;
    %         temptempBoolIndex_disRangeC(~temptempBoolIndexC,:) = false;
    %         temptempBoolIndex_disRangeC(:,~temptempBoolIndexC) = false;
    %
    %
    %         temp_pairNumA = sum(temptempBoolIndex_disRangeA,'all');
    %         temp_pairNumB = sum(temptempBoolIndex_disRangeB,'all');
    %         temp_pairNumC = sum(temptempBoolIndex_disRangeC,'all');
    %
    %         temp_r2_diff_memoryPrecision = r2_diff_memoryPrecision(temptempBoolIndex_disRangeA);
    %         temp_r2_diff_choiceMemory = r2_diff_choiceMemory(temptempBoolIndex_disRangeB);
    %         temp_r2_diff_choiceMemory_baseline = r2_diff_choiceMemory_baseline(temptempBoolIndex_disRangeC);
    %
    %         r2_diff_memoryPrecision_disBin{tempi} = temp_r2_diff_memoryPrecision;
    %         r2_diff_choiceMemory_disBin{tempi} = temp_r2_diff_choiceMemory;
    %         r2_diff_choiceMemory_baseline_disBin{tempi} = temp_r2_diff_choiceMemory_baseline;
    %
    %         clusterIndex_memoryPrecision_disBin_mean(tempi) = 1./mean(temp_r2_diff_memoryPrecision);
    %         clusterIndex_choiceMemory_disBin_mean(tempi) = 1./mean(temp_r2_diff_choiceMemory);
    %         clusterIndex_choiceMemory_baseline_disBin_mean(tempi) = 1./mean(temp_r2_diff_choiceMemory_baseline);
    %
    %         clusterIndex_memoryPrecision_disBin_sem(tempi) = 1./(std(temp_r2_diff_memoryPrecision)/sqrt(temp_pairNumA));
    %         clusterIndex_choiceMemory_disBin_sem(tempi) = 1./(std(temp_r2_diff_choiceMemory)/sqrt(temp_pairNumB));
    %         clusterIndex_choiceMemory_baseline_disBin_sem(tempi) = 1./(std(temp_r2_diff_choiceMemory_baseline)/sqrt(temp_pairNumC));
    %
    % %         clusterIndex_memoryPrecision_disBin_mean(tempi) = mean(temp_r2_diff_memoryPrecision);
    % %         clusterIndex_choiceMemory_disBin_mean(tempi) = mean(temp_r2_diff_choiceMemory);
    % %         clusterIndex_memoryPrecision_disBin_sem(tempi) = (std(temp_r2_diff_memoryPrecision)/sqrt(temp_pairNumA));
    % %         clusterIndex_choiceMemory_disBin_sem(tempi) = (std(temp_r2_diff_choiceMemory)/sqrt(temp_pairNumB));
    %
    %     end
    
    
    
    % shuffled
    temp_t0 = tic;
    
    %r2_diff_memoryPrecision_disBin_shuffled = cell(size(dis_range,1),temp_shuffleNum_B);
    %r2_diff_choiceMemory_disBin_shuffled = cell(size(dis_range,1),temp_shuffleNum_B);
    clusterIndex_memoryPrecision_disBin_shuffled_mean = nan(size(dis_range,1),temp_shuffleNum_B);
    clusterIndex_choiceMemory_disBin_shuffled_mean = nan(size(dis_range,1),temp_shuffleNum_B);
    clusterIndex_choiceMemory_baseline_disBin_shuffled_mean = nan(size(dis_range,1),temp_shuffleNum_B);
    clusterIndex_memoryPrecision_disBin_shuffled_sem = nan(size(dis_range,1),temp_shuffleNum_B);
    clusterIndex_choiceMemory_disBin_shuffled_sem = nan(size(dis_range,1),temp_shuffleNum_B);
    clusterIndex_choiceMemory_baseline_disBin_shuffled_sem = nan(size(dis_range,1),temp_shuffleNum_B);
    
    dis_shuffled  = nan(temp_roiNum,temp_roiNum,temp_shuffleNum_B);
    % for tempShuffledIndex=1:temp_shuffleNum_B
    parfor tempShuffledIndex=1:temp_shuffleNum_B
        for tempi=1:temp_roiNum
            temp_medA = roi_med_shuffled(tempi,:,tempShuffledIndex);
            for tempj=1:temp_roiNum
                temp_medB = roi_med_shuffled(tempj,:,tempShuffledIndex); %#ok<*PFBNS>
                temp_dis = norm(temp_medA-temp_medB);
                %dis_shuffled(tempi,tempj,tempShuffledIndex) = temp_dis;
                dis_shuffled(tempi,tempj,tempShuffledIndex) = temp_dis*700/512;
            end
        end
    end
    
    % for tempShuffledIndex=1:temp_shuffleNum_B
    parfor tempShuffledIndex=1:temp_shuffleNum_B
        
        temp_dis_shuffled = dis_shuffled(:,:,tempShuffledIndex);
        
        temp_r2_diff_memoryPrecision_disBin_shuffled_mean = nan(size(dis_range,1),1);
        temp_r2_diff_choiceMemory_disBin_shuffled_mean = nan(size(dis_range,1),1);
        temp_r2_diff_choiceMemory_baseline_disBin_shuffled_mean = nan(size(dis_range,1),1);
        temp_r2_diff_memoryPrecision_disBin_shuffled_sem = nan(size(dis_range,1),1);
        temp_r2_diff_choiceMemory_disBin_shuffled_sem = nan(size(dis_range,1),1);
        temp_r2_diff_choiceMemory_baseline_disBin_shuffled_sem = nan(size(dis_range,1),1);
        
        for tempi=1:size(dis_range,1)
            temptempBoolIndex_disRange = (temp_dis_shuffled>(dis_range(tempi,1))) & (temp_dis_shuffled<(dis_range(tempi,2)));
            temp_pairNum = sum(temptempBoolIndex_disRange,'all');
            
            temp_r2_diff_memoryPrecision = r2_diff_memoryPrecision(temptempBoolIndex_disRange);
            temp_r2_diff_choiceMemory = r2_diff_choiceMemory(temptempBoolIndex_disRange);
            temp_r2_diff_choiceMemory_baseline = r2_diff_choiceMemory_baseline(temptempBoolIndex_disRange);
            
            temp_r2_diff_memoryPrecision = temp_r2_diff_memoryPrecision(~isnan(temp_r2_diff_memoryPrecision));
            temp_r2_diff_choiceMemory = temp_r2_diff_choiceMemory(~isnan(temp_r2_diff_choiceMemory));
            temp_r2_diff_choiceMemory_baseline = temp_r2_diff_choiceMemory_baseline(~isnan(temp_r2_diff_choiceMemory_baseline));
            
            %             temp_r2_diff_memoryPrecision_disBin_shuffled_mean(tempi) = mean(temp_r2_diff_memoryPrecision);
            %             temp_r2_diff_choiceMemory_disBin_shuffled_mean(tempi) = mean(temp_r2_diff_choiceMemory);
            %             temp_r2_diff_choiceMemory_baseline_disBin_shuffled_mean(tempi) = mean(temp_r2_diff_choiceMemory_baseline);
            temp_r2_diff_memoryPrecision_disBin_shuffled_mean(tempi) = median(temp_r2_diff_memoryPrecision);
            temp_r2_diff_choiceMemory_disBin_shuffled_mean(tempi) = median(temp_r2_diff_choiceMemory);
            temp_r2_diff_choiceMemory_baseline_disBin_shuffled_mean(tempi) = median(temp_r2_diff_choiceMemory_baseline);
            
            temp_r2_diff_memoryPrecision_disBin_shuffled_sem(tempi) = std(temp_r2_diff_memoryPrecision)/sqrt(temp_pairNum);
            temp_r2_diff_choiceMemory_disBin_shuffled_sem(tempi) = std(temp_r2_diff_choiceMemory)/sqrt(temp_pairNum);
            temp_r2_diff_choiceMemory_baseline_disBin_shuffled_sem(tempi) = std(temp_r2_diff_choiceMemory_baseline)/sqrt(temp_pairNum);
            
        end
        clusterIndex_memoryPrecision_disBin_shuffled_mean(:,tempShuffledIndex) = 1./temp_r2_diff_memoryPrecision_disBin_shuffled_mean;
        clusterIndex_choiceMemory_disBin_shuffled_mean(:,tempShuffledIndex) = 1./temp_r2_diff_choiceMemory_disBin_shuffled_mean;
        clusterIndex_choiceMemory_baseline_disBin_shuffled_mean(:,tempShuffledIndex) = 1./temp_r2_diff_choiceMemory_baseline_disBin_shuffled_mean;
        
        clusterIndex_memoryPrecision_disBin_shuffled_sem(:,tempShuffledIndex) = 1./temp_r2_diff_memoryPrecision_disBin_shuffled_sem;
        clusterIndex_choiceMemory_disBin_shuffled_sem(:,tempShuffledIndex) = 1./temp_r2_diff_choiceMemory_disBin_shuffled_sem;
        clusterIndex_choiceMemory_baseline_disBin_shuffled_sem(:,tempShuffledIndex) = 1./temp_r2_diff_choiceMemory_baseline_disBin_shuffled_sem;
        
        %         clusterIndex_memoryPrecision_disBin_shuffled_mean(:,tempShuffledIndex) = temp_r2_diff_memoryPrecision_disBin_shuffled_mean;
        %         clusterIndex_choiceMemory_disBin_shuffled_mean(:,tempShuffledIndex) = temp_r2_diff_choiceMemory_disBin_shuffled_mean;
        %         clusterIndex_memoryPrecision_disBin_shuffled_sem(:,tempShuffledIndex) = temp_r2_diff_memoryPrecision_disBin_shuffled_sem;
        %         clusterIndex_choiceMemory_disBin_shuffled_sem(:,tempShuffledIndex) = temp_r2_diff_choiceMemory_disBin_shuffled_sem;
        
    end
    
    
    
    
    %     % for tempShuffledIndex=1:temp_shuffleNum_B
    %     parfor tempShuffledIndex=1:temp_shuffleNum_B
    %
    %         temp_dis_shuffled = dis_shuffled(:,:,tempShuffledIndex);
    %
    %         temp_r2_diff_memoryPrecision_disBin_shuffled_mean = nan(size(dis_range,1),1);
    %         temp_r2_diff_choiceMemory_disBin_shuffled_mean = nan(size(dis_range,1),1);
    %         temp_r2_diff_choiceMemory_baseline_disBin_shuffled_mean = nan(size(dis_range,1),1);
    %         temp_r2_diff_memoryPrecision_disBin_shuffled_sem = nan(size(dis_range,1),1);
    %         temp_r2_diff_choiceMemory_disBin_shuffled_sem = nan(size(dis_range,1),1);
    %         temp_r2_diff_choiceMemory_baseline_disBin_shuffled_sem = nan(size(dis_range,1),1);
    %
    %         for tempi=1:size(dis_range,1)
    %             temptempBoolIndex_disRange = (temp_dis_shuffled>(dis_range(tempi,1))) & (temp_dis_shuffled<(dis_range(tempi,2)));
    %
    %             temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow);
    %             temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
    %             temptempBoolIndex12A = temptempBoolIndex1 | temptempBoolIndex2;
    %
    %             temptempBoolIndex_disRangeA = temptempBoolIndex_disRange;
    %             temptempBoolIndex_disRangeA(~temptempBoolIndex12A,:) = false;
    %             temptempBoolIndex_disRangeA(:,~temptempBoolIndex12A) = false;
    %
    %
    %             temptempBoolIndex1 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh);
    %             temptempBoolIndex2 = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh);
    %             temptempBoolIndex12B = temptempBoolIndex1 | temptempBoolIndex2;
    %
    %             temptempBoolIndex_disRangeB = temptempBoolIndex_disRange;
    %             temptempBoolIndex_disRangeB(~temptempBoolIndex12B,:) = false;
    %             temptempBoolIndex_disRangeB(:,~temptempBoolIndex12B) = false;
    %
    %             temptempBoolIndexC = ismember(tempMappingCellIndex_suite2p,cellIndex_suite2p_choiceMemoryBaselineHigh);
    %
    %             temptempBoolIndex_disRangeC = temptempBoolIndex_disRange;
    %             temptempBoolIndex_disRangeC(~temptempBoolIndexC,:) = false;
    %             temptempBoolIndex_disRangeC(:,~temptempBoolIndexC) = false;
    %
    %
    %             temp_pairNumA = sum(temptempBoolIndex_disRangeA,'all');
    %             temp_pairNumB = sum(temptempBoolIndex_disRangeB,'all');
    %             temp_pairNumC = sum(temptempBoolIndex_disRangeC,'all');
    %
    %             temp_r2_diff_memoryPrecision = r2_diff_memoryPrecision(temptempBoolIndex_disRangeA);
    %             temp_r2_diff_choiceMemory = r2_diff_choiceMemory(temptempBoolIndex_disRangeB);
    %             temp_r2_diff_choiceMemory_baseline = r2_diff_choiceMemory_baseline(temptempBoolIndex_disRangeC);
    %
    %
    %             temp_r2_diff_memoryPrecision_disBin_shuffled_mean(tempi) = mean(temp_r2_diff_memoryPrecision);
    %             temp_r2_diff_choiceMemory_disBin_shuffled_mean(tempi) = mean(temp_r2_diff_choiceMemory);
    %             temp_r2_diff_choiceMemory_baseline_disBin_shuffled_mean(tempi) = mean(temp_r2_diff_choiceMemory_baseline);
    %             temp_r2_diff_memoryPrecision_disBin_shuffled_sem(tempi) = std(temp_r2_diff_memoryPrecision)/sqrt(temp_pairNumA);
    %             temp_r2_diff_choiceMemory_disBin_shuffled_sem(tempi) = std(temp_r2_diff_choiceMemory)/sqrt(temp_pairNumB);
    %             temp_r2_diff_choiceMemory_baseline_disBin_shuffled_sem(tempi) = std(temp_r2_diff_choiceMemory_baseline)/sqrt(temp_pairNumC);
    %         end
    %
    %         %clusterIndex_memoryPrecision_disBin_shuffled_mean(:,tempShuffledIndex) = 1./temp_r2_diff_memoryPrecision_disBin_shuffled_mean;
    %         %clusterIndex_choiceMemory_disBin_shuffled_mean(:,tempShuffledIndex) = 1./temp_r2_diff_choiceMemory_disBin_shuffled_mean;
    %         %clusterIndex_memoryPrecision_disBin_shuffled_sem(:,tempShuffledIndex) = 1./temp_r2_diff_memoryPrecision_disBin_shuffled_sem;
    %         %clusterIndex_choiceMemory_disBin_shuffled_sem(:,tempShuffledIndex) = 1./temp_r2_diff_choiceMemory_disBin_shuffled_sem;
    %
    %         clusterIndex_memoryPrecision_disBin_shuffled_mean(:,tempShuffledIndex) = temp_r2_diff_memoryPrecision_disBin_shuffled_mean;
    %         clusterIndex_choiceMemory_disBin_shuffled_mean(:,tempShuffledIndex) = temp_r2_diff_choiceMemory_disBin_shuffled_mean;
    %         clusterIndex_choiceMemory_baseline_disBin_shuffled_mean(:,tempShuffledIndex) = temp_r2_diff_choiceMemory_baseline_disBin_shuffled_mean;
    %         clusterIndex_memoryPrecision_disBin_shuffled_sem(:,tempShuffledIndex) = temp_r2_diff_memoryPrecision_disBin_shuffled_sem;
    %         clusterIndex_choiceMemory_disBin_shuffled_sem(:,tempShuffledIndex) = temp_r2_diff_choiceMemory_disBin_shuffled_sem;
    %         clusterIndex_choiceMemory_baseline_disBin_shuffled_sem(:,tempShuffledIndex) = temp_r2_diff_choiceMemory_baseline_disBin_shuffled_sem;
    %     end
    
    
    
    
    fprintf('t_clusterB_shuffled=%.3f s\n',toc(temp_t0));
    
    clusterIndex_memoryPrecision_disBin_shuffled_meanB = mean(clusterIndex_memoryPrecision_disBin_shuffled_mean,2);
    clusterIndex_memoryPrecision_disBin_shuffled_semB = std(clusterIndex_memoryPrecision_disBin_shuffled_mean,1,2)./sqrt(temp_shuffleNum_B);
    clusterIndex_memoryPrecision_disBin_shuffled_stdB = std(clusterIndex_memoryPrecision_disBin_shuffled_mean,1,2);
    
    clusterIndex_choiceMemory_disBin_shuffled_meanB = mean(clusterIndex_choiceMemory_disBin_shuffled_mean,2);
    clusterIndex_choiceMemory_disBin_shuffled_semB = std(clusterIndex_choiceMemory_disBin_shuffled_mean,1,2)./sqrt(temp_shuffleNum_B);
    clusterIndex_choiceMemory_disBin_shuffled_stdB = std(clusterIndex_choiceMemory_disBin_shuffled_mean,1,2);
    
    clusterIndex_choiceMemory_baseline_disBin_shuffled_meanB = mean(clusterIndex_choiceMemory_baseline_disBin_shuffled_mean,2);
    clusterIndex_choiceMemory_baseline_disBin_shuffled_semB = std(clusterIndex_choiceMemory_baseline_disBin_shuffled_mean,1,2)./sqrt(temp_shuffleNum_B);
    clusterIndex_choiceMemory_baseline_disBin_shuffled_stdB = std(clusterIndex_choiceMemory_baseline_disBin_shuffled_mean,1,2);
    
    % temp_p = nan(size(dis_range,1),1);
    % for tempi=1:size(dis_range,1)
    %     [~,temp_p(tempi)] = ttest2(r2_diff_memoryPrecision_disBin{tempi},clusterIndex_memoryPrecision_disBin_shuffled_mean(tempi,:));
    % end
    
    
    % [~,temp_p] = ttest2(clusterIndex_memoryPrecision_disBin_shuffled_mean,);
    
    temp_prctile_high_memoryPrecision = prctile(clusterIndex_memoryPrecision_disBin_shuffled_mean,threshold_prctile_clusterB_high,2);
    temp_prctile_low_memoryPrecision = prctile(clusterIndex_memoryPrecision_disBin_shuffled_mean,threshold_prctile_clusterB_low,2);
    clusterIndex_memoryPrecision_sigBoolIndex_high = clusterIndex_memoryPrecision_disBin_mean>temp_prctile_high_memoryPrecision;
    clusterIndex_memoryPrecision_sigBoolIndex_low = clusterIndex_memoryPrecision_disBin_mean<temp_prctile_low_memoryPrecision;
    %clusterIndex_memoryPrecision_sigBoolIndex = clusterIndex_memoryPrecision_sigBoolIndex_high | clusterIndex_memoryPrecision_sigBoolIndex_low;
    clusterIndex_memoryPrecision_sigBoolIndex = clusterIndex_memoryPrecision_sigBoolIndex_high;
    
    temp_prctile_high_choiceMemory = prctile(clusterIndex_choiceMemory_disBin_shuffled_mean,threshold_prctile_clusterB_high,2);
    temp_prctile_low_choiceMemory = prctile(clusterIndex_choiceMemory_disBin_shuffled_mean,threshold_prctile_clusterB_low,2);
    clusterIndex_choiceMemory_sigBoolIndex_high = clusterIndex_choiceMemory_disBin_mean>temp_prctile_high_choiceMemory;
    clusterIndex_choiceMemory_sigBoolIndex_low = clusterIndex_choiceMemory_disBin_mean<temp_prctile_low_choiceMemory;
    %clusterIndex_choiceMemory_sigBoolIndex = clusterIndex_choiceMemory_sigBoolIndex_high | clusterIndex_choiceMemory_sigBoolIndex_low;
    clusterIndex_choiceMemory_sigBoolIndex = clusterIndex_choiceMemory_sigBoolIndex_high;
    
    temp_prctile_high_choiceMemory_baseline = prctile(clusterIndex_choiceMemory_baseline_disBin_shuffled_mean,threshold_prctile_clusterB_high,2);
    temp_prctile_low_choiceMemory_baseline = prctile(clusterIndex_choiceMemory_baseline_disBin_shuffled_mean,threshold_prctile_clusterB_low,2);
    clusterIndex_choiceMemory_baseline_sigBoolIndex_high = clusterIndex_choiceMemory_baseline_disBin_mean>temp_prctile_high_choiceMemory_baseline;
    clusterIndex_choiceMemory_baseline_sigBoolIndex_low = clusterIndex_choiceMemory_baseline_disBin_mean<temp_prctile_low_choiceMemory_baseline;
    %clusterIndex_choiceMemory_baseline_sigBoolIndex = clusterIndex_choiceMemory_baseline_sigBoolIndex_high | clusterIndex_choiceMemory_baseline_sigBoolIndex_low;
    clusterIndex_choiceMemory_baseline_sigBoolIndex = clusterIndex_choiceMemory_baseline_sigBoolIndex_high;
    
end


if if_plot_twoFeatureSplitB == 1
    
    if if_plot_clusterA == 1
        %% Plot centriod distance distribution
        fig = figure('Name',' ','NumberTitle','off');
        %set(gcf,'Position',[1200 120 240 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[600 520 264*2 260*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        
        
        %% memory precision
        nexttile
        
        h_line = [];
        
        x1 = KList(2:end);
        y1 = CriterionValues_memoryPrecision(2:end);
        h_line = [h_line plot(x1,y1,'LineWidth',1,'color',[0.25 0.25 0.25])];
        hold on
        
        x2 = KList(2:end);
        y2 = CriterionValues_memoryPrecision_shuffled(2:end,:);
        
        y2_mean = mean(y2,2)';
        y2_sem = (std(y2,1,2)./sqrt(size(y2,2)))';
        
        h_line = [h_line plot(x2,y2_mean,'LineWidth',1,'color',[0.75 0.75 0.75])];
        hold on
        
        patch([x2(:); flipud(x2(:))]', [y2_mean(:)+y2_sem(:); flipud(y2_mean(:)-y2_sem(:))]',[0.75 0.75 0.75],'FaceAlpha',0.05,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        x_min = min([x1,x2]);
        x_max = max([x1,x2]);
        y_min = min([y1,y2_mean-y2_sem]);
        y_max = max([y1,y2_mean+y2_sem]);
        
        scatter(x1(CriterionValues_memoryPrecision_sigBoolIndex(2:end)),y_max+(y_max-y_min)*0.08,20,[0 0 0],'*');
        hold on
        
        le = legend(h_line,'Data','Shuffled','Location','southeast','fontsize',8);%10
        le.ItemTokenSize = ones(1,2)*10;
        
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.15]);
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Cluster number', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Clustering scores', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Memory precision'),'FontSize',11);
        temp_title.Interpreter = 'none';
        
        
        
        %% meta-memory
        nexttile
        
        h_line = [];
        
        x1 = KList(2:end);
        y1 = CriterionValues_choiceMemory(2:end);
        h_line = [h_line plot(x1,y1,'LineWidth',1,'color',[0.25 0.25 0.25])];
        hold on
        
        x2 = KList(2:end);
        y2 = CriterionValues_choiceMemory_shuffled(2:end,:);
        
        y2_mean = mean(y2,2)';
        y2_sem = (std(y2,1,2)./sqrt(size(y2,2)))';
        
        h_line = [h_line plot(x2,y2_mean,'LineWidth',1,'color',[0.75 0.75 0.75])];
        hold on
        
        patch([x2(:); flipud(x2(:))]', [y2_mean(:)+y2_sem(:); flipud(y2_mean(:)-y2_sem(:))]',[0.75 0.75 0.75],'FaceAlpha',0.05,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        x_min = min([x1,x2]);
        x_max = max([x1,x2]);
        y_min = min([y1,y2_mean-y2_sem]);
        y_max = max([y1,y2_mean+y2_sem]);
        
        scatter(x1(CriterionValues_choiceMemory_sigBoolIndex(2:end)),y_max+(y_max-y_min)*0.08,20,[0 0 0],'*');
        hold on
        
        le = legend(h_line,'Data','Shuffled','Location','southeast','fontsize',8);%10
        le.ItemTokenSize = ones(1,2)*10;
        
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.15]);
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Cluster number', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Clustering scores', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Meta-memory'),'FontSize',11);
        temp_title.Interpreter = 'none';
        
    end
    
    
    
    
    
    %% Plot clustering analysis B
    if if_plot_clusterB == 1
        clusterIndex_memoryPrecision_disBin_mean;
        clusterIndex_memoryPrecision_disBin_sem;
        clusterIndex_memoryPrecision_disBin_shuffled_meanB;
        clusterIndex_memoryPrecision_disBin_shuffled_semB;
        
        clusterIndex_choiceMemory_disBin_mean;
        clusterIndex_choiceMemory_disBin_sem;
        clusterIndex_choiceMemory_disBin_shuffled_meanB;
        clusterIndex_choiceMemory_disBin_shuffled_semB;
        
        dis_range;
        
        fig = figure('Name',' ','NumberTitle','off');
        %set(gcf,'Position',[1200 120 240 260]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[600 520 264*2 260*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[50 520 720 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
        
        
        %% memory precision
        nexttile
        
        h_line = [];
        
        x1 = 1:size(dis_range,1);
        y1 = clusterIndex_memoryPrecision_disBin_mean;
        h_line = [h_line plot(x1,y1,'LineWidth',1,'color',[0.25 0.25 0.25])];
        hold on
        
        x2 = 1:size(dis_range,1);
        y2 = clusterIndex_memoryPrecision_disBin_shuffled_meanB;
        
        y2_sem = clusterIndex_memoryPrecision_disBin_shuffled_semB;
        %y2_sem = clusterIndex_memoryPrecision_disBin_shuffled_stdB;
        
        
        h_line = [h_line plot(x2,y2,'LineWidth',1,'color',[0.75 0.75 0.75])];
        hold on
        
        patch([x2(:); flipud(x2(:))]', [y2(:)+y2_sem(:); flipud(y2(:)-y2_sem(:))]',[0.75 0.75 0.75],'FaceAlpha',0.05,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        x_min = min([x1,x2]);
        x_max = max([x1,x2]);
        y_min = min([y1;y2-y2_sem]);
        y_max = max([y1;y2+y2_sem]);
        
        scatter(x1(clusterIndex_memoryPrecision_sigBoolIndex),y_max+(y_max-y_min)*0.08,20,[0 0 0],'*');
        hold on
        
        le = legend(h_line,'Data','Shuffled','Location','southwest','fontsize',8);%10
        le.ItemTokenSize = ones(1,2)*10;
        
        xticks(1:size(dis_range,1));
        
        
        tempStr = string;
        for tempi=1:size(dis_range,1)
            temp1 = [num2str(dis_range(tempi,1)),'-',num2str(dis_range(tempi,2))];
            tempStr = [tempStr; temp1]; %#ok<*AGROW>
            
        end
        tempStr(1) = [];
        set(gca,'xticklabel',tempStr);
        xtickangle(20);
        
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.25 y_max+(y_max-y_min)*0.15]);
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Spatial size (um)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Clustering index', 'FontSize', 12, 'FontWeight', 'bold');
        if if_tuning_location0_precision1 == 1
            temp_title = title(sprintf('Memory precision'),'FontSize',11);
        elseif if_tuning_location0_precision1 == 0
            temp_title = title(sprintf('Memory'),'FontSize',11);
        end
        temp_title.Interpreter = 'none';
        
        
        
        %% meta-memory
        nexttile
        
        h_line = [];
        
        x1 = 1:size(dis_range,1);
        y1 = clusterIndex_choiceMemory_disBin_mean;
        h_line = [h_line plot(x1,y1,'LineWidth',1,'color',[0.25 0.25 0.25])];
        hold on
        
        x2 = 1:size(dis_range,1);
        y2 = clusterIndex_choiceMemory_disBin_shuffled_meanB;
        
        y2_sem = clusterIndex_choiceMemory_disBin_shuffled_semB;
        %y2_sem = clusterIndex_choiceMemory_disBin_shuffled_stdB;
        
        h_line = [h_line plot(x2,y2,'LineWidth',1,'color',[0.75 0.75 0.75])];
        hold on
        
        patch([x2(:); flipud(x2(:))]', [y2(:)+y2_sem(:); flipud(y2(:)-y2_sem(:))]',[0.75 0.75 0.75],'FaceAlpha',0.05,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        x_min = min([x1,x2]);
        x_max = max([x1,x2]);
        y_min = min([y1;y2-y2_sem]);
        y_max = max([y1;y2+y2_sem]);
        
        scatter(x1(clusterIndex_choiceMemory_sigBoolIndex),y_max+(y_max-y_min)*0.08,20,[0 0 0],'*');
        hold on
        
        le = legend(h_line,'Data','Shuffled','Location','southwest','fontsize',8);%10
        le.ItemTokenSize = ones(1,2)*10;
        
        xticks(1:size(dis_range,1));
        
        tempStr = string;
        for tempi=1:size(dis_range,1)
            temp1 = [num2str(dis_range(tempi,1)),'-',num2str(dis_range(tempi,2))];
            tempStr = [tempStr; temp1]; %#ok<*AGROW>
            
        end
        tempStr(1) = [];
        set(gca,'xticklabel',tempStr);
        xtickangle(20);
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.25 y_max+(y_max-y_min)*0.15]);
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Spatial size (um)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Clustering index', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Meta-memory'),'FontSize',11);
        temp_title.Interpreter = 'none';
        
        
        
        %% baseline meta-memory
        nexttile
        
        h_line = [];
        
        x1 = 1:size(dis_range,1);
        y1 = clusterIndex_choiceMemory_baseline_disBin_mean;
        h_line = [h_line plot(x1,y1,'LineWidth',1,'color',[0.25 0.25 0.25])];
        hold on
        
        x2 = 1:size(dis_range,1);
        y2 = clusterIndex_choiceMemory_baseline_disBin_shuffled_meanB;
        
        y2_sem = clusterIndex_choiceMemory_baseline_disBin_shuffled_semB;
        %y2_sem = clusterIndex_choiceMemory_baseline_disBin_shuffled_stdB;
        
        h_line = [h_line plot(x2,y2,'LineWidth',1,'color',[0.75 0.75 0.75])];
        hold on
        
        patch([x2(:); flipud(x2(:))]', [y2(:)+y2_sem(:); flipud(y2(:)-y2_sem(:))]',[0.75 0.75 0.75],'FaceAlpha',0.05,'EdgeColor',[0.62 0.62 0.62]);
        hold on
        
        x_min = min([x1,x2]);
        x_max = max([x1,x2]);
        y_min = min([y1;y2-y2_sem]);
        y_max = max([y1;y2+y2_sem]);
        
        scatter(x1(clusterIndex_choiceMemory_baseline_sigBoolIndex),y_max+(y_max-y_min)*0.08,20,[0 0 0],'*');
        hold on
        
        le = legend(h_line,'Data','Shuffled','Location','southwest','fontsize',8);%10
        le.ItemTokenSize = ones(1,2)*10;
        
        xticks(1:size(dis_range,1));
        
        tempStr = string;
        for tempi=1:size(dis_range,1)
            temp1 = [num2str(dis_range(tempi,1)),'-',num2str(dis_range(tempi,2))];
            tempStr = [tempStr; temp1]; %#ok<*AGROW>
            
        end
        tempStr(1) = [];
        set(gca,'xticklabel',tempStr);
        xtickangle(20);
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.25 y_max+(y_max-y_min)*0.15]);
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Spatial size (um)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Clustering index', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Baseline Meta'),'FontSize',11);
        temp_title.Interpreter = 'none';
        
    end
    
    
    
    
    
    
    
    
end




%% new
if true
    temptempBoolIndex_1 = r2_memoryPrecision_valid>0;
    temptempBoolIndex_2 = r2_choiceMemory_valid>0;
    temptempBoolIndex_3 = r2_choiceMemory_baseline_valid>0;
    
    roi_med;
    temp_roi_med = roi_med(:,2:-1:1);

    
    disNew  = nan(temp_roiNum,temp_roiNum);
    for tempi=1:temp_roiNum
        temp_medA = temp_roi_med(tempi,:);
        for tempj=1:temp_roiNum
            temp_medB = temp_roi_med(tempj,:);
            temp_dis = norm(temp_medA-temp_medB);
            disNew(tempi,tempj) = temp_dis*700/512;
        end
    end
    
    
    
    % Dis within neuron groups
    temp_disNew_1to1 = disNew(temptempBoolIndex_1,temptempBoolIndex_1);
    temp_disNew_2to2 = disNew(temptempBoolIndex_2,temptempBoolIndex_2);
    temp_disNew_3to3 = disNew(temptempBoolIndex_3,temptempBoolIndex_3);
    
    % Dis between neuron groups
    temp_disNew_1to2 = disNew(temptempBoolIndex_1,temptempBoolIndex_2);
    temp_disNew_1to3 = disNew(temptempBoolIndex_1,temptempBoolIndex_3);
    temp_disNew_1toOhters = disNew(temptempBoolIndex_1,temptempBoolIndex_2|temptempBoolIndex_3);
    temp_disNew_2to1 = disNew(temptempBoolIndex_2,temptempBoolIndex_1);
    temp_disNew_2to3 = disNew(temptempBoolIndex_2,temptempBoolIndex_3);
    temp_disNew_2toOhters = disNew(temptempBoolIndex_2,temptempBoolIndex_1|temptempBoolIndex_3);
    temp_disNew_3to1 = disNew(temptempBoolIndex_3,temptempBoolIndex_1);
    temp_disNew_3to2 = disNew(temptempBoolIndex_3,temptempBoolIndex_2);
    temp_disNew_3toOhters = disNew(temptempBoolIndex_3,temptempBoolIndex_1|temptempBoolIndex_2);
    
    
    temp_disNew_1to1_mean = mean(temp_disNew_1to1,2);
    temp_disNew_2to2_mean = mean(temp_disNew_2to2,2);
    temp_disNew_3to3_mean = mean(temp_disNew_3to3,2);
    temp_disNew_1to2_mean = mean(temp_disNew_1to2,2);
    temp_disNew_1to3_mean = mean(temp_disNew_1to3,2);
    temp_disNew_2to1_mean = mean(temp_disNew_2to1,2);
    temp_disNew_2to3_mean = mean(temp_disNew_2to3,2);
    temp_disNew_3to1_mean = mean(temp_disNew_3to1,2);
    temp_disNew_3to2_mean = mean(temp_disNew_3to2,2);
    temp_disNew_1toOthers_mean = mean(temp_disNew_1toOhters,2);
    temp_disNew_2toOthers_mean = mean(temp_disNew_2toOhters,2);
    temp_disNew_3toOthers_mean = mean(temp_disNew_3toOhters,2);
    
    temp_disNew_1to1_mean2 = mean(temp_disNew_1to1_mean);
    temp_disNew_2to2_mean2 = mean(temp_disNew_2to2_mean);
    temp_disNew_3to3_mean2 = mean(temp_disNew_3to3_mean);
    temp_disNew_1to2_mean2 = mean(temp_disNew_1to2_mean);
    temp_disNew_1to3_mean2 = mean(temp_disNew_1to3_mean);
    temp_disNew_2to1_mean2 = mean(temp_disNew_2to1_mean);
    temp_disNew_2to3_mean2 = mean(temp_disNew_2to3_mean);
    temp_disNew_3to1_mean2 = mean(temp_disNew_3to1_mean);
    temp_disNew_3to2_mean2 = mean(temp_disNew_3to2_mean);
    temp_disNew_1toOthers_mean2 = mean(temp_disNew_1toOthers_mean);
    temp_disNew_2toOthers_mean2 = mean(temp_disNew_2toOthers_mean);
    temp_disNew_3toOthers_mean2 = mean(temp_disNew_3toOthers_mean);
    
    [~,temp_p_11_12] = ttest(temp_disNew_1to1_mean,temp_disNew_1to2_mean);
    [~,temp_p_11_13] = ttest(temp_disNew_1to1_mean,temp_disNew_1to3_mean);
    
    [~,temp_p_22_21] = ttest(temp_disNew_2to2_mean,temp_disNew_2to1_mean);
    [~,temp_p_22_23] = ttest(temp_disNew_2to2_mean,temp_disNew_2to3_mean);
    
    [~,temp_p_33_31] = ttest(temp_disNew_3to3_mean,temp_disNew_3to1_mean);
    [~,temp_p_33_32] = ttest(temp_disNew_3to3_mean,temp_disNew_3to2_mean);
        
    if if_plot_distance_pairwise0_centroid1 == 0
        %[~,temp_p_11_1Others] = ttest(temp_disNew_1to1_mean,temp_disNew_1toOthers_mean);
        %[~,temp_p_22_2Others] = ttest(temp_disNew_2to2_mean,temp_disNew_2toOthers_mean);
        %[~,temp_p_33_3Others] = ttest(temp_disNew_3to3_mean,temp_disNew_3toOthers_mean);
        [~,temp_p_11_1Others] = ttest2(temp_disNew_1to1_mean,temp_disNew_1toOthers_mean);
        [~,temp_p_22_2Others] = ttest2(temp_disNew_2to2_mean,temp_disNew_2toOthers_mean);
        [~,temp_p_33_3Others] = ttest2(temp_disNew_3to3_mean,temp_disNew_3toOthers_mean);
        
        fprintf('dis_1to1 = %.1f, dis_1toOthers = %.1f, p = %.3f.\n',...
            temp_disNew_1to1_mean2,temp_disNew_1toOthers_mean2,temp_p_11_1Others);
    
        fprintf('dis_2to2 = %.1f, dis_2toOthers = %.1f, p = %.3f.\n',...
            temp_disNew_2to2_mean2,temp_disNew_2toOthers_mean2,temp_p_22_2Others);
    
        fprintf('dis_3to3 = %.1f, dis_3toOthers = %.1f, p = %.3f.\n',...
            temp_disNew_3to3_mean2,temp_disNew_3toOthers_mean2,temp_p_33_3Others);        
    end
    
    
    
    
    centriod_memoryPrecision_all;
    centriod_meta_all;
    centriod_meta_baseline_all;
    
    
    centriodDis_1toOthers = norm(centriod_memoryPrecision_all-(centriod_meta_all+centriod_meta_baseline_all)/2)*700/512;
    centriodDis_2toOthers = norm(centriod_meta_all-(centriod_memoryPrecision_all+centriod_meta_baseline_all)/2)*700/512;
    centriodDis_3toOthers = norm(centriod_meta_baseline_all-(centriod_memoryPrecision_all+centriod_meta_all)/2)*700/512;
    
    if if_plot_distance_pairwise0_centroid1 == 1
        [~,temp_p_11_1Others] = ttest(temp_disNew_1to1_mean,centriodDis_1toOthers);
        [~,temp_p_22_2Others] = ttest(temp_disNew_2to2_mean,centriodDis_2toOthers);
        [~,temp_p_33_3Others] = ttest(temp_disNew_3to3_mean,centriodDis_3toOthers);
        
        fprintf('dis_1to1 = %.1f, centriodDis_1toOthers = %.1f, p = %.3f.\n',...
            temp_disNew_1to1_mean2,centriodDis_1toOthers,temp_p_11_1Others);
        
        fprintf('dis_2to2 = %.1f, centriodDis_2toOthers = %.1f, p = %.3f.\n',...
            temp_disNew_2to2_mean2,centriodDis_2toOthers,temp_p_22_2Others);
        
        fprintf('dis_3to3 = %.1f, centriodDis_3toOthers = %.1f, p = %.3f.\n',...
            temp_disNew_3to3_mean2,centriodDis_3toOthers,temp_p_33_3Others);
    end
    
    
end


%% new
if true
    temptempBoolIndex_1 = r2_memoryPrecision_valid>0;
    temptempBoolIndex_2 = r2_choiceMemory_valid>0;
    temptempBoolIndex_3 = r2_choiceMemory_baseline_valid>0;
    
    roi_med;
    temp_roi_med = roi_med(:,2:-1:1);
    

    centriod_memoryPrecision_all;
    centriod_meta_all;
    centriod_meta_baseline_all;
    
    disNew2_neuronToCentriod  = nan(temp_roiNum,3);
    for tempi=1:temp_roiNum
        temp_medA = temp_roi_med(tempi,:);
        for tempj=1:3
            %temp_medB = temp_roi_med(tempj,:);
            
            if tempj == 1
                temp_medB = centriod_memoryPrecision_all;
            elseif tempj == 2
                temp_medB = centriod_meta_all;
            elseif tempj == 3
                temp_medB = centriod_meta_baseline_all;
            end
            temp_dis = norm(temp_medA-temp_medB);
            disNew2_neuronToCentriod(tempi,tempj) = temp_dis*700/512;
        end
    end
    
    
    % Dis within neuron groups
    temp_disNew2_1to1 = disNew2_neuronToCentriod(temptempBoolIndex_1,1);
    temp_disNew2_2to2 = disNew2_neuronToCentriod(temptempBoolIndex_2,2);
    temp_disNew2_3to3 = disNew2_neuronToCentriod(temptempBoolIndex_3,3);

    % Dis between neuron groups
    temp_disNew2_1to2 = disNew2_neuronToCentriod(temptempBoolIndex_1,2);
    temp_disNew2_1to3 = disNew2_neuronToCentriod(temptempBoolIndex_1,3);
    temp_disNew2_1toOthers = (temp_disNew2_1to2+temp_disNew2_1to3)/2;
    
    temp_disNew2_2to1 = disNew2_neuronToCentriod(temptempBoolIndex_2,1);
    temp_disNew2_2to3 = disNew2_neuronToCentriod(temptempBoolIndex_2,3);
    temp_disNew2_2toOthers = (temp_disNew2_2to1+temp_disNew2_2to3)/2;
    
    temp_disNew2_3to1 = disNew2_neuronToCentriod(temptempBoolIndex_3,1);
    temp_disNew2_3to2 = disNew2_neuronToCentriod(temptempBoolIndex_3,2);
    temp_disNew2_3toOthers = (temp_disNew2_3to1+temp_disNew2_3to2)/2;
    
    
    [~,temp_p_disNew2_11_1Others] = ttest(temp_disNew2_1to1,temp_disNew2_1toOthers);
    [~,temp_p_disNew2_22_2Others] = ttest(temp_disNew2_2to2,temp_disNew2_2toOthers);
    [~,temp_p_disNew2_33_3Others] = ttest(temp_disNew2_3to3,temp_disNew2_3toOthers);
    
%     fprintf('disNew2_1to1 = %.1f, disNew2_1toOthers = %.1f, p = %.3f.\n',...
%         mean(temp_disNew2_1to1),mean(temp_disNew2_1toOthers),temp_p_disNew2_11_1Others);
%     
%     fprintf('disNew2_2to2 = %.1f, disNew2_2toOthers = %.1f, p = %.3f.\n',...
%         mean(temp_disNew2_2to2),mean(temp_disNew2_2toOthers),temp_p_disNew2_22_2Others);
%     
%     fprintf('disNew2_3to3 = %.1f, disNew2_3toOthers = %.1f, p = %.3f.\n',...
%         mean(temp_disNew2_3to3),mean(temp_disNew2_3toOthers),temp_p_disNew2_33_3Others);

    
end


if false
    temptempBoolIndex_1 = r2_memoryPrecision_valid>0;
    temptempBoolIndex_2 = r2_choiceMemory_valid>0;
    temptempBoolIndex_3 = r2_choiceMemory_baseline_valid>0;
    
    roi_med;
    temp_roi_med = roi_med(:,2:-1:1);
    %temp_roi_med_um = temp_roi_med*700/512;
    temp_roi_med_um = temp_roi_med;
    
    temp_xPos_1 = temp_roi_med_um(temptempBoolIndex_1,1);
    temp_xPos_2 = temp_roi_med_um(temptempBoolIndex_2,1);
    temp_xPos_3 = temp_roi_med_um(temptempBoolIndex_3,1);
    
    [~,temp_p_xPos_12] = ttest2(temp_xPos_1,temp_xPos_2);
    
end

%% Plot Within groups VS. Between groups' distance distribution in one panel
if true
    fig = figure('Name','Centroid distance','NumberTitle','off');
    %set(gcf,'Position',[50+0 400+0 240*0.80 336*1.11*0.67*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50+0 50+0 240*0.80*0.9 336*1.11*0.67*0.9*1.05*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_1 = temp_disNew_1to1_mean';
    temp_2 = temp_disNew_2to2_mean';
    temp_3 = temp_disNew_3to3_mean';
    
    if if_plot_distance_pairwise0_centroid1 == 0
        temp_x1_mean = median(temp_disNew_1toOthers_mean);
        temp_x2_mean = median(temp_disNew_2toOthers_mean);
        temp_x3_mean = median(temp_disNew_3toOthers_mean);
    elseif if_plot_distance_pairwise0_centroid1 == 1
        temp_x1_mean = centriodDis_1toOthers;
        temp_x2_mean = centriodDis_2toOthers;
        temp_x3_mean = centriodDis_3toOthers;
    end
    
    
    temp_p_1 = temp_p_11_1Others;
    temp_p_2 = temp_p_22_2Others;
    temp_p_3 = temp_p_33_3Others;
    
    temp_y_min = min([temp_1 temp_2 temp_3 temp_x1_mean temp_x2_mean temp_x3_mean]);
    temp_y_max = max([temp_1 temp_2 temp_3 temp_x1_mean temp_x2_mean temp_x3_mean]);
    
    temp_data = [temp_3';temp_1';temp_2'];
    
    g1 = repmat({'A'},length(temp_3),1);
    g2 = repmat({'B'},length(temp_1),1);
    g3 = repmat({'C'},length(temp_2),1);
    
    temp_label = [g1;g2;g3];
    
    temptemp_color1 = [1 1 1]*0.5;
    temptemp_color2 = repmat(temptemp_color1, 3, 1);
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',[{'A'};{'B'};{'C'}]);
    h(1).ViolinPlot.FaceAlpha = 0.1;
    h(2).ViolinPlot.FaceAlpha = 0.1;
    h(3).ViolinPlot.FaceAlpha = 0.1;
    
    
    plot([1-0.45 1+0.45],temp_x3_mean*[1 1],'color', [0.25 0.25 0.25],'linewidth',4);%[0.6350 0.0780 0.1840]
    hold on
    
    plot([2-0.45 2+0.45],temp_x1_mean*[1 1],'color', [0.25 0.25 0.25],'linewidth',4);
    hold on
    
    plot([3-0.45 3+0.45],temp_x2_mean*[1 1],'color', [0.25 0.25 0.25],'linewidth',4);
    hold on
    
    tempTxt = sprintf('');
    if temp_p_3 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p_3 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p_3 < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_min-(temp_y_max-temp_y_min)*0.15,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    tempTxt = sprintf('');
    if temp_p_1 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p_1 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p_1 < 0.05
        tempTxt = sprintf('*');
    end
    text(2,temp_y_min-(temp_y_max-temp_y_min)*0.15,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');

    
    tempTxt = sprintf('');
    if temp_p_2 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p_2 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p_2 < 0.05
        tempTxt = sprintf('*');
    end
    text(3,temp_y_min-(temp_y_max-temp_y_min)*0.15,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0.3 3.7]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
    set(gca, 'FontSize', 8) %14
    set(gca,'box','off');% 取消右、上边框
    
    
    %xtl = ["WM", "Meta", "Baseline"];
    xtl = ["BSL", "WM", "Meta"];    
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.38;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',8);%9-->12
    
    set(gca,'xticklabel','');
    
    %yticks([0 40 80]);
    
    ylabel('Distance (um)', 'FontSize', 9, 'FontWeight', 'normal');
    %temp_title = title(sprintf('Centroid distance'),'fontsize',9);
    %temp_title.Interpreter = 'none';
    
end



%% End