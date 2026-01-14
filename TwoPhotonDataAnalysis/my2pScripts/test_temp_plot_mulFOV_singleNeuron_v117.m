% Chuan's 3rd script (20251214)
% This script: To plot multi-FOVs' results, related to many figures in the paper.
%% Initialization
% clear
close all


if_monkey_D0_Z1 = 1;

if_plotTwoMonkey_subspace = 0;%run Ding first, 0,1

if if_monkey_D0_Z1 == 0
    if_plotTwoMonkey_subspace = 0;
end

if if_plotTwoMonkey_subspace == 1
    temp_flag = 1;
    if exist('if_monkey_D0_Z1_old','var') == 1
        if if_monkey_D0_Z1_old == 1
            temp_flag = 0;
        end
    end
    
    if temp_flag == 1
        temp_degree_ab_priorMemory_summary_Ding = temp_degree_ab_priorMemory_summary;
        temp_degree_ac_priorMeta_summary_Ding = temp_degree_ac_priorMeta_summary;
        temp_degree_bc_memoryMeta_summary_Ding = temp_degree_bc_memoryMeta_summary;
        
        temp_degree_aa_prior_summary_Ding = temp_degree_aa_prior_summary;
        temp_degree_bb_memory_summary_Ding = temp_degree_bb_memory_summary;
        temp_degree_cc_meta_summary_Ding = temp_degree_cc_meta_summary;
        
        
        temp_VAF_ratio_ab_priorMemory_summary_Ding = temp_VAF_ratio_ab_priorMemory_summary;
        temp_VAF_ratio_ac_priorMeta_summary_Ding = temp_VAF_ratio_ac_priorMeta_summary;
        temp_VAF_ratio_bc_memoryMeta_summary_Ding = temp_VAF_ratio_bc_memoryMeta_summary;
        
        temp_VAF_ratio_aa_prior_summary_Ding = temp_VAF_ratio_aa_prior_summary;
        temp_VAF_ratio_bb_memory_summary_Ding = temp_VAF_ratio_bb_memory_summary;
        temp_VAF_ratio_cc_meta_summary_Ding = temp_VAF_ratio_cc_meta_summary;
        
        fprintf('Get subspace info from Ding.\n');
    end
end

if_monkey_D0_Z1_old = if_monkey_D0_Z1;

a = 1;

if if_monkey_D0_Z1 == 0
    %fileName = 'AllPart_summary_12FOV_Ding_20240409A_n11nPrecision.mat';
    %fileName = 'AllPart_summary_12FOV_Ding_20240419B_n11nPrecision.mat';
    %fileName = 'AllPart_summary_12FOV_Ding_20240424A_n11nPrecision.mat';
    %fileName = 'AllPart_summary_12FOV_Ding_20240425A_n11nPrecision.mat';
    %fileName = 'AllPart_summary_12FOV_Ding_20240513A_n11nPrecision.mat';
    %fileName = 'AllPart_summary_12FOV_Ding_20240621A_trialBehavior_historyWeight.mat';% trial tuning (behavior)
    %fileName = 'AllPart_summary_12FOV_Ding_20240625A_resampledMeta_historyWeight.mat';% wrong in weighted reward
    %fileName = 'AllPart_summary_12FOV_Ding_20240625B_resampledMeta_historyWeight.mat';
    %fileName = 'AllPart_summary_12FOV_Ding_20240625C_resampledMeta_historyWeight.mat';%temp_cellIndexMapping_AB_mulFOV
    %fileName = 'AllPart_summary_12FOV_Ding_20240625D_resampledMeta_historyWeight.mat';%temp_cellIndexMapping_AB_mulFOV
    %fileName = 'AllPart_summary_12FOV_Ding_20240708A_locationTuning.mat';
    %fileName = 'AllPart_summary_12FOV_Ding_20240710A_locationTuning_newBeta.mat';
    %fileName = 'AllPart_summary_12FOV_Ding_20240712A_locationTuning_newBeta_expHistory.mat';
    %fileName = 'AllPart_summary_12FOV_Ding_20240712B_locationTuning_newBeta_expHistory.mat';% with r_mean_correct_theoretical
    %fileName = 'AllPart_summary_12FOV_Ding_20240712C_locationTuning_newBeta_expHistory.mat';% with correlation between meta-memory and memory precision within each sequence
    %fileName = 'AllPart_summary_12FOV_Ding_20241103A.mat';% relate to metaMonitoringFigure_20241103A
    %fileName = 'AllPart_summary_12FOV_Ding_20241110A.mat';% entropy-based, relate to metaMonitoringFigure_20241111A
    %fileName = 'AllPart_summary_12FOV_Ding_20241111B.mat';% (new) entropy-based, relate to metaMonitoringFigure_20241111A
    %fileName = 'AllPart_summary_12FOV_Ding_20241120A_newFig4.mat';% (new) entropy-based, new fig4
    %fileName = 'AllPart_summary_12FOV_Ding_20241121A_newFig4.mat';% (new) entropy-based, new fig4, new fig5. (Good latency)
    %fileName = 'AllPart_summary_12FOV_Ding_20241201A_newAUROC.mat';% (new all above), new AUROC. (Bad latency)
    %fileName = 'AllPart_summary_12FOV_Ding_20241203A_newAUROC2.mat';% (new all above), new AUROC fot trial history. (choice in meta regress, Good latency)
    %fileName = 'AllPart_summary_12FOV_Ding_20250102A_new.mat';% (new all above), new Fig.4 proj shift, new Fig.5 spatial distance
    %fileName = 'AllPart_summary_12FOV_Ding_20250405A_new.mat';% (new all above), new linear axis angles
    %fileName = 'AllPart_summary_12FOV_Ding_20250406A_new.mat';% (new all above), new unit view chance
    %fileName = 'AllPart_summary_12FOV_Ding_20250407B_VAF.mat';% (new all above), Angle & VAF (if_resample0_twoSession1=1)
    %fileName = 'AllPart_summary_12FOV_Ding_20250408A_lasso.mat';% (new all above), lasso, Angle & VAF (if_resample0_twoSession1=0)
    %fileName = 'AllPart_summary_12FOV_Ding_20250408B_EV.mat';% (new all above), EV, pca
    %fileName = 'AllPart_summary_12FOV_Ding_20250416A_subspaceSpatial.mat';% (new all above), subspace-based spatial
    %fileName = 'AllPart_summary_12FOV_Ding_20250420A_subSpatial_subR2.mat';% (new all above), subspace-based spatial, subspace r2
    fileName = 'AllPart_summary_12FOV_Ding_20250826A_unitSpatial_review.mat';% (new all above), unit-based spatial, review comments update
    
elseif if_monkey_D0_Z1 == 1
    %fileName = 'AllPart_summary_16FOV_Zelku_20240415A_n11nPrecision.mat';
    %fileName = 'AllPart_summary_16FOV_Zelku_20240420A_n11nPrecision.mat';
    %fileName = 'AllPart_summary_16FOV_Zelku_20240424A_n11nPrecision.mat';
    %fileName = 'AllPart_summary_16FOV_Zelku_20240426A_n11nPrecision_mixAndtTrial.mat';% mix for beta, trial for sptial
    %fileName = 'AllPart_summary_16FOV_Zelku_20240514A_n11nPrecision.mat';
    %fileName = 'AllPart_summary_16FOV_Zelku_20240621A_trialBehavior_historyWeight.mat';% trial tuning (behavior)
    %fileName = 'AllPart_summary_16FOV_Zelku_20240625A_resampledMeta_historyWeight.mat';% wrong in weighted reward
    %fileName = 'AllPart_summary_16FOV_Zelku_20240625B_resampledMeta_historyWeight.mat';
    %fileName = 'AllPart_summary_16FOV_Zelku_20240625C_resampledMeta_historyWeight.mat';%temp_cellIndexMapping_AB_mulFOV
    %fileName = 'AllPart_summary_16FOV_Zelku_20240625D_resampledMeta_historyWeight.mat';%temp_cellIndexMapping_AB_mulFOV
    %fileName = 'AllPart_summary_16FOV_Zelku_20240708A_locationTuning.mat';
    %fileName = 'AllPart_summary_16FOV_Zelku_20240710A_locationTuning_newBeta.mat';
    %fileName = 'AllPart_summary_16FOV_Zelku_20240712A_locationTuning_newBeta_expHistory.mat';
    %fileName = 'AllPart_summary_16FOV_Zelku_20240712B_locationTuning_newBeta_expHistory.mat';% with r_mean_correct_theoretical
    %fileName = 'AllPart_summary_16FOV_Zelku_20240712C_locationTuning_newBeta_expHistory.mat';% with correlation between meta-memory and memory precision within each sequence
    %fileName = 'AllPart_summary_16FOV_Zelku_20241103A.mat';% relate to metaMonitoringFigure_20241103A
    %fileName = 'AllPart_summary_16FOV_Zelku_20241110A.mat';% entropy-based, relate to metaMonitoringFigure_20241111A
    %fileName = 'AllPart_summary_16FOV_Zelku_20241111B.mat';% (new) entropy-based, relate to metaMonitoringFigure_20241111A
    %fileName = 'AllPart_summary_16FOV_Zelku_20241120A_newFig4.mat';% (new) entropy-based, new fig4
    %fileName = 'AllPart_summary_16FOV_Zelku_20241121A_newFig4.mat';% (new) entropy-based, new fig4, new fig5. (Good latency)
    %fileName = 'AllPart_summary_16FOV_Zelku_20241201A_newAUROC.mat';% (new all above), new AUROC. (Bad latency)
    %fileName = 'AllPart_summary_16FOV_Zelku_20241203A_newAUROC2.mat';% (new all above), new AUROC fot trial history. (choice in meta regress, Good latency)
    %fileName = 'AllPart_summary_16FOV_Zelku_20250102A_new.mat';% (new all above), new Fig.4 proj shift, new Fig.5 spatial distance
    %fileName = 'AllPart_summary_16FOV_Zelku_20250405A_new.mat';% (new all above), new linear axis angles
    %fileName = 'AllPart_summary_16FOV_Zelku_20250406A_new.mat';% (new all above), new unit view chance
    %fileName = 'AllPart_summary_16FOV_Zelku_20250407B_VAF.mat';% (new all above), Angle & VAF (if_resample0_twoSession1=1)
    %fileName = 'AllPart_summary_16FOV_Zelku_20250408B_EV.mat';% (new all above), EV, pca
    %fileName = 'AllPart_summary_16FOV_Zelku_20250416A_subspaceSpatial.mat';% (new all above), subspace-based spatial
    %fileName = 'AllPart_summary_16FOV_Zelku_20250420A_subSpatial_subR2.mat';% (new all above), subspace-based spatial, subspace r2
    fileName = 'AllPart_summary_16FOV_Zelku_20250826A_unitSpatial_review.mat';% (new all above), unit-based spatial, review comments update
    
end


if_compute_part6 = 0;
if_ePRIRS = 0;

if_plot = 1;

if_plot_part2 = 0;
if_plot_part3 = 0;
if_plot_part4 = 0;
if_plot_part5 = 1;
if_plot_part6 = 0;
if_plot_cellIndexMapping = 0;

if_cellIndexMapping = 1;

if_entropy_behaviorInverted = 0;


if_plot_distance_pairwise0_centroid1 = 0;%1

if_PCA_feature0_neuron1 = 1;

if_PCA_singleFOV0_multiFOV1 = 1;%1

if_locTuningFit_lasso0_fitglm1 = 1;

if_tuning_location0_precision1 = 0;

if_tuning_r20_beta1 = 1;%0


if_useSelectiveNeuron = 0;

if_plot_3d_beta = 0;

if_plotBeta_trialTuning0_seqTuning1_mix2 = 2;%1,2

if_plot_violinplot0_pairline1 = 1;

if_plot_additionalSmooth = 1;


% allPart_summary_mulFOV;

targetPATHA = 'C:\ASDROOT\STUDY\temp\data\AllPart_summary_mulFOV'; %#ok<*UNRCH>
cd(targetPATHA)


fprintf('Load %s.\n',fileName);
load(fileName);

allPart_summary_mulFOV = allPart_summary_mulFOV; %#ok<*ASGSL>

% if false
targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)
% end

color_choiceMemoryLow = [1,133,113]/255; %[5,113,176]/255
color_choiceMemoryHigh = [128,205,193]/255; %[146,197,222]/255
color_choiceOffloadLow = [223,194,125]/255; %[244,165,130]/255
color_choiceOffloadHigh = [166,97,26]/255; %[202,0,32]/255

color_memoryQuality = [252,141,89]/255;
color_meta = [145,191,219]/255;
color_prior = [217,217,162]/255;


r2_seqPrecision_summary = allPart_summary_mulFOV.part3_summary_mulFOV.r2_seqPrecision_summary;
r2_rProb_summary = allPart_summary_mulFOV.part3_summary_mulFOV.r2_rProb_summary;

beta_memoryPrecision_summary = allPart_summary_mulFOV.part6_summary_mulFOV.beta_memoryPrecision_summary;
beta_choiceMemory_summary = allPart_summary_mulFOV.part6_summary_mulFOV.beta_choiceMemory_summary;
beta_choiceMemory_baseline_summary = allPart_summary_mulFOV.part6_summary_mulFOV.beta_choiceMemory_baseline_summary;
% beta_choiceMemory_baseline_summary = std_beta_6loc;

beta_seqPrecision_summary = allPart_summary_mulFOV.part6_summary_mulFOV.beta_seqPrecision_summary;
beta_gProb_summary = allPart_summary_mulFOV.part6_summary_mulFOV.beta_gProb_summary;
beta_gProb_baseline_summary = allPart_summary_mulFOV.part6_summary_mulFOV.beta_gProb_baseline_summary;

% beta_seqPrecision_summary = std_beta_6loc;
% beta_gProb_summary = std_beta_6loc;

% r2_memoryPrecision_summary = allPart_summary_mulFOV.part6_summary_mulFOV.r2_memoryPrecision_summary;
% r2_choiceMemory_summary = allPart_summary_mulFOV.part6_summary_mulFOV.r2_choiceMemory_summary;
r2_choiceMemory_baseline_summary = allPart_summary_mulFOV.part6_summary_mulFOV.r2_choiceMemory_baseline_summary;
% % r2_seqPrecision_summary = allPart_summary_mulFOV.part6_summary_mulFOV.r2_seqPrecision_summary;
r2_gProb_summary = allPart_summary_mulFOV.part6_summary_mulFOV.r2_gProb_summary;
% r2_gProb_baseline_summary = allPart_summary_mulFOV.part6_summary_mulFOV.r2_gProb_baseline_summary;
%
%
r2_6loc_summary = allPart_summary_mulFOV.part6_summary_mulFOV.r2_6loc_summary;


beta_6loc_summary = allPart_summary_mulFOV.part6_summary_mulFOV.beta_6loc_summary;
%
% if false
%     beta_6loc_summary = beta_6loc_raw;
% end
%

if true
    temp1 = max(beta_6loc_summary,[],2);
    temp2 = min(beta_6loc_summary,[],2);
    temp3 = nan(length(beta_6loc_summary),1);
    for tempi=1:length(beta_6loc_summary)
        if abs(temp1(tempi)) > abs(temp2(tempi))
            temp3(tempi) = temp1(tempi);
        else
            temp3(tempi) = temp2(tempi);
        end
    end
    beta_6loc_summary_peak = temp3;
end


if false
    % if true
    beta_6loc_summary_peak = std(beta_6loc_summary,0,2);
end

if true
    temp1 = beta_6loc_summary;
    temp2 = norm(temp1);
    
    beta_6loc_summary_norm = nan(size(beta_6loc_summary,1),1);
    
    for tempi=1:size(beta_6loc_summary,1)
        temp1 = beta_6loc_summary(tempi,:);
        beta_6loc_summary_norm(tempi) = norm(temp1);
    end
    
    beta_gProb_summary_norm = abs(beta_gProb_summary);
    beta_choiceMemory_baseline_summary_norm = abs(beta_choiceMemory_baseline_summary);
end


% r2_seqPrecision_summary_n11n = r2_seqPrecision_summary;
% r2_seqPrecision_summary_n11n(r2_seqPrecision_summary_n11n<0) = 0;
% % r2_seqPrecision_summary_n11n = rescale(r2_seqPrecision_summary_n11n,0,1);
% r2_seqPrecision_summary_n11n = rescale(r2_seqPrecision_summary_n11n,-1,1);
%
% r2_gProb_summary_n11n = r2_gProb_summary;
% r2_gProb_summary_n11n(r2_gProb_summary_n11n<0) = 0;
% % r2_gProb_summary_n11n = rescale(r2_gProb_summary_n11n,0,1);
% r2_gProb_summary_n11n = rescale(r2_gProb_summary_n11n,-1,1);
%
% r2_choiceMemory_baseline_summary_n11n = r2_choiceMemory_baseline_summary;
% r2_choiceMemory_baseline_summary_n11n(r2_choiceMemory_baseline_summary_n11n<0) = 0;
% % r2_choiceMemory_baseline_summary_n11n = rescale(r2_choiceMemory_baseline_summary_n11n,0,1);
% r2_choiceMemory_baseline_summary_n11n = rescale(r2_choiceMemory_baseline_summary_n11n,-1,1);
%
% r2_6loc_summary_n11n = r2_6loc_summary;
% r2_6loc_summary_n11n(r2_6loc_summary_n11n<0) = 0;
% % r2_6loc_summary_n11n = rescale(r2_6loc_summary_n11n,0,1);
% r2_6loc_summary_n11n = rescale(r2_6loc_summary_n11n,-1,1);


% beta_seqPrecision_summary = r2_6loc_summary;
% beta_gProb_summary = r2_gProb_summary;
% beta_choiceMemory_baseline_summary = r2_choiceMemory_baseline_summary;

% beta_seqPrecision_summary = r2_6loc_summary_n11n;

% beta_seqPrecision_summary = r2_seqPrecision_summary_n11n;
% beta_gProb_summary = r2_gProb_summary_n11n;
% beta_choiceMemory_baseline_summary = r2_choiceMemory_baseline_summary_n11n;


temp1_raw = beta_choiceMemory_baseline_summary;
temp1 = temp1_raw-mean(temp1_raw);
if abs(max(temp1)) > abs(min(temp1))
    temp1 = temp1./abs(max(temp1));
else
    temp1 = temp1./abs(min(temp1));
end
beta_choiceMemory_baseline_summary = temp1;
% beta_choiceMemory_baseline_summary = temp1_raw.^2;


% beta_seqPrecision_summary = beta_6loc_summary_peak;

temp1_raw = beta_seqPrecision_summary;
temp1 = temp1_raw-mean(temp1_raw);
if abs(max(temp1)) > abs(min(temp1))
    temp1 = temp1./abs(max(temp1));
else
    temp1 = temp1./abs(min(temp1));
end
beta_seqPrecision_summary = temp1;
% beta_seqPrecision_summary = temp1_raw.^2;

temp1_raw = beta_gProb_summary;
temp1 = temp1_raw-mean(temp1_raw);
if abs(max(temp1)) > abs(min(temp1))
    temp1 = temp1./abs(max(temp1));
else
    temp1 = temp1./abs(min(temp1));
end
beta_gProb_summary = temp1;
% beta_gProb_summary = temp1_raw.^2;


% temp1_raw = beta_6loc_summary_peak;
% temp1 = temp1_raw-mean(temp1_raw);
% if abs(max(temp1)) > abs(min(temp1))
%     temp1 = temp1./abs(max(temp1));
% else
%     temp1 = temp1./abs(min(temp1));
% end
% beta_6loc_summary_peak = temp1;
% % beta_6loc_summary_peak = temp1_raw.^2;


%% plot_cellIndexMapping
if if_plot == 1 && if_cellIndexMapping == 1
    temp_cellIndexMapping_AB_mulFOV = allPart_summary_mulFOV.temp_cellIndexMapping_AB_mulFOV;
    
    numFOV = length(temp_cellIndexMapping_AB_mulFOV);
    
    temp_cellIndexMapping_AB_mulFOV_collapsed = [];
    for tempi=1:numFOV
        temp_cellIndexMapping_AB_mulFOV_collapsed = [temp_cellIndexMapping_AB_mulFOV_collapsed; temp_cellIndexMapping_AB_mulFOV{tempi}];
    end
    
    temp_FOVIndex_mulFOV = nan(size(temp_cellIndexMapping_AB_mulFOV_collapsed));
    for tempi=1:numFOV
        temp1 = sum(~isnan(temp_FOVIndex_mulFOV))+1;
        temp2 = temp1-1 + length(temp_cellIndexMapping_AB_mulFOV{tempi});
        temp_range = temp1:temp2;
        
        temp_FOVIndex_mulFOV(temp_range) = tempi;
    end
    
    p_tuning = allPart_summary_mulFOV.p_tuning;
    
    p_6loc_mulFOV_collapsed = [];
    p_memoryPrecision_mulFOV_collapsed = [];
    p_choiceMemory_mulFOV_collapsed = [];
    p_choiceMemory_baseline_mulFOV_collapsed = [];
    
    p_seqPrecision_mulFOV_collapsed = [];
    p_gProb_mulFOV_collapsed = [];
    p_gProb_baseline_mulFOV_collapsed = [];
    
    boolIndex_loc_mulFOV_collapsed = [];
    
    for tempi=1:numFOV
        p_6loc_mulFOV_collapsed = [p_6loc_mulFOV_collapsed; p_tuning.p_6loc_mulFOV{tempi}];
        
        p_memoryPrecision_mulFOV_collapsed = [p_memoryPrecision_mulFOV_collapsed; p_tuning.p_memoryPrecision_mulFOV{tempi}];
        p_choiceMemory_mulFOV_collapsed = [p_choiceMemory_mulFOV_collapsed; p_tuning.p_choiceMemory_mulFOV{tempi}];
        p_choiceMemory_baseline_mulFOV_collapsed = [p_choiceMemory_baseline_mulFOV_collapsed; p_tuning.p_choiceMemory_baseline_mulFOV{tempi}];
        
        p_seqPrecision_mulFOV_collapsed = [p_seqPrecision_mulFOV_collapsed; p_tuning.p_seqPrecision_mulFOV{tempi}];
        p_gProb_mulFOV_collapsed = [p_gProb_mulFOV_collapsed; p_tuning.p_gProb_mulFOV{tempi}];
        p_gProb_baseline_mulFOV_collapsed = [p_gProb_baseline_mulFOV_collapsed; p_tuning.p_gProb_baseline_mulFOV{tempi}];
        
        boolIndex_loc_mulFOV_collapsed = [boolIndex_loc_mulFOV_collapsed; p_tuning.boolIndex_loc_mulFOV{tempi}];
    end
    
    p_memoryPrecision_mulFOV_collapsed;
    p_choiceMemory_mulFOV_collapsed;
    p_choiceMemory_baseline_mulFOV_collapsed;
    
    p_seqPrecision_mulFOV_collapsed;
    p_gProb_mulFOV_collapsed;
    p_gProb_baseline_mulFOV_collapsed;
    
    boolIndex_loc_mulFOV_collapsed;
    
    p_threshold = 0.01;%0.01,0.001
    
    numRoi = length(p_memoryPrecision_mulFOV_collapsed);
    
    temptempBoolIndex_1 = p_memoryPrecision_mulFOV_collapsed<p_threshold;
    temptempBoolIndex_2 = p_choiceMemory_mulFOV_collapsed<p_threshold;
    temptempBoolIndex_3 = p_choiceMemory_baseline_mulFOV_collapsed<p_threshold;
    
    temptempBoolIndex_4 = p_seqPrecision_mulFOV_collapsed<p_threshold;
    temptempBoolIndex_5 = p_gProb_mulFOV_collapsed<p_threshold;
    temptempBoolIndex_6 = p_gProb_baseline_mulFOV_collapsed<p_threshold;
    
    fprintf('if_locTuningFit_lasso0_fitglm1 = %d.\n',if_locTuningFit_lasso0_fitglm1);
    if if_locTuningFit_lasso0_fitglm1 == 0
        temptempBoolIndex_7 = boolIndex_loc_mulFOV_collapsed;
    elseif if_locTuningFit_lasso0_fitglm1 == 1
        temptempBoolIndex_7 = p_6loc_mulFOV_collapsed<p_threshold;
    end
    
    %if true
    if false
        temp1 = selevtivity_multiFOV.selectiveCellBoolIndex_rProb_glm_delay1_multiFOV(temp_cellIndexMapping_AB_mulFOV_collapsed);
        temptempBoolIndex_5 = temp1;
    end
    
    if false
        %% Seq-level glm for location tuning
        
        x = boolIndex_location_seq_T(1:41,:);
        y = dff_seqMean_delay1_multiFOV(temp_cellIndexMapping_AB_mulFOV_collapsed,1:41);
        
        temp_glm_beta = zeros(size(y,1),numFrames);
        temp_glm_pValue = zeros(size(y,1),numFrames);
        %temp_glm_pValue = zeros(size(y,1),1);
        temp_glm_r2 = zeros(size(y,1),1);
        parfor tempi=1:size(y,1)
            %for tempi=1:1
            warning('off');
            temp_x = x;
            temp_y = y(tempi,:)';
            temp_mdl= fitglm(temp_x,temp_y,'linear','Distribution','normal','Intercept',false);
            temp_glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
            temp_glm_pValue(tempi,:) = temp_mdl.Coefficients.pValue;
            %temp_glm_pValue(tempi) = coefTest(temp_mdl);
            temp_glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
            warning('on');
        end
        sum(temp_glm_r2>0.05)
        
    end
    
    
    if if_useSelectiveNeuron == 1
        numRoi_raw = numRoi;
        %selectiveNeuronBoolIndex = false(numRoi,1);
        selectiveNeuronBoolIndex = temptempBoolIndex_3 | temptempBoolIndex_5 | temptempBoolIndex_7;
        numRoi = sum(selectiveNeuronBoolIndex);
        
        temptempBoolIndex_1 = temptempBoolIndex_1(selectiveNeuronBoolIndex);
        temptempBoolIndex_2 = temptempBoolIndex_2(selectiveNeuronBoolIndex);
        temptempBoolIndex_3 = temptempBoolIndex_3(selectiveNeuronBoolIndex);
        
        temptempBoolIndex_4 = temptempBoolIndex_4(selectiveNeuronBoolIndex);
        temptempBoolIndex_5 = temptempBoolIndex_5(selectiveNeuronBoolIndex);
        temptempBoolIndex_6 = temptempBoolIndex_6(selectiveNeuronBoolIndex);
        
        temptempBoolIndex_7 = temptempBoolIndex_7(selectiveNeuronBoolIndex);
        
    elseif if_useSelectiveNeuron == 0
        selectiveNeuronBoolIndex = true(numRoi,1);
    end
    
    temptempBoolIndex_45 = temptempBoolIndex_4 & temptempBoolIndex_5;
    temptempBoolIndex_4_not5 = temptempBoolIndex_4 & ~temptempBoolIndex_5;
    temptempBoolIndex_5_not4 = temptempBoolIndex_5 & ~temptempBoolIndex_4;
    
    temptempBoolIndex_57 = temptempBoolIndex_5 & temptempBoolIndex_7;
    temptempBoolIndex_5_not7 = temptempBoolIndex_5 & ~temptempBoolIndex_7;
    temptempBoolIndex_7_not5 = temptempBoolIndex_7 & ~temptempBoolIndex_5;
    
    temptempBoolIndex_35 = temptempBoolIndex_3 & temptempBoolIndex_5;
    temptempBoolIndex_37 = temptempBoolIndex_3 & temptempBoolIndex_7;
    temptempBoolIndex_357 = temptempBoolIndex_3 & temptempBoolIndex_5 & temptempBoolIndex_7;
    
    %fprintf('memoryPrecision selective %d/%d=%.1f%%.\n',sum(temptempBoolIndex_1),numRoi,(sum(temptempBoolIndex_1)/numRoi)*100);
    %fprintf('gProb_baseline selective %d/%d=%.0f%%.\n',sum(temptempBoolIndex_6),numRoi,(sum(temptempBoolIndex_6)/numRoi)*100);
    
    
    fprintf('Location selective %d/%d=%.1f%%.\n',sum(temptempBoolIndex_7),numRoi,(sum(temptempBoolIndex_7)/numRoi)*100);
    fprintf('choiceMemory selective %d/%d=%.1f%%.\n',sum(temptempBoolIndex_2),numRoi,(sum(temptempBoolIndex_2)/numRoi)*100);
    
    %fprintf('seqPrecision selective %d/%d=%.1f%%.\n',sum(temptempBoolIndex_4),numRoi,(sum(temptempBoolIndex_4)/numRoi)*100);
    %fprintf('gProb selective %d/%d=%.1f%%.\n',sum(temptempBoolIndex_5),numRoi,(sum(temptempBoolIndex_5)/numRoi)*100);
    
    %fprintf('seqPrecision, gProb [mixAB,pureA,pureB] selective [%d,%d,%d] [%.1f%%,%.1f%%,%.1f%%].\n',...
    %    sum(temptempBoolIndex_45),sum(temptempBoolIndex_4_not5),sum(temptempBoolIndex_5_not4),sum(temptempBoolIndex_45)/numRoi*100,sum(temptempBoolIndex_4_not5)/numRoi*100,sum(temptempBoolIndex_5_not4)/numRoi*100);
    
    
    fprintf('gProb selective %d/%d=%.1f%%.\n',sum(temptempBoolIndex_5),numRoi,(sum(temptempBoolIndex_5)/numRoi)*100);
    
    fprintf('Location, gProb [mixAB,pureA,pureB] selective [%d,%d,%d] [%.1f%%,%.1f%%,%.1f%%].\n',...
        sum(temptempBoolIndex_57),sum(temptempBoolIndex_7_not5),sum(temptempBoolIndex_5_not7),sum(temptempBoolIndex_57)/numRoi*100,sum(temptempBoolIndex_7_not5)/numRoi*100,sum(temptempBoolIndex_5_not7)/numRoi*100);
    
    fprintf('choiceMemory_baseline selective %d/%d=%.1f%%.\n',sum(temptempBoolIndex_3),numRoi,(sum(temptempBoolIndex_3)/numRoi)*100);
    
    fprintf('choiceMemory_baseline & Location %d/%d=%.1f%%.\n',sum(temptempBoolIndex_37),numRoi,(sum(temptempBoolIndex_37)/numRoi)*100);
    fprintf('choiceMemory_baseline & gProbselective %d/%d=%.1f%%.\n',sum(temptempBoolIndex_35),numRoi,(sum(temptempBoolIndex_35)/numRoi)*100);
    fprintf('choiceMemory_baseline & gProbselective & Location %d/%d=%.1f%%.\n',sum(temptempBoolIndex_357),numRoi,(sum(temptempBoolIndex_357)/numRoi)*100);
    
    temptempBoolIndex_7;
    temptempBoolIndex_2;
    temp_FOVIndex_mulFOV;
    
    sum_temptempBoolIndex_7_eachFOV = nan(numFOV,1);
    sum_temptempBoolIndex_2_eachFOV = nan(numFOV,1);
    numRoi_eachFOV = nan(numFOV,1);
    for tempi=1:numFOV
        temp1 = temp_FOVIndex_mulFOV==tempi;
        
        sum_temptempBoolIndex_7_eachFOV(tempi) = sum(temptempBoolIndex_7 & temp1);
        sum_temptempBoolIndex_2_eachFOV(tempi) = sum(temptempBoolIndex_2 & temp1);
        numRoi_eachFOV(tempi) = sum(temp1);
    end
    
    
    if if_plot_cellIndexMapping
        
        fig = figure('Name','population tuning std','NumberTitle','off');
        set(gcf,'Position',[10 50 233*0.7*0.9 190*0.85*0.7*0.95*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        temp_1 = sum_temptempBoolIndex_7_eachFOV./numRoi_eachFOV*100';
        temp_2 = sum_temptempBoolIndex_2_eachFOV./numRoi_eachFOV*100';
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        
        temp_data = [temp_1;temp_2];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        
        temp_label = [g1;g2];
        
        temptemp_color1 = [1 1 1]*0.5;%0.5
        temptemp_color2 = repmat(temptemp_color1, 2, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        
        
        
        set(gca,'linewidth',1.5)
        %xlim([-0.5 1.5])
        xlim([0.5 2.5]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);
        set(gca, 'FontSize', 8)
        %set(gca,'XTickLabel', ["Location"; "Choice"],'FontSize', 8);%给坐标加标签
        set(gca,'XTickLabel', ["WM"; "Decision"],'FontSize', 8);%给坐标加标签
        xtickangle(0);
        
        %yticklabels('');
        %yticks([]);
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Selective proportion', 'FontSize', 8, 'FontWeight', 'normal');
        
        if if_monkey_D0_Z1 == 0
            temp_str = 'Monkey D';
        elseif if_monkey_D0_Z1 == 1
            temp_str = 'Monkey Z';
        end
        title(sprintf('%s, n=%d',temp_str,numRoi),'fontsize',9);
        
    end
    
    if if_plot_cellIndexMapping
        
        fig = figure('Name','population tuning std','NumberTitle','off');
        set(gcf,'Position',[210 50 233*0.7*0.9 190*0.85*0.7*0.95*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        temp_1 = (sum(temptempBoolIndex_7)/numRoi)*100;
        temp_2 = (sum(temptempBoolIndex_2)/numRoi)*100;
        
        
        temp_y_min = 0;
        temp_y_max = max([temp_1;temp_2]);
        
        
        bar([0 1], [temp_1 temp_2],0.5, ...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        
        tempTxt = sprintf('%.1f%%',temp_1);
        text(0,temp_1+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color','black','FontSize',8,'FontWeight','normal',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('%.1f%%',temp_2);
        text(1,temp_2+(temp_y_max-temp_y_min)*0.13,tempTxt,'Color','black','FontSize',8,'FontWeight','normal',...
            'HorizontalAlignment','center');
        
        
        
        set(gca,'linewidth',1.5)
        xlim([-0.5 1.5])
        ylim([temp_y_min (temp_y_max-temp_y_min)*1.3]);
        set(gca, 'FontSize', 8)
        %set(gca,'XTickLabel', ["Location"; "Choice"],'FontSize', 8);%给坐标加标签
        set(gca,'XTickLabel', ["WM"; "Decision"],'FontSize', 8);%给坐标加标签
        xtickangle(0);
        
        yticklabels('');
        yticks([]);
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Selective proportion', 'FontSize', 8, 'FontWeight', 'normal');
        
        if if_monkey_D0_Z1 == 0
            temp_str = 'Monkey D';
        elseif if_monkey_D0_Z1 == 1
            temp_str = 'Monkey Z';
        end
        title(sprintf('%s, n=%d',temp_str,numRoi),'fontsize',9);
        
    end
    
end

%% Part 2G: Location correlation of correct, stimuliError and responseError trials
r_mean_correct_summary = allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_summary;
r_mean_stimuliError_summary = allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_stimuliError_summary;
r_mean_responseError_summary = allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_responseError_summary;
r_mean_chanceLevel_summary = allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_chanceLevel_summary;
r_mean_correct_theoretical_summary = allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_theoretical_summary;

temp_1 = r_mean_correct_summary;
temp_2 = r_mean_stimuliError_summary;
temp_3 = r_mean_responseError_summary;
temp_1_chanceLevel = r_mean_chanceLevel_summary;

[~,temp_p12,~,~] = ttest(r_mean_correct_summary,r_mean_stimuliError_summary);
[~,temp_p13,~,~] = ttest(r_mean_correct_summary,r_mean_responseError_summary); %#ok<*ASGLU>
[~,temp_p23,~,~] = ttest(r_mean_stimuliError_summary,r_mean_responseError_summary);


if if_monkey_D0_Z1 == 1
    temptempBoolIndex = temp_1>0.33;
    [~,temp_p12,~,~] = ttest(temp_1(temptempBoolIndex),temp_2(temptempBoolIndex));
    [~,temp_p13,~,~] = ttest(temp_1(temptempBoolIndex),temp_3(temptempBoolIndex));
    [~,temp_p23,~,~] = ttest(temp_3(temptempBoolIndex),mean(temp_2(temptempBoolIndex)),'Tail','right');
end

if false
    
    temptempBoolIndex = temp_1>0.33;
    [~,temp_p12,~,~] = ttest(temp_1(temptempBoolIndex),temp_2(temptempBoolIndex));
    [~,temp_p23,~,~] = ttest(temp_3(temptempBoolIndex),temp_2(temptempBoolIndex));
    [~,temp_p12,~,~] = ttest(temp_1(temptempBoolIndex),temp_2(temptempBoolIndex),'Tail','right');
    [~,temp_p23,~,~] = ttest(temp_3(temptempBoolIndex),temp_2(temptempBoolIndex),'Tail','right');
    [~,temp_p23,~,~] = ttest(temp_3(temptempBoolIndex),mean(temp_2(temptempBoolIndex)),'Tail','right');
    
    
    [~,temp_p12,~,~] = ttest(r_mean_correct_summary,r_mean_stimuliError_summary,'Tail','right');
    [~,temp_p13,~,~] = ttest(r_mean_correct_summary,r_mean_responseError_summary,'Tail','right');
    [~,temp_p23,~,~] = ttest(r_mean_responseError_summary,r_mean_stimuliError_summary,'Tail','right');
    
    [~,temp_p12,~,~] = ttest2(r_mean_correct_summary,r_mean_stimuliError_summary);
    [~,temp_p13,~,~] = ttest2(r_mean_correct_summary,r_mean_responseError_summary);
    [~,temp_p23,~,~] = ttest2(r_mean_stimuliError_summary,r_mean_responseError_summary);
    
    [~,temp_p12,~,~] = ttest2(r_mean_correct_summary,r_mean_stimuliError_summary,'Tail','right');
    [~,temp_p13,~,~] = ttest2(r_mean_correct_summary,r_mean_responseError_summary,'Tail','right');
    [~,temp_p23,~,~] = ttest2(r_mean_responseError_summary,r_mean_stimuliError_summary,'Tail','right');
    
    [~,temp_p13_2,~,~] = ttest2([r_mean_correct_summary r_mean_responseError_summary],r_mean_stimuliError_summary);
    [~,temp_p13_2,~,~] = ttest2([r_mean_correct_summary r_mean_responseError_summary],r_mean_stimuliError_summary,'Tail','right');
    
    
    x = [1*ones(size(temp_1,2),1);2*ones(size(temp_3,2),1);3*ones(size(temp_2,2),1)];
    y = [temp_1,temp_3,temp_2]';
    
    temp_mdl = fitglm(x,y,'linear');
    p_linear = temp_mdl.Coefficients.pValue(2);
    
    
    temp_p12 = ranksum(r_mean_correct_summary,r_mean_stimuliError_summary,'Tail','right');
    temp_p13 = ranksum(r_mean_correct_summary,r_mean_responseError_summary,'Tail','right');
    temp_p23 = ranksum(r_mean_responseError_summary,r_mean_stimuliError_summary,'Tail','right');
    
end


if if_plot == 1 && if_plot_part2 == 1
    fig = figure('Name','locDistri','NumberTitle','off');
    set(gcf,'Position',[50+0 50+0 240*0.80 336*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    temp3_SEM = std(temp_3)/sqrt(length(temp_3));
    
    temp_y_min = mean(temp_1_chanceLevel);
    %temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
    temp_y_max = max([temp_1 temp_2 temp_3]);
    
    temp_y_max12 = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM]);
    temp_y_max13 = max([mean(temp_1)+temp1_SEM,mean(temp_3)+temp3_SEM]);
    temp_y_max23 = max([mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM]);
    
    temp_data = [temp_1';temp_2';temp_3'];
    
    g1 = repmat({'A'},length(temp_1),1);
    g2 = repmat({'B'},length(temp_2),1);
    g3 = repmat({'C'},length(temp_3),1);
    
    temp_label = [g1;g2;g3];
    
    temptemp_color1 = [1 1 1]*0.5;%0.5
    temptemp_color2 = repmat(temptemp_color1, 3, 1);
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',[{'A'};{'B'};{'C'}]);
    h(1).ViolinPlot.FaceAlpha = 0.1;
    h(2).ViolinPlot.FaceAlpha = 0.1;
    h(3).ViolinPlot.FaceAlpha = 0.1;
    
    plot([0 4],mean(temp_1_chanceLevel)*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    tempTxt = sprintf('');
    if temp_p12 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p12 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p12 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.06,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    %     tempTxt = sprintf('');
    %     if temp_p13 < 0.001
    %         tempTxt = sprintf('***');
    %     elseif temp_p13 < 0.01
    %         tempTxt = sprintf('**');
    %     elseif temp_p13 < 0.05
    %         tempTxt = sprintf('*');
    %     end
    %     text(2,temp_y_max+(temp_y_max-temp_y_min)*0.15,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
    %         'HorizontalAlignment','center');
    %     plot([1.1 2.9],temp_y_max+(temp_y_max-temp_y_min)*0.13*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    %     hold on
    
    tempTxt = sprintf('');
    if temp_p23 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p23 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p23 < 0.05
        tempTxt = sprintf('*');
    end
    text(2.5,temp_y_max+(temp_y_max-temp_y_min)*0.06,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([2.1 2.9],temp_y_max+(temp_y_max-temp_y_min)*0.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0 4]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.3
    set(gca, 'FontSize', 12) %14
    set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = ["Correct", "Stimuli-error", "Response-error"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.17;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%9-->12
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',10);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('Location correlation'),'fontsize',12);
    temp_title.Interpreter = 'none';
    
end


temp_1 = r_mean_correct_theoretical_summary;
temp_1_chanceLevel = mean(r_mean_chanceLevel_summary);

[~,temp_p] = ttest(temp_1,temp_1_chanceLevel);

if if_plot == 1 && if_plot_part2 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1 336*1.11*0.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.85 336*1.11*0.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35 295*0.78*0.94*0.6*0.965]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_y_min = min([temp_1 temp_1_chanceLevel]);
    %temp_y_max = max(temp_1);
    temp_y_max = 1;
    
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    hold on
    
    
    plot([0 2],temp_1_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    set(gca,'YTick',0:0.2:1,'FontSize', 9);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('Location correlation'),'fontsize',9);
    %subtitle(sprintf('Correct trials'),'FontSize',7,'FontWeight','normal');
    subtitle(sprintf('theoretical & decoder'),'FontSize',7,'FontWeight','normal');
    
    temp_title.Interpreter = 'none';
    
end



temp_1 = r_mean_correct_summary;
temp_1_chanceLevel = mean(r_mean_chanceLevel_summary);

[~,temp_p] = ttest(temp_1,temp_1_chanceLevel);

if if_plot == 1 && if_plot_part2 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1 336*1.11*0.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.85 336*1.11*0.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35 295*0.78*0.94*0.6*0.965]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_y_min = min([temp_1 temp_1_chanceLevel]);
    %temp_y_max = max(temp_1);
    temp_y_max = 1;
    
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    hold on
    
    
    plot([0 2],temp_1_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    set(gca,'YTick',0:0.2:1,'FontSize', 9);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('Location correlation'),'fontsize',9);
    %subtitle(sprintf('Correct trials'),'FontSize',7,'FontWeight','normal');
    subtitle(sprintf('behavior & decoder'),'FontSize',7,'FontWeight','normal');
    
    temp_title.Interpreter = 'none';
    
end


%% Part 2 additional: Compare r_mean_correct_theoretical_summary & r_mean_correct_summary
temp_1 = r_mean_correct_theoretical_summary';
temp_2 = r_mean_correct_summary';

[~,temp_p] = ttest(r_mean_correct_theoretical_summary,r_mean_correct_summary);

if if_plot == 1 && if_plot_part2 == 1
    
    if_temp_violinplot0_pairline1 = 1;
    
    fig = figure('Name','theoretical & correct','NumberTitle','off');
    
    set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35 295*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    %temp_y_min = min([temp_1;temp_2]);
    %temp_y_max = max([temp_1;temp_2]);
    temp_y_min = 0;
    temp_y_max = 1;
    
    if if_temp_violinplot0_pairline1 == 0
        temp_data = [temp_1;temp_2];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        
        temp_label = [g1;g2];
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 2, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        
    elseif if_temp_violinplot0_pairline1 == 1
        
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);%30
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);%30
        hold on
        
    end
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);
    set(gca, 'FontSize', 8);%10
    
    xticks([1 2]);
    
    xtl = ["Theoretical"; "Behavior"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.22;%0.56,0.4
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9
    set(gca,'xticklabel','');
    
    %yticks([0 1]);
    yticks([0:0.2:1]);
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');%10
    
    title(sprintf('Location correlation'),'FontSize',9,'FontWeight','normal');
    subtitle(sprintf('With decoder'),'FontSize',7,'FontWeight','normal');
    
end


%% Part 2H: Sequence-level memory precision evidence

% memory correct trials
temptemp1 = allPart_summary_mulFOV.part2_summary_mulFOV.r_precision_summary;

% memory correct + error trials
temptemp2 = allPart_summary_mulFOV.part2_summary_mulFOV.r_precision_allMemory_summary;

% noChoice correct + error trials
temptemp3 = allPart_summary_mulFOV.part2_summary_mulFOV.r_precision_allNoChoice_summary;

r_precision_summary = temptemp3;

if if_entropy_behaviorInverted == 0
    r_precision_summary = abs(r_precision_summary) * -1;
end

temp_1 = r_precision_summary;

temp_1_chanceLevel = 0;

[~,temp_p] = ttest(temp_1,temp_1_chanceLevel);

if if_plot == 1 && if_plot_part2 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1 336*1.11*0.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.85 336*1.11*0.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35 295*0.78*0.94*0.6*0.965]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    %temp_1 = r_precision_summary;
    
    %temp_y_min = min(temp_1);
    %temp_y_max = max(temp_1);
    
    if if_entropy_behaviorInverted == 1
        temp_y_min = min([temp_1 temp_1_chanceLevel]);
        temp_y_max = 1;
    elseif if_entropy_behaviorInverted == 0
        temp_y_min = -1;
        temp_y_max = 0;
    end
    
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    
    
    plot([0 2],temp_1_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    if if_entropy_behaviorInverted == 1
        text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
    elseif if_entropy_behaviorInverted == 0
        text(1,temp_y_min-(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
    end
    
    set(gca,'linewidth',1.5)
    
    if if_entropy_behaviorInverted == 1
        set(gca,'YTick',0:0.2:1,'FontSize', 9);%12
        
    elseif if_entropy_behaviorInverted == 0
        set(gca,'ydir','reverse');
        set(gca,'YTick',-1:0.2:0,'FontSize', 9);%12
    end
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.1:1,'FontSize', 12);%12
    %set(gca,'YTick',0:0.2:1,'FontSize', 9);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('WM strength'),'fontsize',9);
    subtitle(sprintf('Across seqs'),'FontSize',7,'FontWeight','normal');
    
    temp_title.Interpreter = 'none';
    
end


%% Part 2I, Right: Trial-level memory precision evidence
seqAccuracy_memoryPrecisionLow_allMemoryStimuli_mean_summary = allPart_summary_mulFOV.part2_summary_mulFOV.seqAccuracy_memoryPrecisionLow_allMemoryStimuli_mean_summary;
seqAccuracy_memoryPrecisionHigh_allMemoryStimuli_mean_summary = allPart_summary_mulFOV.part2_summary_mulFOV.seqAccuracy_memoryPrecisionHigh_allMemoryStimuli_mean_summary;

temp_1 = seqAccuracy_memoryPrecisionLow_allMemoryStimuli_mean_summary';
temp_2 = seqAccuracy_memoryPrecisionHigh_allMemoryStimuli_mean_summary';
[~,temp_p,~,~] = ttest(temp_1,temp_2);

if if_plot == 1 && if_plot_part2 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 450 355*0.5*0.95 200*1.8*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35 295*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    %temp_y_min = min([temp_1;temp_2]);
    %temp_y_max = max([temp_1;temp_2]);
    
    temp_y_min = 0;
    temp_y_max = 1;
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.25 2.75])%[0.15 2.85]
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    xtl = ["Low-strength"; "High-strength"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %     if if_monkey_D0_Z1 == 0
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.185;%0.56,0.4
    %     elseif if_monkey_D0_Z1 == 1
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.425;%0.56,0.4
    %     end
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.21;
    
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    set(gca,'YTick',0:0.2:1,'FontSize', 9);%12
    
    ylabel('Accuracy', 'FontSize', 9, 'FontWeight', 'normal');
    %title(sprintf('Stimuli-labeled'),'fontsize',9);
    %subtitle(sprintf('Seq accuracy'),'FontSize',7,'FontWeight','normal');
    title(sprintf('Within seqs'),'fontsize',9);
    subtitle(sprintf('All seqs'),'FontSize',7,'FontWeight','normal');
    
end



seqAccuracy_memoryPrecisionLow_allMemoryResponse_mean_summary = allPart_summary_mulFOV.part2_summary_mulFOV.seqAccuracy_memoryPrecisionLow_allMemoryResponse_mean_summary;
seqAccuracy_memoryPrecisionHigh_allMemoryResponse_mean_summary = allPart_summary_mulFOV.part2_summary_mulFOV.seqAccuracy_memoryPrecisionHigh_allMemoryResponse_mean_summary;

temp_1 = seqAccuracy_memoryPrecisionLow_allMemoryResponse_mean_summary';
temp_2 = seqAccuracy_memoryPrecisionHigh_allMemoryResponse_mean_summary';
[~,temp_p,~,~] = ttest(temp_1,temp_2);

if if_plot == 1 && if_plot_part2 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 450 355*0.5*0.95 200*1.8*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35 295*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    %temp_y_min = min([temp_1;temp_2]);
    %temp_y_max = max([temp_1;temp_2]);
    
    temp_y_min = 0;
    temp_y_max = 1;
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.25 2.75])%[0.15 2.85]
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    xtl = ["Low-strength"; "High-strength"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.21;%0.56,0.4
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
    
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    set(gca,'YTick',0:0.2:1,'FontSize', 9);%12
    
    ylabel('Accuracy', 'FontSize', 9, 'FontWeight', 'normal');
    title(sprintf('Response-labeled'),'fontsize',9);
    subtitle(sprintf('Seq accuracy'),'FontSize',7,'FontWeight','normal');
    
end


%% (New) Part 2I, Right: Trial-level memory precision evidence
if if_plot == 1 && if_plot_part2 == 1
    
    memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary = allPart_summary_mulFOV.part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary;
    memoryPrecision_allSeq_allMemoryStimuliErrorMean_mean_summary = allPart_summary_mulFOV.part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryStimuliErrorMean_mean_summary;
    
    temp_1 = memoryPrecision_allSeq_allMemoryStimuliErrorMean_mean_summary';
    temp_2 = memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary';
    % [~,temp_p,~,~] = ttest(temp_1,temp_2);
    [~,temp_p,~,~] = ttest(temp_1,temp_2,'Tail','left');
end
if if_plot == 1 && if_plot_part2 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35 295*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 10)
    
    xticks([1 2]);
    xtl = ["Error"; "Correct"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    if if_monkey_D0_Z1 == 0
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.205;%0.56,0.4
    elseif if_monkey_D0_Z1 == 1
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.555;%0.56,0.4
    end
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Memory strength', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Stimuli-labeled'),'fontsize',10);
end

if if_plot == 1 && if_plot_part2 == 1
    
    memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary = allPart_summary_mulFOV.part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary;
    memoryPrecision_allSeq_allMemoryResponseErrorMean_mean_summary = allPart_summary_mulFOV.part2_summary_mulFOV.memoryPrecision_allSeq_allMemoryResponseErrorMean_mean_summary;
    
    temp_1 = memoryPrecision_allSeq_allMemoryResponseErrorMean_mean_summary';
    temp_2 = memoryPrecision_allSeq_allMemoryCorrectMean_mean_summary';
    % [~,temp_p,~,~] = ttest(temp_1,temp_2);
    [~,temp_p,~,~] = ttest(temp_1,temp_2,'Tail','left');
end
if if_plot == 1 && if_plot_part2 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35 295*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 10)
    
    xticks([1 2]);
    xtl = ["Error"; "Correct"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    if if_monkey_D0_Z1 == 0
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.205;%0.56,0.4
    elseif if_monkey_D0_Z1 == 1
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.555;%0.56,0.4
    end
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Memory strength', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Response-labeled'),'fontsize',10);
end




%% Part 3C: Sequence-level meta-memory evidence
AUROC_meta_delay1_summary = allPart_summary_mulFOV.part3_summary_mulFOV.AUROC_meta_delay1_summary;
AUROC_chanceLevel = 0.5;

[~,temp_p] = ttest(AUROC_meta_delay1_summary,AUROC_chanceLevel);

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.85*0.9 336*1.11*0.5*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 161.57*0.87 336*1.11*0.5*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50+0 50+0 161.57*0.87*0.85 336*1.11*0.5*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_1 = AUROC_meta_delay1_summary;
    
    temp_y_min = min([temp_1 AUROC_chanceLevel]);
    %temp_y_max = max(temp_1);
    temp_y_max = 1;
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    hold on
    
    
    plot([0 2],AUROC_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    set(gca,'YTick',0:0.1:1,'FontSize', 9);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('AUROC', 'FontSize', 9, 'FontWeight', 'normal');
    %temp_title = title(sprintf('Meta-memory AUROC'),'fontsize',9);
    temp_title = title(sprintf('Meta-WM\ndecoder'),'fontsize',8);
    temp_title.Interpreter = 'none';
    
end


r_meta_seqLevel_choice_delay1_summary = allPart_summary_mulFOV.part3_summary_mulFOV.r_meta_seqLevel_choice_delay1_summary;
r_meta_seqLevel_noChoice_delay1_summary = allPart_summary_mulFOV.part3_summary_mulFOV.r_meta_seqLevel_noChoice_delay1_summary;

temp_1 = r_meta_seqLevel_choice_delay1_summary;
temp_2 = r_meta_seqLevel_noChoice_delay1_summary;

r_chanceLevel = 0;

[~,temp_p1] = ttest(temp_1,r_chanceLevel);
[~,temp_p2] = ttest(temp_2,r_chanceLevel);

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','locDistri','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.85 336*1.11*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.85*0.87 336*1.11*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.85*0.87*0.8 336*1.11*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    %temp_y_min = min([temp_1 temp_2]);
    %temp_y_max = max([temp_1 temp_2]);
    %temp_y_min = min([temp_1 temp_2 r_chanceLevel]);
    temp_y_min = -1;
    temp_y_max = max([temp_1 temp_2 r_chanceLevel]);
    
    temp_data = [temp_1';temp_2'];
    
    g1 = repmat({'A'},length(temp_1),1);
    g2 = repmat({'B'},length(temp_2),1);
    
    temp_label = [g1;g2];
    
    temptemp_color1 = [1 1 1]*0.5;
    temptemp_color2 = repmat(temptemp_color1, 2, 1);
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',[{'A'};{'B'}]);
    h(1).ViolinPlot.FaceAlpha = 0.1;
    h(2).ViolinPlot.FaceAlpha = 0.1;
    
    plot([0 3],r_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p1 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p1 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p1 < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p2 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p2 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p2 < 0.05
        tempTxt = sprintf('*');
    end
    text(2,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0 3]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    set(gca,'ydir','reverse');
    
    xtl = ["Free-choice", "Forced-to-test"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %     if if_monkey_D0_Z1 == 0
    %         %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.35;
    %         xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.32;
    %     elseif if_monkey_D0_Z1 == 1
    %         xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.34;
    %     end
    xtext_yp=temp_y_max*ones(1,length(xt))+(temp_y_max-temp_y_min)*0.15;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('Meta-WM\nAcross seqs'),'fontsize',8.5);
    temp_title.Interpreter = 'none';
    
end


%% Part 3: Trial-level correlation between WM strength & meta-WM
trialCorr_precisionMeta_choice_summary = allPart_summary_mulFOV.part3_summary_mulFOV.trialCorr_precisionMeta_choice_summary;
trialCorr_precisionMeta_choiceMemory_summary = allPart_summary_mulFOV.part3_summary_mulFOV.trialCorr_precisionMeta_choiceMemory_summary;
trialCorr_precisionMeta_choiceOffload_summary = allPart_summary_mulFOV.part3_summary_mulFOV.trialCorr_precisionMeta_choiceOffload_summary;

temp_1 = trialCorr_precisionMeta_choice_summary;
temp_2 = trialCorr_precisionMeta_choiceMemory_summary;
temp_3 = trialCorr_precisionMeta_choiceOffload_summary;

% temptempBoolIndex = ~isoutlier(allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_summary,'percentiles',[20 100]);
% temptempBoolIndex = ~isoutlier(allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_theoretical_summary,'mean','ThresholdFactor',1.2);
temptempBoolIndex = ~isoutlier(allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_theoretical_summary,'median','ThresholdFactor',1.5);
% temptempBoolIndex = ~isoutlier(allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_theoretical_summary,'quartiles','ThresholdFactor',0.7);
% temptempBoolIndex = ~isoutlier(allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_theoretical_summary,'grubbs','ThresholdFactor',0.8);
% temptempBoolIndex = ~isoutlier(allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_theoretical_summary,'gesd','ThresholdFactor',0.8);


% find(~temptempBoolIndex)

temp_1 = temp_1(temptempBoolIndex);
temp_2 = temp_2(temptempBoolIndex);
temp_3 = temp_3(temptempBoolIndex);

r_chanceLevel = 0;

[~,temp_p1] = ttest(temp_1,r_chanceLevel,'tail','right');
[~,temp_p2] = ttest(temp_2,r_chanceLevel,'tail','right');
[~,temp_p3] = ttest(temp_3,r_chanceLevel,'tail','right');

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','locDistri','NumberTitle','off');
    set(gcf,'Position',[50+250 50+0 125*1.4 336*1.11*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(3,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    temp_y_min = -0.15;%0,-0.15
    %temp_y_max = max([temp_1 temp_2 temp_3 r_chanceLevel]);
    temp_y_max = 0.375;%0.25
    
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
    
    plot([0 4],r_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p1 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p1 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p1 < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p2 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p2 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p2 < 0.05
        tempTxt = sprintf('*');
    end
    text(2,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p3 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p3 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p3 < 0.05
        tempTxt = sprintf('*');
    end
    text(3,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 4]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.2 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    %set(gca,'ydir','reverse');
    
    xtl = ["Choice","Memory","Offload"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %     if if_monkey_D0_Z1 == 0
    %         %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.35;
    %         xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.32;
    %     elseif if_monkey_D0_Z1 == 1
    %         xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.34;
    %     end
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.3;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('Across trials\n Strength VS. Meta'),'fontsize',8.5);
    temp_title.Interpreter = 'none';
    
end


%% Part 3: Trial-level correlation between WM strength & meta-WM (20000+ trials)
memoryPrecision_trialLevel_summary = allPart_summary_mulFOV.part5_summary_mulFOV.memoryPrecision_trialLevel_summary;
meta_trialLevel_summary = allPart_summary_mulFOV.part5_summary_mulFOV.meta_trialLevel_summary;

temptempBoolIndex_CMC_summary = allPart_summary_mulFOV.part4_summary_mulFOV.temptempBoolIndex_CMC_summary;
temptempBoolIndex_CF_summary = allPart_summary_mulFOV.part4_summary_mulFOV.temptempBoolIndex_CF_summary;
temptempBoolIndex_CME_summary = allPart_summary_mulFOV.part4_summary_mulFOV.temptempBoolIndex_CME_summary;

numFOV;
temptempBoolIndex_CMC_summary_collapsed = [];
temptempBoolIndex_CF_summary_collapsed = [];
temptempBoolIndex_CME_summary_collapsed = [];
for tempi=1:numFOV
    temptempBoolIndex_CMC_summary_collapsed = [temptempBoolIndex_CMC_summary_collapsed; temptempBoolIndex_CMC_summary{tempi}];
    temptempBoolIndex_CF_summary_collapsed = [temptempBoolIndex_CF_summary_collapsed; temptempBoolIndex_CF_summary{tempi}];
    temptempBoolIndex_CME_summary_collapsed = [temptempBoolIndex_CME_summary_collapsed; temptempBoolIndex_CME_summary{tempi}];   
end
temptempBoolIndex_CMC_summary_collapsed = temptempBoolIndex_CMC_summary_collapsed == 1;
temptempBoolIndex_CF_summary_collapsed = temptempBoolIndex_CF_summary_collapsed == 1;
temptempBoolIndex_CME_summary_collapsed = temptempBoolIndex_CME_summary_collapsed == 1;

temptempBoolIndex_CM_summary_collapsed = temptempBoolIndex_CMC_summary_collapsed | temptempBoolIndex_CME_summary_collapsed;

temptempBoolIndex_Choice_summary_collapsed = temptempBoolIndex_CM_summary_collapsed | temptempBoolIndex_CF_summary_collapsed;

[temp_r_choice,temp_p_choice] = corr(memoryPrecision_trialLevel_summary(temptempBoolIndex_Choice_summary_collapsed),meta_trialLevel_summary(temptempBoolIndex_Choice_summary_collapsed));
[temp_r_memory,temp_p_memory] = corr(memoryPrecision_trialLevel_summary(temptempBoolIndex_CM_summary_collapsed),meta_trialLevel_summary(temptempBoolIndex_CM_summary_collapsed));
[temp_r_offload,temp_p_offload] = corr(memoryPrecision_trialLevel_summary(temptempBoolIndex_CF_summary_collapsed),meta_trialLevel_summary(temptempBoolIndex_CF_summary_collapsed));


%% Part 3: Across seq correlation between offloading rate & meta-WM
r_meta_seqLevel_choice_delay1_summary = allPart_summary_mulFOV.part3_summary_mulFOV.r_meta_seqLevel_choice_delay1_summary;
r_meta_seqLevel_choiceMemory_delay1_summary = allPart_summary_mulFOV.part3_summary_mulFOV.r_meta_seqLevel_choiceMemory_delay1_summary;
r_meta_seqLevel_choiceOffload_delay1_summary = allPart_summary_mulFOV.part3_summary_mulFOV.r_meta_seqLevel_choiceOffload_delay1_summary;

temp_1 = r_meta_seqLevel_choice_delay1_summary;
temp_2 = r_meta_seqLevel_choiceMemory_delay1_summary;
temp_3 = r_meta_seqLevel_choiceOffload_delay1_summary;

% temptempBoolIndex = ~isoutlier(allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_summary,'percentiles',[20 100]);
% temp_1 = temp_1(temptempBoolIndex);
% temp_2 = temp_2(temptempBoolIndex);
% temp_3 = temp_3(temptempBoolIndex);

r_chanceLevel = 0;

% [~,temp_p1] = ttest(temp_1,r_chanceLevel,'tail','right');
% [~,temp_p2] = ttest(temp_2,r_chanceLevel,'tail','right');
% [~,temp_p3] = ttest(temp_3,r_chanceLevel,'tail','right');
[~,temp_p1] = ttest(temp_1,r_chanceLevel);
[~,temp_p2] = ttest(temp_2,r_chanceLevel);
[~,temp_p3] = ttest(temp_3,r_chanceLevel);

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','locDistri','NumberTitle','off');
    set(gcf,'Position',[50+250*2 50+0 125*1.4 336*1.11*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(3,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    temp_y_min = -1.1;%0,-0.15
    %temp_y_max = max([temp_1 temp_2 temp_3 r_chanceLevel]);
    temp_y_max = 0.1;%0.25
    
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
    
    plot([0 4],r_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p1 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p1 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p1 < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p2 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p2 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p2 < 0.05
        tempTxt = sprintf('*');
    end
    text(2,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p3 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p3 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p3 < 0.05
        tempTxt = sprintf('*');
    end
    text(3,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 4]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.2 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    set(gca,'ydir','reverse');
    
    xtl = ["Choice","Memory","Offload"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %     if if_monkey_D0_Z1 == 0
    %         %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.35;
    %         xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.32;
    %     elseif if_monkey_D0_Z1 == 1
    %         xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.34;
    %     end
    xtext_yp=temp_y_max*ones(1,length(xt))+(temp_y_max-temp_y_min)*0.3;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('Meta-WM\nAcross seqs'),'fontsize',8.5);
    temp_title.Interpreter = 'none';
    
end



%% Part 3D, Right: Trial-level meta-memory evidence
seqRProb_metaLow_choice_mean_summary = allPart_summary_mulFOV.part3_summary_mulFOV.seqRProb_metaLow_choice_mean_summary;
seqRProb_metaHigh_choice_mean_summary = allPart_summary_mulFOV.part3_summary_mulFOV.seqRProb_metaHigh_choice_mean_summary;

temp_1 = seqRProb_metaLow_choice_mean_summary';
temp_2 = seqRProb_metaHigh_choice_mean_summary';
[~,temp_p,~,~] = ttest(temp_1,temp_2);

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5*0.85 200*1.8*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*0.5*0.85*0.87 200*1.8*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    %temp_y_min = min([temp_1;temp_2]);
    %temp_y_max = max([temp_1;temp_2]);
    temp_y_min = 0;
    temp_y_max = 1;
    
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    xticks([1 2]);
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xtl = ["Low-meta"; "High-meta"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %     if if_monkey_D0_Z1 == 0
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.18;%0.56,0.4
    %     elseif if_monkey_D0_Z1 == 1
    %         %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.45;%0.56,0.4
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.24;%0.56,0.4
    %     end
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;
    
    
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%25
    
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Offloading rate', 'FontSize', 9, 'FontWeight', 'normal');
    %title(sprintf('Free-choice\nSeq Offloading rate'),'fontsize',8.5);
    title(sprintf('All seqs'),'fontsize',8.5);
    
end


%% (New) Part 3D, Right: Trial-level meta-memory evidence
if if_plot == 1 && if_plot_part3 == 1
    
    meta_allSeq_choiceOffloadMean_mean_summary = allPart_summary_mulFOV.part3_summary_mulFOV.meta_allSeq_choiceOffloadMean_mean_summary;
    meta_allSeq_choiceMemoryStimuliMean_mean_summary = allPart_summary_mulFOV.part3_summary_mulFOV.meta_allSeq_choiceMemoryStimuliMean_mean_summary;
    
    temp_1 = meta_allSeq_choiceOffloadMean_mean_summary';
    temp_2 = meta_allSeq_choiceMemoryStimuliMean_mean_summary';
    [~,temp_p,~,~] = ttest(temp_1,temp_2);
    %[~,temp_p,~,~] = ttest(temp_1,temp_2,'Tail','left');
    
end

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*0.5*0.85*0.87 200*1.8*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    %     temp_y_min = min([temp_1;temp_2]);
    %     temp_y_max = max([temp_1;temp_2]);
    temp_y_min = 0;
    temp_y_max = 1;
    
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    %xtl = ["Choice-offload"; "Choice-memory"];
    xtl = ["Offload"; "Memory"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %     if if_monkey_D0_Z1 == 0
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.655;%0.56,0.4
    %     elseif if_monkey_D0_Z1 == 1
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.425;%0.56,0.4
    %     end
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Meta-WM', 'FontSize', 9, 'FontWeight', 'normal');
    %title(sprintf('Choice'),'fontsize',8.5);
    title(sprintf('Within seqs'),'fontsize',8.5);
end



%% (New) Part 3D, Right: Trial-level meta-memory evidence
if if_plot == 1 && if_plot_part3 == 1
    
    memoryPrecision_allSeq_choiceOffloadMean_mean_summary = allPart_summary_mulFOV.part3_summary_mulFOV.memoryPrecision_allSeq_choiceOffloadMean_mean_summary;
    memoryPrecision_allSeq_choiceMemoryStimuliMean_mean_summary = allPart_summary_mulFOV.part3_summary_mulFOV.memoryPrecision_allSeq_choiceMemoryStimuliMean_mean_summary;
    
    temp_1 = memoryPrecision_allSeq_choiceOffloadMean_mean_summary';
    temp_2 = memoryPrecision_allSeq_choiceMemoryStimuliMean_mean_summary';
    [~,temp_p,~,~] = ttest(temp_1,temp_2);
    %[~,temp_p,~,~] = ttest(temp_1,temp_2,'Tail','left');
    
end

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[10+150 450 355*0.5*0.85*0.87 200*1.8*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    %     temp_y_min = min([temp_1;temp_2]);
    %     temp_y_max = max([temp_1;temp_2]);
    temp_y_min = 0;
    temp_y_max = 1;
    if if_monkey_D0_Z1 == 1
        temp_y_max = 0.3;
    end
    
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    xtl = ["Offload"; "Memory"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;

    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('WM strength', 'FontSize', 9, 'FontWeight', 'normal');
    title(sprintf('Within seqs'),'fontsize',8.5);
end



%% Part 3F: Single neuron tuning property of memory precision and meta-memory
if true
    if if_plot == 1 && if_plot_part3 == 1
        fig = figure('Name','tuning property (seq-level)','NumberTitle','off');
        %set(gcf,'Position',[400 50 308.4*0.85 200*0.92*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[400 50 308.4*0.85*0.72 200*0.92*0.94*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        x = r2_6loc_summary;
        y = r2_rProb_summary;
        
        [temp_r, temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        
        x_min = -0.032;
        x_max = 0.4834;
        y_min = -0.0185;
        y_max = 0.7;
        
        temptempBoolIndex_A = ~temptempBoolIndex_7;
        temptempBoolIndex_B = temptempBoolIndex_7;
        temptempBoolIndex_C = ~temptempBoolIndex_5;
        temptempBoolIndex_D =  temptempBoolIndex_5;
        
        temptempBoolIndex_AC = temptempBoolIndex_A & temptempBoolIndex_C;% no selective
        temptempBoolIndex_AD = temptempBoolIndex_A & temptempBoolIndex_D;% pure meta
        temptempBoolIndex_BC = temptempBoolIndex_B & temptempBoolIndex_C;% pure location
        temptempBoolIndex_BD = temptempBoolIndex_B & temptempBoolIndex_D;% mixed selective
        
        
        temp_MarkerFaceAlpha = 0.5;%0.5,1,0.75,1,0.5
        temp_LineWidth = 0.75;%0.5
        
        temp_size = 4;%10,6
        
        %         scatter(x(temptempBoolIndex_AC),y(temptempBoolIndex_AC),temp_size,...
        %            'filled','MarkerFaceColor',[1 1 1]*0.85,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        scatter(x(temptempBoolIndex_AC),y(temptempBoolIndex_AC),temp_size,...
            'filled','MarkerFaceColor',[1 1 1]*0.4,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        hold on
        
        %         scatter(x(temptempBoolIndex_AD),y(temptempBoolIndex_AD),temp_size,...
        %            'filled','MarkerFaceColor',[1 1 1]*0.5,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        scatter(x(temptempBoolIndex_AD),y(temptempBoolIndex_AD),temp_size,...
            'filled','MarkerFaceColor',[1 1 1]*0.4,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        hold on
        
        %         scatter(x(temptempBoolIndex_BC),y(temptempBoolIndex_BC),temp_size,...
        %            'filled','MarkerFaceColor',[1 1 1]*1,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth,'MarkerEdgeColor',[1 1 1]*0);
        scatter(x(temptempBoolIndex_BC),y(temptempBoolIndex_BC),temp_size,...
            'filled','MarkerFaceColor',[1 1 1]*0.4,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        hold on
        
        %         scatter(x(temptempBoolIndex_BD),y(temptempBoolIndex_BD),temp_size,...
        %            'filled','MarkerFaceColor',[1 1 1]*0,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth,'MarkerEdgeColor',[1 1 1]*0);
        scatter(x(temptempBoolIndex_BD),y(temptempBoolIndex_BD),temp_size,...
            'filled','MarkerFaceColor',[1 1 1]*0.4,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        hold on
        
        
        n = 1;
        [p_mapping,S] = polyfit(x,y,n);
        x_fit = x_min:0.01:x_max;
        y_fit = polyval(p_mapping,x_fit);
        
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        
        
        
        text(x_max-(x_max-x_min)*0.245,y_min+(y_max-y_min)*(0.95-0*0.13),sprintf('n = %d',length(x)), 'fontsize',6.5,'HorizontalAlignment','left');
        
        
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        text(x_max-(x_max-x_min)*0.245,y_min+(y_max-y_min)*(0.95-1*0.13),sprintf('r = %.3f',temp_r), 'fontsize',6.5,'HorizontalAlignment','left');
        text(x_max-(x_max-x_min)*0.245,y_min+(y_max-y_min)*(0.95-2*0.13),sprintf('%s',tempTxt), 'fontsize',6.5,'HorizontalAlignment','left');
        
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        %xlabel('r2 (Memory)', 'FontSize', 9, 'FontWeight', 'normal');
        %ylabel('r2 (Meta-memory)', 'FontSize', 9, 'FontWeight', 'normal');
        xlabel('Memory (r2)', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Meta-memory (r2)', 'FontSize', 9, 'FontWeight', 'normal');
        %title(sprintf('Regression explained variance, n = %d',length(x)),'FontSize',6.5);
        title(sprintf('Regression explained variance'),'FontSize',8);
        %subtitle(sprintf('n = %d',length(x)),'FontSize',7);
        
    end
    
end


%% Part 3F: Single neuron tuning property of memory precision and meta-memory
if true
    if if_plot == 1 && if_plot_part3 == 1
        fig = figure('Name','tuning property (seq-level)','NumberTitle','off');
        %set(gcf,'Position',[400 50 308.4*0.85 200*0.92*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[600 50 308.4*0.85*0.72 200*0.92*0.94*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[600 50 308.4*0.85*0.72 200*0.92*0.94*0.8*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        temp_if_plot_selectiveNeurons = 0;
        
        x = beta_6loc_summary_norm;
        y = beta_gProb_summary_norm;
        %y = beta_gProb_summary;
        
        temptempBoolIndex_A = ~temptempBoolIndex_7;
        temptempBoolIndex_B = temptempBoolIndex_7;
        temptempBoolIndex_C = ~temptempBoolIndex_5;
        temptempBoolIndex_D =  temptempBoolIndex_5;
        
        temptempBoolIndex_AC = temptempBoolIndex_A & temptempBoolIndex_C;% no selective
        temptempBoolIndex_AD = temptempBoolIndex_A & temptempBoolIndex_D;% pure meta
        temptempBoolIndex_BC = temptempBoolIndex_B & temptempBoolIndex_C;% pure location
        temptempBoolIndex_BD = temptempBoolIndex_B & temptempBoolIndex_D;% mixed selective
        
        
        if temp_if_plot_selectiveNeurons == 1
            temptempBoolIndex = ~temptempBoolIndex_AC;
            x_valid = x(temptempBoolIndex);
            y_valid = y(temptempBoolIndex);
        elseif temp_if_plot_selectiveNeurons == 0
            x_valid = x;
            y_valid = y;
        end
        
        %[temp_r, temp_p] = corr(x,y);
        [temp_r, temp_p] = corr(x_valid,y_valid);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        
        x_min = 0;
        x_max = 0.8;
        y_min = 0;
        y_max = 2.3;%2.6
        
        
        temp_MarkerFaceAlpha = 0.5;%0.5,1,0.75,1,0.5
        temp_LineWidth = 0.75;%0.5
        
        temp_size = 4;%10,6
        
        scatter(x_valid,y_valid,temp_size,...
            'filled','MarkerFaceColor',[1 1 1]*0.4,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        hold on
        
        %         %         scatter(x(temptempBoolIndex_AC),y(temptempBoolIndex_AC),temp_size,...
        %         %            'filled','MarkerFaceColor',[1 1 1]*0.85,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        %         scatter(x(temptempBoolIndex_AC),y(temptempBoolIndex_AC),temp_size,...
        %             'filled','MarkerFaceColor',[1 1 1]*0.4,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        %         hold on
        %
        %         %         scatter(x(temptempBoolIndex_AD),y(temptempBoolIndex_AD),temp_size,...
        %         %            'filled','MarkerFaceColor',[1 1 1]*0.5,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        %         scatter(x(temptempBoolIndex_AD),y(temptempBoolIndex_AD),temp_size,...
        %             'filled','MarkerFaceColor',[1 1 1]*0.4,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        %         hold on
        %
        %         %         scatter(x(temptempBoolIndex_BC),y(temptempBoolIndex_BC),temp_size,...
        %         %            'filled','MarkerFaceColor',[1 1 1]*1,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth,'MarkerEdgeColor',[1 1 1]*0);
        %         scatter(x(temptempBoolIndex_BC),y(temptempBoolIndex_BC),temp_size,...
        %             'filled','MarkerFaceColor',[1 1 1]*0.4,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        %         hold on
        %
        %         %         scatter(x(temptempBoolIndex_BD),y(temptempBoolIndex_BD),temp_size,...
        %         %            'filled','MarkerFaceColor',[1 1 1]*0,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth,'MarkerEdgeColor',[1 1 1]*0);
        %         scatter(x(temptempBoolIndex_BD),y(temptempBoolIndex_BD),temp_size,...
        %             'filled','MarkerFaceColor',[1 1 1]*0.4,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        %         hold on
        
        
        n = 1;
        [p_mapping,S] = polyfit(x_valid,y_valid,n);
        x_fit = x_min:0.01:x_max;
        y_fit = polyval(p_mapping,x_fit);
        
        plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
        hold on
        
        
        
        %text(x_max-(x_max-x_min)*1.005,y_min+(y_max-y_min)*(0.95-0*0.13),sprintf('n = %d',length(x)), 'fontsize',6.5,'HorizontalAlignment','left');
        text(x_max-(x_max-x_min)*1.005,y_min+(y_max-y_min)*(0.95-0*0.13),sprintf('n = %d',length(x_valid)), 'fontsize',6.5,'HorizontalAlignment','left');
        
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        text(x_max-(x_max-x_min)*1.005,y_min+(y_max-y_min)*(0.95-1*0.13),sprintf('r = %.3f',temp_r), 'fontsize',6.5,'HorizontalAlignment','left');
        text(x_max-(x_max-x_min)*1.005,y_min+(y_max-y_min)*(0.95-2*0.13),sprintf('%s',tempTxt), 'fontsize',6.5,'HorizontalAlignment','left');
        
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('WM (beta)', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Meta-WM (beta)', 'FontSize', 9, 'FontWeight', 'normal');
        %title(sprintf('Regression coefficients'),'FontSize',8);
        
    end
    
end


if false
    if if_plot == 1 && if_plot_part3 == 1
        fig = figure('Name','tuning property (seq-level)','NumberTitle','off');
        %set(gcf,'Position',[400 50 445*0.8*1.02 379*0.7*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 50 720 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 50 445*0.8*1.02*0.9 379*0.7*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[400 50 445*0.8*1.02*0.9*0.9 379*0.7*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        if if_tuning_location0_precision1 == 1
            x = r2_seqPrecision_summary;
        elseif if_tuning_location0_precision1 == 0
            x = r2_6loc_summary;
        end
        y = r2_rProb_summary;
        
        [temp_r, temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        
        temp_size = 3;%10
        scatter(x,y,temp_size,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        if if_tuning_location0_precision1 == 1
            text(x_min+(x_max-x_min)*0.00,y_min+(y_max-y_min)*0.999,sprintf('r = %.3f',temp_r), 'fontsize',8,'FontWeight','normal');
            text(x_min+(x_max-x_min)*0.00,y_min+(y_max-y_min)*0.879,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        end
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        if if_tuning_location0_precision1 == 1
            xlabel('r2 of memory strength', 'FontSize', 12, 'FontWeight', 'bold');
        elseif if_tuning_location0_precision1 == 0
            xlabel('r2 of memory', 'FontSize', 12, 'FontWeight', 'bold');
        end
        ylabel('r2 of meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Multi-FOV (seq-level), n=%d\n Single neuron tuning property',length(x)),'FontSize',11);
        temp_title.Interpreter = 'none';
        
    end
end

%% Part 3G: Correlation between memory precision and meta-memory
r_precisionMeta_seqLevel_summary = allPart_summary_mulFOV.part3_summary_mulFOV.r_precisionMeta_seqLevel_summary;


r_chanceLevel = 0;
[~,temp_p] = ttest(r_precisionMeta_seqLevel_summary,r_chanceLevel);

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.15*0.85*0.9 336*1.11*0.5*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.15*0.85*0.9*0.87 336*1.11*0.5*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50+0 50+0 240*0.80*1.15*0.85*0.9*0.87*0.9 336*1.11*0.5*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_1 = r_precisionMeta_seqLevel_summary;
    
    %temp_y_min = min(temp_1);
    %temp_y_max = max(temp_1);
    temp_y_min = 0;
    temp_y_max = 1;
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    
    
    plot([0 2],r_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.1:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    %temp_title = title(sprintf('Sequence-level evidence \n Memory precision & Meta '),'fontsize',8.5);
    %temp_title = title(sprintf('Strength & Meta\nAcross seqs'),'fontsize',8.5);
    temp_title = title(sprintf('Across seqs\nStrength VS. Meta'),'fontsize',8.5);
    temp_title.Interpreter = 'none';
    
end

%% (New) Part 3: memoryPrecision & meta trialLevelEvidence in stimuli & response-labeled trials
temp_1 = allPart_summary_mulFOV.part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceStimuli_mean_summary;
temp_2 = allPart_summary_mulFOV.part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceStimuli_mean_summary;
[~,temp_p12] = ttest(temp_1,temp_2);

temp_3 = allPart_summary_mulFOV.part3_summary_mulFOV.seqMeta_memoryPrecisionLow_choiceResponse_mean_summary;
temp_4 = allPart_summary_mulFOV.part3_summary_mulFOV.seqMeta_memoryPrecisionHigh_choiceResponse_mean_summary;
[~,temp_p34] = ttest(temp_3,temp_4);

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','memory precsion pdf (mildSeq)','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*1.5*1.05*0.7*0.9 295]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2 295*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,2,'TileSpacing','compact','Padding','loose');
    
    
    % Compare meta-memory of lowPrecision & highPrecision, stimuli-labeled trials
    nexttile
    
    temp_p = temp_p12;
    
    %temp_y_min = min([temp_1;temp_2]);
    %temp_y_max = max([temp_1;temp_2]);
    temp_y_min = 0.375;%0
    temp_y_max = 0.625;%1
    
    if if_plot_violinplot0_pairline1 == 0
        temp_data = [temp_1;temp_2];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        
        temp_label = [g1;g2];
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 2, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        
    elseif if_plot_violinplot0_pairline1 == 1
        
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        
    end
    
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.65])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    xtl = ["Low-strength"; "High-strength"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    % xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.405;%0.56,0.4
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;%0.56,0.4
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
    set(gca,'xticklabel','');
    
    
    % yticks([0 1]);
    yticks([0:0.2:1]);
    
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
    % title(sprintf('Stimuli-labeled\nAll seqs'),'fontsize',9);
    title(sprintf('Stimuli-labeled'),'fontsize',9);
    
    
    
    % Compare meta-memory of lowPrecision & highPrecision, response-labeled trials
    nexttile
    
    temp_p = temp_p34;
    
    temp_1 = temp_3;
    temp_2 = temp_4;
    
    %temp_y_min = min([temp_1;temp_2]);
    %temp_y_max = max([temp_1;temp_2]);
    temp_y_min = 0.375;%0
    temp_y_max = 0.625;%1
    
    if if_plot_violinplot0_pairline1 == 0
        temp_data = [temp_1;temp_2];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        
        temp_label = [g1;g2];
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 2, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        
    elseif if_plot_violinplot0_pairline1 == 1
        
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        
    end
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.65])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    xtl = ["Low-strength"; "High-strength"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    % xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.405;%0.56,0.4
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;%0.56,0.4
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
    set(gca,'xticklabel','');
    
    
    % yticks([0 1]);
    yticks([0:0.2:1]);
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
    % title(sprintf('Response-labeled\nAll seqs'),'fontsize',9);
    title(sprintf('Response-labeled'),'fontsize',9);
    
end


%% (New) Part 3: Correlation, memoryPrecision & meta trialLevelEvidence in stimuli & response-labeled trials
temp_1 = allPart_summary_mulFOV.part3_summary_mulFOV.seqCorr_precisionMeta_choiceStimuli_mean_summary';
temp_2 = allPart_summary_mulFOV.part3_summary_mulFOV.seqCorr_precisionMeta_choiceResponse_mean_summary';

temp_1_chanceLevel = 0;
temp_2_chanceLevel = 0;

[~,temp_p1] = ttest(temp_1,temp_1_chanceLevel);
[~,temp_p2] = ttest(temp_2,temp_2_chanceLevel);

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2*1.02 295*0.78*0.94*1.14]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*0.35*2*1.02*0.75*0.9 295*0.78*0.94*1.14]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 197.6*1.1*1.15*1.05 246.6*1.1*1.15*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,2,'TileSpacing','tight','Padding','loose');
    
    %t.Title.String = sprintf('Strength VS. Meta, seq-level');
    t.Title.String = sprintf('Strength VS. Meta, within seqs');
    t.Title.FontSize = 9;
    t.Title.Interpreter = 'none';
    
    %stimuli-labeled trials
    nexttile
    
    %temp_y_min = min([temp_1 temp_1_chanceLevel]);
    %temp_y_max = max(temp_1);
    temp_y_min = -0.05;%-0.5
    temp_y_max = 0.30;%1,0.25
    
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    hold on
    
    
    plot([0 2],temp_1_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p1 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p1 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p1 < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14,12
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    %set(gca,'YTick',0:1,'FontSize', 10);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    %title(sprintf('Precision VS. Meta'),'fontsize',9);
    %subtitle(sprintf('Seq-level, stimuli-labeled'),'fontsize',7);
    title(sprintf('Stimuli-labeled'),'fontsize',7);
    
    
    
    %response-labeled trials
    nexttile
    
    %temp_y_min = min([temp_2 temp_2_chanceLevel]);
    %temp_y_max = max(temp_2);
    temp_y_min = -0.05;%-0.5
    temp_y_max = 0.30;%1,0.25
    
    temp_data = temp_2';
    temp_label = repmat({'A'},length(temp_2),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    hold on
    
    
    plot([0 2],temp_2_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p2 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p2 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p2 < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14,12
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    %set(gca,'YTick',0:1,'FontSize', 10);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    %ylabel('Correlation', 'FontSize', 9, 'FontWeight', 'normal');
    %title(sprintf('Precision VS. Meta'),'fontsize',9);
    %subtitle(sprintf('Seq-level, response-labeled'),'fontsize',7);
    title(sprintf('Response-labeled'),'fontsize',7);
    
end

%% Part 3H, Right: Trial-level memory precision VS. offloading rate evidence
seqRProb_memoryPrecisionLow_choiceResponse_mean_summary = allPart_summary_mulFOV.part3_summary_mulFOV.seqRProb_memoryPrecisionLow_choiceResponse_mean_summary;
seqRProb_memoryPrecisionHigh_choiceResponse_mean_summary = allPart_summary_mulFOV.part3_summary_mulFOV.seqRProb_memoryPrecisionHigh_choiceResponse_mean_summary;

temp_1 = seqRProb_memoryPrecisionLow_choiceResponse_mean_summary';
temp_2 = seqRProb_memoryPrecisionHigh_choiceResponse_mean_summary';
[~,temp_p,~,~] = ttest(temp_1,temp_2);

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5*0.85 200*1.8*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*0.5*0.85*0.87 200*1.8*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    %temp_y_min = min([temp_1;temp_2]);
    %temp_y_max = max([temp_1;temp_2]);
    temp_y_min = 0;
    temp_y_max = 1;
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    xticks([1 2]);
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xtl = ["Low-strength"; "High-strength"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.225;%0.56,0.4
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
    
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Offloading rate', 'FontSize', 9, 'FontWeight', 'normal');
    %title(sprintf('Response-labeled\nSeq Offloading rate'),'fontsize',8.5);
    title(sprintf('Response-labeled'),'fontsize',8.5);
end


%% (New) Part 3H, Right: Trial-level memory precision VS. offloading rate evidence
if if_plot == 1 && if_plot_part3 == 1
    
    memoryPrecision_allSeq_choiceOffloadMean_mean_summary = allPart_summary_mulFOV.part3_summary_mulFOV.memoryPrecision_allSeq_choiceOffloadMean_mean_summary;
    memoryPrecision_allSeq_choiceMemoryResponseMean_mean_summary = allPart_summary_mulFOV.part3_summary_mulFOV.memoryPrecision_allSeq_choiceMemoryResponseMean_mean_summary;
    
    temp_1 = memoryPrecision_allSeq_choiceOffloadMean_mean_summary';
    temp_2 = memoryPrecision_allSeq_choiceMemoryResponseMean_mean_summary';
    % [~,temp_p,~,~] = ttest(temp_1,temp_2);
    [~,temp_p,~,~] = ttest(temp_1,temp_2,'Tail','left');
    
end

if if_plot == 1 && if_plot_part3 == 1
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[10 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 10)
    
    xticks([1 2]);
    xtl = ["Choice-offload"; "Choice-memory"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    if if_monkey_D0_Z1 == 0
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.295;%0.56,0.4
    elseif if_monkey_D0_Z1 == 1
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.385;%0.56,0.4
    end
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Memory strength', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Response-labeled'),'fontsize',10);
end


%% Part 4B, Left: Trial number of OverMismatch, HighMatch, LowMatch, UnderMismatch
trialNum_memoryPrecisionLow_metaLow_choice_summary = allPart_summary_mulFOV.part4_summary_mulFOV.trialNum_memoryPrecisionLow_metaLow_choice_summary;
trialNum_memoryPrecisionHigh_metaHigh_choice_summary = allPart_summary_mulFOV.part4_summary_mulFOV.trialNum_memoryPrecisionHigh_metaHigh_choice_summary;
trialNum_memoryPrecisionLow_metaHigh_choice_summary = allPart_summary_mulFOV.part4_summary_mulFOV.trialNum_memoryPrecisionLow_metaHigh_choice_summary;
trialNum_memoryPrecisionHigh_metaLow_choice_summary = allPart_summary_mulFOV.part4_summary_mulFOV.trialNum_memoryPrecisionHigh_metaLow_choice_summary;

temp_1 = trialNum_memoryPrecisionLow_metaHigh_choice_summary;
temp_2 = trialNum_memoryPrecisionHigh_metaHigh_choice_summary;
temp_3 = trialNum_memoryPrecisionLow_metaLow_choice_summary;
temp_4 = trialNum_memoryPrecisionHigh_metaLow_choice_summary;

trialNum_all = temp_1+temp_2+temp_3+temp_4;

temp_1 = temp_1./trialNum_all;
temp_2 = temp_2./trialNum_all;
temp_3 = temp_3./trialNum_all;
temp_4 = temp_4./trialNum_all;


temp_chanceLevel = 0;

[~,temp_p1] = ttest(temp_1,temp_chanceLevel);
[~,temp_p2] = ttest(temp_2,temp_chanceLevel);
[~,temp_p3] = ttest(temp_3,temp_chanceLevel);
[~,temp_p4] = ttest(temp_4,temp_chanceLevel);


color_choiceMemoryLow = [1,133,113]/255; %[5,113,176]/255
color_choiceMemoryHigh = [128,205,193]/255; %[146,197,222]/255
color_choiceOffloadLow = [223,194,125]/255; %[244,165,130]/255
color_choiceOffloadHigh = [166,97,26]/255; %[202,0,32]/255

if false
    if if_plot == 1 && if_plot_part4 == 1
        fig = figure('Name','locDistri','NumberTitle','off');
        set(gcf,'Position',[50+0 50+0 240*0.80 336*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        %temp_y_min = min([temp_1 temp_2 temp_3 temp_4]);
        %temp_y_max = max([temp_1 temp_2 temp_3 temp_4]);
        temp_y_min = 0;
        temp_y_max = 0.4;
        
        temp_data = [temp_1';temp_2';temp_3';temp_4'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        g4 = repmat({'D'},length(temp_4),1);
        
        temp_label = [g1;g2;g3;g4];
        
        %temptemp_color1 = [1 1 1]*0.5;
        %temptemp_color2 = repmat(temptemp_color1, 4, 1);
        
        temptemp_color2 = ...
            [color_choiceMemoryLow;
            color_choiceMemoryHigh;
            color_choiceOffloadLow;
            color_choiceOffloadHigh];
        
        %violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'ViolinAlpha',0.5,...
        %    'BoxColor',[1 1 1]*0.2,...
        %    'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.5,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
        h(1).ViolinPlot.FaceAlpha = 0.3;
        h(2).ViolinPlot.FaceAlpha = 0.3;
        h(3).ViolinPlot.FaceAlpha = 0.3;
        h(4).ViolinPlot.FaceAlpha = 0.3;
        
        plot([0 5],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
        hold on
        
        tempTxt = sprintf('');
        if temp_p1 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p1 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p1 < 0.05
            tempTxt = sprintf('*');
        end
        text(1,temp_y_max+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p2 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p2 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p2 < 0.05
            tempTxt = sprintf('*');
        end
        text(2,temp_y_max+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p3 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p3 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p3 < 0.05
            tempTxt = sprintf('*');
        end
        text(3,temp_y_max+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        tempTxt = sprintf('');
        if temp_p4 < 0.001
            tempTxt = sprintf('***');
        elseif temp_p4 < 0.01
            tempTxt = sprintf('**');
        elseif temp_p4 < 0.05
            tempTxt = sprintf('*');
        end
        text(4,temp_y_max+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        
        xlim([0 5]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.3
        set(gca, 'FontSize', 8) %14
        %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = ["Over-mismatch", "High-match", "Low-match", "Under-mismatch"];
        %xtl = ["", "", "", ""];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
        xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.15;
        
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',6.5);%9-->12
        
        set(gca,'xticklabel','');
        
        
        ylabel('Trial proportion', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('Location correlation'),'fontsize',12);
        %temp_title.Interpreter = 'none';
        
    end
end

trialNum_memoryPrecisionLow_metaHigh_choice_summary;% Over-mismatch
trialNum_memoryPrecisionHigh_metaHigh_choice_summary;% High-match
trialNum_memoryPrecisionLow_metaLow_choice_summary;% Low-match
trialNum_memoryPrecisionHigh_metaLow_choice_summary;% Under-mismatch

% Match
temp_1 = trialNum_memoryPrecisionHigh_metaHigh_choice_summary + ...
    trialNum_memoryPrecisionLow_metaLow_choice_summary;
temp_1 = temp_1';

% Mismatch
temp_2 = trialNum_memoryPrecisionLow_metaHigh_choice_summary+ ...
    trialNum_memoryPrecisionHigh_metaLow_choice_summary;
temp_2 = temp_2';


trialNum_all = temp_1+temp_2;

temp_1 = temp_1./trialNum_all;
temp_2 = temp_2./trialNum_all;


[~,temp_p,~,~] = ttest(temp_1,temp_2);

if false
    if if_plot == 1 && if_plot_part4 == 1
        fig = figure('Name','','NumberTitle','off');
        set(gcf,'Position',[10 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        
        nexttile
        
        %temp_y_min = min([temp_1;temp_2]);
        %temp_y_max = max([temp_1;temp_2]);
        temp_y_min = 0.3;
        temp_y_max = 0.7;
        
        
        if if_plot_violinplot0_pairline1 == 0
            violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
            
        elseif if_plot_violinplot0_pairline1 == 1
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
        end
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 8)
        
        xticks([1 2]);
        xtl = ["Match"; "Mismatch"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %     if if_monkey_D0_Z1 == 0
        %         %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.405;%0.56,0.4
        %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;%0.56,0.4
        %     elseif if_monkey_D0_Z1 == 1
        %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.185;%0.56,0.4
        %     end
        xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;
        
        %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
        
        set(gca,'xticklabel','');
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Trial proportion', 'FontSize', 9, 'FontWeight', 'normal');
        %title(sprintf('Stimuli-labeled trial\nSeq accuracy'),'fontsize',10);
    end
end

%% new Part 4, signifTimeStampCount
signifTimeStampCount_memory_mean_summary = allPart_summary_mulFOV.part4_summary_mulFOV.signifTimeStampCount_memory_mean_summary;
signifTimeStampCount_meta_mean_summary = allPart_summary_mulFOV.part4_summary_mulFOV.signifTimeStampCount_meta_mean_summary;

temp_1 = signifTimeStampCount_memory_mean_summary';
temp_2 = signifTimeStampCount_meta_mean_summary';

temp_1 = (temp_1-1).*33.3;
temp_2 = (temp_2-1).*33.3;

% temptempBoolIndex = (~isoutlier(temp_1,'median')) & (~isoutlier(temp_2,'median'));%good
% temptempBoolIndex = (~isoutlier(temp_1,'mean')) & (~isoutlier(temp_2,'mean'));
temptempBoolIndex = (~isoutlier(temp_1,'quartiles')) & (~isoutlier(temp_2,'quartiles'));%good!!!
% temptempBoolIndex = (~isoutlier(temp_1,'grubbs')) & (~isoutlier(temp_2,'grubbs'));
% temptempBoolIndex = (~isoutlier(temp_1,'gesd')) & (~isoutlier(temp_2,'gesd'));%good!!!

% sum(temptempBoolIndex)

temp_1 = temp_1(temptempBoolIndex);
temp_2 = temp_2(temptempBoolIndex);

% [~,temp_p] = ttest(temp_1,temp_2,'tail','left'); %#ok<*UNRCH>
[~,temp_p] = ttest(temp_1,temp_2); %#ok<*UNRCH>
% temp_p

if if_plot_part4 == 1
    fprintf('signifTimeStampCount_memory: mean = %.3f, std = %.3f.\n',mean(temp_1),std(temp_1));
    fprintf('signifTimeStampCount_meta: mean = %.3f, std = %.3f.\n',mean(temp_2),std(temp_2));
end

if if_plot == 1 && if_plot_part4 == 1
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[210 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    %temp_y_min = min([temp_1;temp_2]);
    %temp_y_max = max([temp_1;temp_2]);
    temp_y_min = 0;
    temp_y_max = 700;%600
    
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.35 2.65])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    xtl = ["WM"; "Meta"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
    
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Latency time (ms)', 'FontSize', 9, 'FontWeight', 'normal');
    %title(sprintf('asdasd'),'fontsize',9);
end


%% new Part 4, proj comparison
if if_plot == 1 && if_plot_part4 == 1
    
    proj12_summary = allPart_summary_mulFOV.part4_summary_mulFOV.proj12_summary;
    proj3_summary = allPart_summary_mulFOV.part4_summary_mulFOV.proj3_summary;
    
    proj12_summary_1d = nan(1,numFOV);
    proj3_summary_1d = nan(1,numFOV);
    
    for tempi=1:numFOV
        proj12_summary_1d(tempi) = mean(proj12_summary{tempi});
        proj3_summary_1d(tempi) = mean(proj3_summary{tempi});
    end
    
    temp_1 = proj3_summary_1d';
    temp_2 = proj12_summary_1d';
    
    temp_1_mean = mean(temp_1);
    temp_2_mean = mean(temp_2);
    
    [~,temp_p] = ttest(temp_1,temp_2); %#ok<*UNRCH>
    
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[410 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    %temp_y_min = 0.03;%0.15
    %temp_y_max = 0.822;%0.6
    
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.35 2.65])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    xtl = ["Error"; "Correct+offload"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
    
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Mean', 'FontSize', 9, 'FontWeight', 'normal');
    title(sprintf('Projection comparison'),'fontsize',9);
end



%% new Part 4, signifTimeStampCount
overTrialProp_summary = allPart_summary_mulFOV.part4_summary_mulFOV.overTrialProp_summary;
overTrialProp_resampled_mean_summary = allPart_summary_mulFOV.part4_summary_mulFOV.overTrialProp_resampled_mean_summary;

temp_1 = overTrialProp_summary';
temp_2 = overTrialProp_resampled_mean_summary';

[~,temp_p] = ttest(temp_1,temp_2); %#ok<*UNRCH>

if if_plot == 1 && if_plot_part4 == 1
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[610 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    %temp_y_min = min([temp_1;temp_2]);
    %temp_y_max = max([temp_1;temp_2]);
    temp_y_min = 0.03;%0.15
    temp_y_max = 0.822;%0.6
    
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.35 2.65])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    xtl = ["Data"; "Shuffle"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.2;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
    
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Trial proportion', 'FontSize', 9, 'FontWeight', 'normal');
    title(sprintf('Low strength\nmismatch trials'),'fontsize',9);
end


%% Part 5B: Baseline meta-memory evidence
AUROC_meta_baseline_summary = allPart_summary_mulFOV.part5_summary_mulFOV.AUROC_meta_baseline_summary;
AUROC_chanceLevel = 0.5;

[~,temp_p] = ttest(AUROC_meta_baseline_summary,AUROC_chanceLevel);

if if_plot == 1 && if_plot_part5 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1 336*1.11*0.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.9 336*1.11*0.5*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.9*0.9 336*1.11*0.5*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_1 = AUROC_meta_baseline_summary;
    
    temp_y_min = min([temp_1 AUROC_chanceLevel]);
    temp_y_max = max(temp_1);
    
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    
    hold on
    
    
    plot([0 2],AUROC_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    set(gca,'YTick',0:0.1:1);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('AUROC', 'FontSize', 9);
    temp_title = title(sprintf('Baseline Meta\nAUROC'),'fontsize',9);
    temp_title.Interpreter = 'none';
    
end



p_meta_seqLevel_choice_baseline_summary = allPart_summary_mulFOV.part5_summary_mulFOV.p_meta_seqLevel_choice_baseline_summary;
r_meta_seqLevel_choice_baseline_summary = allPart_summary_mulFOV.part5_summary_mulFOV.r_meta_seqLevel_choice_baseline_summary;

p_chanceLevel = 0.05;
[~,temp_p] = ttest(p_meta_seqLevel_choice_baseline_summary,p_chanceLevel);

if if_plot == 1 && if_plot_part5 == 1
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[50+0 50+0 240*0.80*1.1*2*0.8*0.9*0.9 336*1.11*0.5*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,2,'TileSpacing','compact','Padding','loose');
    
    %t.Title.String = sprintf('Sequence-level \n correlation of baseline meta');
    t.Title.String = sprintf('Across seqs \n correlation of baseline meta');
    t.Title.FontSize = 9;
    t.Title.Interpreter = 'none';
    
    
    nexttile
    
    temp_1 = r_meta_seqLevel_choice_baseline_summary;
    
    temp_y_min = min(temp_1);
    temp_y_max = max(temp_1);
    
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    
    hold on
    
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.1:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Correlation', 'FontSize', 9);
    %temp_title = title(sprintf('Sequence-level evidence \n Baseline Meta'),'fontsize',11);
    %temp_title.Interpreter = 'none';
    
    
    nexttile
    
    temp_1 = p_meta_seqLevel_choice_baseline_summary;
    
    %temp_y_min = min(temp_1);
    temp_y_min = min([temp_1 p_chanceLevel]);
    temp_y_max = max(temp_1);
    
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    
    hold on
    
    plot([0 2],p_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    %text(0.3,p_chanceLevel+0.1,sprintf('%.2f',p_chanceLevel),'Color','black','FontSize',10,'FontWeight','normal',...
    %    'HorizontalAlignment','center');
    
    yticks([p_chanceLevel 0.5 1]);
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.1:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('P value', 'FontSize', 9);
    %temp_title = title(sprintf('Sequence-level evidence \n Baseline Meta'),'fontsize',11);
    %temp_title.Interpreter = 'none';
    
    
    
end


%% Part 5B, Right: Baseline & meta
baselineMeta_highMeta_summary = allPart_summary_mulFOV.part5_summary_mulFOV.baselineMeta_highMeta_summary;
baselineMeta_lowMeta_summary = allPart_summary_mulFOV.part5_summary_mulFOV.baselineMeta_lowMeta_summary;

temp_1 = baselineMeta_highMeta_summary';
temp_2 = baselineMeta_lowMeta_summary';

[~,temp_p,~,~] = ttest(temp_1,temp_2);
if if_plot == 1 && if_plot_part5 == 1
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[10 450 355*0.5*0.9 200*1.8*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    xtl = ["High-meta"; "Low-meta"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    if if_monkey_D0_Z1 == 0
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.53;%0.56,0.4
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.09;%0.56,0.4
    elseif if_monkey_D0_Z1 == 1
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.245;%0.56,0.4
    end
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%25
    
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Baseline meta', 'FontSize', 9);
    %title(sprintf('Meta-memory'),'fontsize',10);
end


%% Part 5B, Right: Baseline & mismatch
baselineMeta_highMatch_summary = allPart_summary_mulFOV.part5_summary_mulFOV.baselineMeta_highMatch_summary;
baselineMeta_lowMatch_summary = allPart_summary_mulFOV.part5_summary_mulFOV.baselineMeta_lowMatch_summary;
baselineMeta_overMismatch_summary = allPart_summary_mulFOV.part5_summary_mulFOV.baselineMeta_overMismatch_summary;
baselineMeta_underMismatch_summary = allPart_summary_mulFOV.part5_summary_mulFOV.baselineMeta_underMismatch_summary;

temp_1 = baselineMeta_highMatch_summary';
temp_2 = baselineMeta_underMismatch_summary';

temp_3 = baselineMeta_lowMatch_summary';
temp_4 = baselineMeta_overMismatch_summary';

[~,temp_p12,~,~] = ttest(temp_1,temp_2);
[~,temp_p34,~,~] = ttest(temp_3,temp_4);

if if_plot == 1 && if_plot_part5 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5*1.8 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 450 355*0.5*1.8*0.9*0.9 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*0.5*1.8*0.9*0.9*0.9 200*1.8*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp_y_min = min([temp_1;temp_2;temp_3;temp_4]);
    temp_y_max = max([temp_1;temp_2;temp_3;temp_4]);
    
    if if_plot_violinplot0_pairline1 == 0
        temp_data = [temp_1';temp_2';temp_3';temp_4'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        g4 = repmat({'D'},length(temp_4),1);
        
        temp_label = [g1;g2;g3;g4];
        
        temptemp_color2 = [...
            color_choiceMemoryHigh;
            color_choiceOffloadHigh;
            color_choiceOffloadLow;
            color_choiceMemoryLow
            ];
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.5,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
        h(1).ViolinPlot.FaceAlpha = 0.3;
        h(2).ViolinPlot.FaceAlpha = 0.3;
        h(3).ViolinPlot.FaceAlpha = 0.3;
        h(4).ViolinPlot.FaceAlpha = 0.3;
        
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([4-0.35 4+0.35],[1 1]*mean(temp_4),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        
        temptemp_color2 = [...
            color_choiceMemoryHigh;
            color_choiceOffloadHigh;
            color_choiceOffloadLow;
            color_choiceMemoryLow
            ];
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);%[1 1 1]*0.4
            hold on
            
            plot([3 4],[temp_3(tempi) temp_4(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',temptemp_color2(1,:),... %[1 1 1]*0.05
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7,'MarkerEdgeColor',[1 1 1]*0.05);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',temptemp_color2(2,:),...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7,'MarkerEdgeColor',[1 1 1]*0.05);
        hold on
        scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',temptemp_color2(3,:),...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7,'MarkerEdgeColor',[1 1 1]*0.05);
        hold on
        scatter(4*ones(1,length(temp_4)),temp_4,15,'filled','MarkerFaceColor',temptemp_color2(4,:),...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7,'MarkerEdgeColor',[1 1 1]*0.05);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p12 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p12 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p12 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p34 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p34 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p34 < 0.05
        tempTxt = sprintf('*');
    end
    text(3.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 4.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xtl = ["High-match"; "Under-mismatch";"Low-match"; "Over-mismatch"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    if if_monkey_D0_Z1 == 0
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.58;%0.56,0.4
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;%0.56,0.4
    elseif if_monkey_D0_Z1 == 1
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.25;%0.56,0.4
    end
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%25
    
    set(gca,'xticklabel','');
    
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Baseline meta', 'FontSize', 9);
    %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
    %temp_title.Interpreter = 'none';
end



%% Part 5C: Linear regression: baseline meta-memory, memory precision, meta-memory
r2_metaRegress_caseA_resampled_mean_summary = allPart_summary_mulFOV.part5_summary_mulFOV.r2_metaRegress_caseA_resampled_mean_summary;
r2_metaRegress_caseB_resampled_mean_summary = allPart_summary_mulFOV.part5_summary_mulFOV.r2_metaRegress_caseB_resampled_mean_summary;
r2_metaRegress_caseC_resampled_mean_summary = allPart_summary_mulFOV.part5_summary_mulFOV.r2_metaRegress_caseC_resampled_mean_summary;

% temp_1 = r2_metaRegress_caseA_resampled_mean_summary;
% temp_2 = r2_metaRegress_caseB_resampled_mean_summary;
% temp_3 = r2_metaRegress_caseC_resampled_mean_summary;


temp_1 = r2_metaRegress_caseB_resampled_mean_summary;
temp_2 = r2_metaRegress_caseA_resampled_mean_summary;
temp_3 = r2_metaRegress_caseC_resampled_mean_summary;

% if true
if false
    %temptempBoolIndex = r_mean_correct_summary>0.33;%0.33
    temptempBoolIndex = temp_2 > 0.0125;
    temp_1 = temp_1(temptempBoolIndex);
    temp_2 = temp_2(temptempBoolIndex);
    temp_3 = temp_3(temptempBoolIndex);
    
    %[~,temp_p12,~,~] = ttest(temp_1,temp_2,'Tail','left');
    
end

% [~,temp_p12,~,~] = ttest(temp_1,temp_2);
[~,temp_p13,~,~] = ttest(temp_1,temp_3);
[~,temp_p23,~,~] = ttest(temp_2,temp_3);

[~,temp_p12,~,~] = ttest(temp_1,temp_2,'Tail','left');


if false
    
    temptempBoolIndex = r_mean_correct_summary>0.38;%0.33
    [~,temp_p12,~,~] = ttest(temp_1(temptempBoolIndex),temp_2(temptempBoolIndex),'Tail','left');
    
    
    temptempBoolIndex = r_mean_correct_summary>-1;%0.33
    x = [];
    y = [];
    for tempi=1:3
        if tempi == 1
            temp_x = temp_1(temptempBoolIndex);
        elseif tempi == 2
            temp_x = temp_2(temptempBoolIndex);
        elseif tempi == 3
            temp_x = temp_3(temptempBoolIndex);
        end
        x = [x;tempi*ones(size(temp_x,2),1)];    %#ok<*AGROW>
        y = [y;temp_x'];
    end
    %temp_mdl = fitglm(x,y,'linear');
    temp_mdl = fitglm(x,y,'linear','intercept',false);
    temp_p123 = temp_mdl.Coefficients.pValue(end)
    %temp_p123 = coefTest(temp_mdl)
    
end

if if_plot == 1 && if_plot_part5 == 1
    
    fig = figure('Name','locDistri','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80*1.2 336*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50+0 50+0 240*0.80*1.2*0.9 336*1.11*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_y_min = min([temp_1 temp_2 temp_3]);
    temp_y_max = max([temp_1 temp_2 temp_3]);
    
    
    if if_plot_violinplot0_pairline1 == 0
        temp_data = [temp_1';temp_2';temp_3'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        
        temp_label = [g1;g2;g3];
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 3, 1);
        
        violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'}]);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2 3],[temp_1(tempi) temp_2(tempi)  temp_3(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p12 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p12 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p12 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.06,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    tempTxt = sprintf('');
    if temp_p13 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p13 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p13 < 0.05
        tempTxt = sprintf('*');
    end
    text(2,temp_y_max+(temp_y_max-temp_y_min)*0.15,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 2.9],temp_y_max+(temp_y_max-temp_y_min)*0.13*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    tempTxt = sprintf('');
    if temp_p23 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p23 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p23 < 0.05
        tempTxt = sprintf('*');
    end
    text(2.5,temp_y_max+(temp_y_max-temp_y_min)*0.06,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([2.1 2.9],temp_y_max+(temp_y_max-temp_y_min)*0.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    
    set(gca,'linewidth',1.5)
    
    xticks([1 2 3]);
    xlim([0 4]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    %xtl = ["Precision", "Baseline", "Both"];
    %xtl = ["Both", "Precision", "Baseline"];
    xtl = ["Baseline", "Strength", "Both"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.17;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('r2', 'FontSize', 9);
    temp_title = title(sprintf('Meta-memory regression'),'fontsize',9);
    temp_title.Interpreter = 'none';
end


%% (New) Part 5C: Linear regression: baseline meta-memory, memory precision, meta-memory
if if_plot_part5 == 1
    beta_precision_metaRegress_caseC_resampled_mean_summary = allPart_summary_mulFOV.part5_summary_mulFOV.beta_precision_metaRegress_caseC_resampled_mean_summary;
    beta_baseline_metaRegress_caseC_resampled_mean_summary = allPart_summary_mulFOV.part5_summary_mulFOV.beta_baseline_metaRegress_caseC_resampled_mean_summary;
    beta_interaction_metaRegress_caseC_resampled_mean_summary = allPart_summary_mulFOV.part5_summary_mulFOV.beta_interaction_metaRegress_caseC_resampled_mean_summary;
    
    temp_r2 = mean(r2_metaRegress_caseC_resampled_mean_summary);
    
    
    temp_1 = beta_baseline_metaRegress_caseC_resampled_mean_summary;
    temp_2 = beta_precision_metaRegress_caseC_resampled_mean_summary;
    temp_3 = beta_interaction_metaRegress_caseC_resampled_mean_summary;
    
    temp_chanceLevel = 0;
    
    % [~,temp_p,~,~] = ttest(temp_1,temp_2);
    [~,temp_p,~,~] = ttest(temp_1,temp_2,'tail','left');
    
    if true
        temp_1 = beta_baseline_metaRegress_caseC_resampled_mean_summary;
        temp_2 = beta_precision_metaRegress_caseC_resampled_mean_summary;
        
        %temptempBoolIndex = r2_metaRegress_caseC_resampled_mean_summary > 0.04;
        %temptempBoolIndex = r_mean_correct_summary>0.24;%0.33
        %temptempBoolIndex = r2_metaRegress_caseA_resampled_mean_summary > 0.0125;
        
        %temp_threshold = 0.10;
        %temp_threshold = 0.08;
        %temptempBoolIndex = (temp_1>temp_threshold) & (temp_2>temp_threshold);
        
        %     temp_threshold = 0;%20,
        %     %temptempBoolIndex = (temp_1>prctile(temp_1,temp_threshold)) & (temp_2>prctile(temp_2,temp_threshold));
        %
        %     temp_threshold2 = prctile([temp_1,temp_2],temp_threshold);
        %     temptempBoolIndex = (temp_1>temp_threshold2) & (temp_2>temp_threshold2);
        
        %temp_1
        %temp_2
        
        %     temptempBoolIndex1 = (~isoutlier(temp_1)) & (~isoutlier(temp_2));
        %     temptempBoolIndex1 = (~isoutlier(temp_1,'mean')) & (~isoutlier(temp_2,'mean'));
        %     temptempBoolIndex1 = (~isoutlier(temp_1,'quartiles')) & (~isoutlier(temp_2,'quartiles'));%good for all (choice+noChoice) trials
        %     temptempBoolIndex1 = (~isoutlier(temp_1,'grubbs')) & (~isoutlier(temp_2,'grubbs'));%
        temptempBoolIndex1 = (~isoutlier(temp_1,'gesd')) & (~isoutlier(temp_2,'gesd'));%good for choice trials
        %temptempBoolIndex1
        
        temp_threshold = 10;%20,10
        %temptempBoolIndex = (temp_1>prctile(temp_1,temp_threshold)) & (temp_2>prctile(temp_2,temp_threshold));
        
        temp_threshold2 = prctile([temp_1,temp_2],temp_threshold);
        temptempBoolIndex2 = (temp_1>temp_threshold2) & (temp_2>temp_threshold2);
        
        
        
        %temp_factor = 1.5;%Z(1.5),D(2.11)
        %temptempBoolIndex1 = ~isoutlier(temp_1,'mean','ThresholdFactor',temp_factor);
        %temptempBoolIndex2 = ~isoutlier(temp_2,'mean','ThresholdFactor',temp_factor);
        
        %temp_factor = 0.1;%Z(0.98),D(0.1)        
        %temptempBoolIndex1 = ~isoutlier(temp_1,'gesd','ThresholdFactor',temp_factor);
        %temptempBoolIndex2 = ~isoutlier(temp_2,'gesd','ThresholdFactor',temp_factor);
        
        %temp_factor = 1;%Z(0.2,0.43),D(0.2,0.43)        
        %temptempBoolIndex1 = ~isoutlier(temp_1,'quartiles','ThresholdFactor',temp_factor);
        %temptempBoolIndex2 = ~isoutlier(temp_2,'quartiles','ThresholdFactor',temp_factor);
        
        
        temp_factor = 1.5;%Z(1,1.5),D(1.5)   
        temp_factor2 = 1.5;%Z(1.5,1.7),D(any)                
        temptempBoolIndex1 = ~isoutlier(temp_1,'quartiles','ThresholdFactor',temp_factor);
        temptempBoolIndex2 = ~isoutlier(temp_2,'quartiles','ThresholdFactor',temp_factor);
        
        %temp_factor = 3.8;%Z(1,1.5,2,3,3.5,3.8),D(1.5,4,3.8)   
        %temp_factor2 = 1.5;%Z(1.5,1.7),D(any)                
        %temptempBoolIndex1 = ~isoutlier(temp_1,'median','ThresholdFactor',temp_factor);
        %temptempBoolIndex2 = ~isoutlier(temp_2,'median','ThresholdFactor',temp_factor);
        
        temptempBoolIndex = temptempBoolIndex1 & temptempBoolIndex2;
                     
        temptempBoolIndex = temptempBoolIndex & ~isoutlier(allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_theoretical_summary,'median','ThresholdFactor',temp_factor2);%1.5,1.7
               
        %temptempBoolIndex = temptempBoolIndex1;
        %     temptempBoolIndex = temptempBoolIndex2;
        
        %temptempBoolIndex = true(size(temptempBoolIndex1));
        
        %sum(temptempBoolIndex)
        temp_1 = temp_1(temptempBoolIndex);
        temp_2 = temp_2(temptempBoolIndex);
        [~,temp_p,~,~] = ttest(temp_1,temp_2,'Tail','left');
        %[~,temp_p,~,~] = ttest2(temp_1,temp_2,'Tail','left');
        %     [~,temp_p,~,~] = ttest(temp_1,temp_2);
        %    temp_p
        
        
        temp_r2 = mean(r2_metaRegress_caseC_resampled_mean_summary(temptempBoolIndex));
        
        temp_3 = temp_3(temptempBoolIndex);
    end
    
    [~,temp_p1,~,~] = ttest(temp_1,temp_chanceLevel,'tail','right');
    [~,temp_p2,~,~] = ttest(temp_2,temp_chanceLevel,'tail','right');
    [~,temp_p3,~,~] = ttest(temp_3,temp_chanceLevel,'tail','right');
    
end

if if_plot == 1 && if_plot_part5 == 1
    
    fig = figure('Name','locDistri','NumberTitle','off');
    set(gcf,'Position',[350+0 50+0 240*0.80*1.2*0.9*0.9*1.02 336*1.11*0.9*0.91*0.95*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_y_min = min([temp_1 temp_2 temp_3]);
    temp_y_max = max([temp_1 temp_2 temp_3]);
    
    temp_y_min = -0.1;
    temp_y_max = 0.4232;
    
    [~,temp_y1_max] = bounds(temp_1);
    [~,temp_y2_max] = bounds(temp_2);
    [~,temp_y3_max] = bounds(temp_3);
    
    
    if if_plot_violinplot0_pairline1 == 0
        temp_data = [temp_1';temp_2';temp_3'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'B'},length(temp_3),1);
        
        temp_label = [g1;g2;g3];
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 3, 1);
        
        violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'}]);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    
    plot([0 4],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.25,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.0 2.0],temp_y_max+(temp_y_max-temp_y_min)*0.22*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    
    
    tempTxt = sprintf('');
    if temp_p1 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p1 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p1 < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y1_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p2 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p2 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p2 < 0.05
        tempTxt = sprintf('*');
    end
    text(2,temp_y2_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p3 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p3 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p3 < 0.05
        tempTxt = sprintf('*');
    end
    text(3,temp_y3_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    xticks([1 2 3]);
    
    xtl = ["Baseline", "WM strength", "Interaction"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.17;
    xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.155;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
    
    set(gca,'xticklabel','');
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0.25 3.75]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.35]);%0.3
    set(gca, 'FontSize', 8) %14
    set(gca,'box','off');% 取消右、上边框
    
    
    
    ylabel('Beta', 'FontSize', 9);
    temp_title = title(sprintf('Meta-WM regression'),'fontsize',9);
    temp_title.Interpreter = 'none';
    subtitle(sprintf('r2 = %.3f',temp_r2),'fontsize',8);%9
end



meta_trialLevel_baseline_summary = allPart_summary_mulFOV.part5_summary_mulFOV.meta_trialLevel_baseline_summary;
memoryPrecision_trialLevel_summary = allPart_summary_mulFOV.part5_summary_mulFOV.memoryPrecision_trialLevel_summary;
meta_trialLevel_summary = allPart_summary_mulFOV.part5_summary_mulFOV.meta_trialLevel_summary;

% meta_trialLevel_baseline_summary = allPart_summary_mulFOV.part5_summary_mulFOV.meta_trialLevel_baseline_cell_summary{8};
% memoryPrecision_trialLevel_summary = allPart_summary_mulFOV.part5_summary_mulFOV.memoryPrecision_trialLevel_cell_summary{8};
% meta_trialLevel_summary = allPart_summary_mulFOV.part5_summary_mulFOV.meta_trialLevel_cell_summary{8};

tempNANBoolIndex1 = isnan(meta_trialLevel_baseline_summary);
tempNANBoolIndex2 = isnan(memoryPrecision_trialLevel_summary);
tempNANBoolIndex3 = isnan(meta_trialLevel_summary);

tempNANBoolIndex123 = tempNANBoolIndex1 | tempNANBoolIndex2 | tempNANBoolIndex3;
tempNONNANBoolIndex123 = ~tempNANBoolIndex123;

tempNONNAN_meta_trialLevel_baseline = meta_trialLevel_baseline_summary(tempNONNANBoolIndex123);
tempNONNAN_memoryPrecision_trialLevel = memoryPrecision_trialLevel_summary(tempNONNANBoolIndex123);
tempNONNAN_meta_trialLevel = meta_trialLevel_summary(tempNONNANBoolIndex123);

[temptemp_r_12,temptemp_p_12] = corr(tempNONNAN_meta_trialLevel_baseline,tempNONNAN_meta_trialLevel);
[temptemp_r_13,temptemp_p_13] = corr(tempNONNAN_memoryPrecision_trialLevel,tempNONNAN_meta_trialLevel);

% Case A
x = tempNONNAN_memoryPrecision_trialLevel;
y = tempNONNAN_meta_trialLevel;
temp_mdl_caseA = fitglm(x,y,'linear');
r2_caseA = temp_mdl_caseA.Rsquared.Adjusted;
beta0_caseA = temp_mdl_caseA.Coefficients.Estimate(1);
beta1_caseA = temp_mdl_caseA.Coefficients.Estimate(2);

% Case B
x = tempNONNAN_meta_trialLevel_baseline;
y = tempNONNAN_meta_trialLevel;
temp_mdl_caseB = fitglm(x,y,'linear');
r2_caseB = temp_mdl_caseB.Rsquared.Adjusted;
beta0_caseB = temp_mdl_caseB.Coefficients.Estimate(1);
beta1_caseB = temp_mdl_caseB.Coefficients.Estimate(2);

% Case C
x = [tempNONNAN_memoryPrecision_trialLevel,tempNONNAN_meta_trialLevel_baseline];
y = tempNONNAN_meta_trialLevel;
temp_mdl_caseC = fitglm(x,y,'linear');
r2_caseC = temp_mdl_caseC.Rsquared.Adjusted;
beta0_caseC = temp_mdl_caseC.Coefficients.Estimate(1);
beta1_caseC= temp_mdl_caseC.Coefficients.Estimate(2);
beta2_caseC= temp_mdl_caseC.Coefficients.Estimate(3);





%% Part 5: Beta of WM strength in linear regression in all trials, memory trials, and offload trials
beta_precision_metaRegress_caseC_resampled_mean_summary = allPart_summary_mulFOV.part5_summary_mulFOV.beta_precision_metaRegress_caseC_resampled_mean_summary;
beta_precision_metaRegress_caseE2_CM_resampled_mean_summary = allPart_summary_mulFOV.part5_summary_mulFOV.beta_precision_metaRegress_caseE2_CM_resampled_mean_summary;
beta_precision_metaRegress_caseF_CF_resampled_mean_summary = allPart_summary_mulFOV.part5_summary_mulFOV.beta_precision_metaRegress_caseF_CF_resampled_mean_summary;

temp_1 = beta_precision_metaRegress_caseC_resampled_mean_summary;
temp_2 = beta_precision_metaRegress_caseE2_CM_resampled_mean_summary;
temp_3 = beta_precision_metaRegress_caseF_CF_resampled_mean_summary;

% temptempBoolIndex = ~isoutlier(allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_summary,'percentiles',[20 100]);
temptempBoolIndex = ~isoutlier(allPart_summary_mulFOV.part2_summary_mulFOV.r_mean_correct_theoretical_summary,'median','ThresholdFactor',1.5);
temp_1 = temp_1(temptempBoolIndex);
temp_2 = temp_2(temptempBoolIndex);
temp_3 = temp_3(temptempBoolIndex);

r_chanceLevel = 0;

[~,temp_p1] = ttest(temp_1,r_chanceLevel,'tail','right');
[~,temp_p2] = ttest(temp_2,r_chanceLevel,'tail','right');
[~,temp_p3] = ttest(temp_3,r_chanceLevel,'tail','right');
% [~,temp_p1] = ttest(temp_1,r_chanceLevel);
% [~,temp_p2] = ttest(temp_2,r_chanceLevel);
% [~,temp_p3] = ttest(temp_3,r_chanceLevel);

if if_plot == 1 && if_plot_part5 == 1
    fig = figure('Name','locDistri','NumberTitle','off');
    set(gcf,'Position',[50+550 50+0 125*1.4 336*1.11*0.85]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(3,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    temp_y_min = -0.15;%0,-0.15
    %temp_y_max = max([temp_1 temp_2 temp_3 r_chanceLevel]);
    temp_y_max = 0.35;%0.2
    
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
    
    plot([0 4],r_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p1 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p1 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p1 < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p2 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p2 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p2 < 0.05
        tempTxt = sprintf('*');
    end
    text(2,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p3 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p3 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p3 < 0.05
        tempTxt = sprintf('*');
    end
    text(3,temp_y_min-(temp_y_max-temp_y_min)*0.05,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 4]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.2 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    %set(gca,'ydir','reverse');
    
    xtl = ["All","Memory","Offload"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %     if if_monkey_D0_Z1 == 0
    %         %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.35;
    %         xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.32;
    %     elseif if_monkey_D0_Z1 == 1
    %         xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min+(temp_y_max-temp_y_min)*2.34;
    %     end
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.3;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Beta', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('WM strength contribution\nto Meta-WM'),'fontsize',8.5);
    temp_title.Interpreter = 'none';
    
end

%% Part 5: Baseline contribution in mismatch regression
r2_mismatchRegress_caseD2B_baseline_resampled_mean_summary = allPart_summary_mulFOV.part5_summary_mulFOV.r2_mismatchRegress_caseD2B_baseline_resampled_mean_summary;

temp_1 = r2_mismatchRegress_caseD2B_baseline_resampled_mean_summary;
temp_chanceLevel = 0;

[~,temp_p] = ttest(temp_1,temp_chanceLevel);

if if_plot == 1 && if_plot_part5 == 1
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.9*0.9 336*1.11*0.5*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    
    %temp_y_min = min([temp_1 temp_chanceLevel]);
    %temp_y_max = max(temp_1);
    temp_y_min = 0;
    temp_y_max = 0.08;%0.1    
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    
    hold on
    
    
    plot([0 2],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.1:1);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Explained variance', 'FontSize', 9);
    temp_title = title(sprintf('Baseline contribution\nto mismatch regression'),'fontsize',9);
    temp_title.Interpreter = 'none';
    
end


%% Part 5: Beta of WM strength in linear regression in all trials, memory trials, and offload trials
memoryPrecision_CF_baseline_median_summary = allPart_summary_mulFOV.part5_summary_mulFOV.memoryPrecision_CF_baseline_median_summary;
memoryPrecision_CME_baseline_median_summary = allPart_summary_mulFOV.part5_summary_mulFOV.memoryPrecision_CME_baseline_median_summary;
meta_CF_baseline_median_summary = allPart_summary_mulFOV.part5_summary_mulFOV.meta_CF_baseline_median_summary;
meta_CME_baseline_median_summary = allPart_summary_mulFOV.part5_summary_mulFOV.meta_CME_baseline_median_summary;

temp_1 = memoryPrecision_CF_baseline_median_summary;
temp_2 = memoryPrecision_CME_baseline_median_summary;
temp_3 = meta_CF_baseline_median_summary;
temp_4 = meta_CME_baseline_median_summary;



[~,temp_p12] = ttest2(temp_1,temp_2,'tail','left');
[~,temp_p34] = ttest2(temp_3,temp_4,'tail','left');
% [~,temp_p12] = ttest2(temp_1,temp_2);
% [~,temp_p34] = ttest2(temp_3,temp_4);

temp_2minus1 = temp_2 - temp_1;
temp_4minus3 = temp_4 - temp_3;

[~,temp_p_diff] = ttest(temp_2minus1,temp_4minus3,'tail','left');


if if_plot == 1 && if_plot_part5 == 1   
    fig = figure('Name','asd','NumberTitle','off');
    set(gcf,'Position',[10 100 391.4*0.7 319.6]);
    t = tiledlayout(2,2,'TileSpacing','Loose','Padding','Loose');
    
    temp_y_min = 0.075;%0.075
    temp_y_max = 0.6;%0.6

    
    nexttile
    
    temp_data = [temp_1';temp_2'];
    
    g1 = repmat({'A'},length(temp_1),1);
    g2 = repmat({'B'},length(temp_2),1);
    
    temp_label = [g1;g2];
    
    temptemp_color1 = [1 1 1]*0.5;
    temptemp_color2 = repmat(temptemp_color1, 2, 1);
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',[{'A'};{'B'}]);
    h(1).ViolinPlot.FaceAlpha = 0.1;
    h(2).ViolinPlot.FaceAlpha = 0.1;
        
    
    tempTxt = sprintf('');
    if temp_p12 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p12 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p12 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.04*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on    
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0 3]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.2 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    %set(gca,'ydir','reverse');
    
    xtl = ["Offload","Low-strength mismatch"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;

    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.3;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('WM strength', 'FontSize', 9, 'FontWeight', 'normal');
    %temp_title = title(sprintf('WM strength contribution\nto Meta-WM'),'fontsize',8.5);
    %temp_title.Interpreter = 'none';
    
    
    
    
    nexttile
    
        
    temp_data = [temp_3';temp_4'];
    
    g1 = repmat({'A'},length(temp_3),1);
    g2 = repmat({'B'},length(temp_4),1);
    
    temp_label = [g1;g2];
    
    temptemp_color1 = [1 1 1]*0.5;
    temptemp_color2 = repmat(temptemp_color1, 2, 1);
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',[{'A'};{'B'}]);
    h(1).ViolinPlot.FaceAlpha = 0.1;
    h(2).ViolinPlot.FaceAlpha = 0.1;
        
    
    tempTxt = sprintf('');
    if temp_p34 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p34 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p34 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.07,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    plot([1.1 1.9],temp_y_max+(temp_y_max-temp_y_min)*0.04*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    hold on
    
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0 3]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.2 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.3
    set(gca, 'FontSize', 8) %14
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    %set(gca,'ydir','reverse');
    
    xtl = ["Offload","Low-strength mismatch"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;

    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.3;
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Baseline meta', 'FontSize', 9, 'FontWeight', 'normal');
    %temp_title = title(sprintf('WM strength contribution\nto Meta-WM'),'fontsize',8.5);
    %temp_title.Interpreter = 'none';    
end


%% Part 5E: Trial history weight
AUROC_optimal_trialHistory_summary = allPart_summary_mulFOV.part5_summary_mulFOV.AUROC_optimal_trialHistory_summary;
AUROC_chanceLevel = 0.5;

[~,temp_p] = ttest(AUROC_optimal_trialHistory_summary,AUROC_chanceLevel);

if if_plot == 1 && if_plot_part5 == 1
    fig = figure('Name','','NumberTitle','off');
    set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.9*0.9 336*1.11*0.5*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_1 = AUROC_optimal_trialHistory_summary;
    
    temp_y_min = min([temp_1 AUROC_chanceLevel]);
    temp_y_max = max(temp_1);
    
    
    temp_data = temp_1';
    temp_label = repmat({'A'},length(temp_1),1);
    
    temptemp_color1 = [1 1 1]*0.5;
    
    h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
        'GroupOrder',{'A'});
    h(1).ViolinPlot.FaceAlpha = 0.1;
    
    hold on
    
    
    plot([0 2],AUROC_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    
    xlim([0 2]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    set(gca,'YTick',0:0.1:1);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = [""]; %#ok<*NBRAK>
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('AUROC', 'FontSize', 9);
    temp_title = title(sprintf('History weight fitting'),'fontsize',9);
    temp_title.Interpreter = 'none';
    
end


normStd_optimal_trialHistory_summary = allPart_summary_mulFOV.part5_summary_mulFOV.normStd_optimal_trialHistory_summary;
if sum(normStd_optimal_trialHistory_summary) > 0
    if if_plot == 1 && if_plot_part5 == 1
        fig = figure('Name','','NumberTitle','off');
        set(gcf,'Position',[50+0 50+0 240*0.80*1.1*0.9*0.9 336*1.11*0.5*0.9*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_1 = normStd_optimal_trialHistory_summary;
        
        temp_y_min = min(temp_1);
        temp_y_max = max(temp_1);
        
        
        temp_data = temp_1';
        temp_label = repmat({'A'},length(temp_1),1);
        
        temptemp_color1 = [1 1 1]*0.5;
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',{'A'});
        h(1).ViolinPlot.FaceAlpha = 0.1;
        
        hold on
        
        set(gca,'linewidth',1.5)
        
        xlim([0 2]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
        set(gca, 'FontSize', 8) %14
        %set(gca,'YTick',0:0.1:1,'FontSize', 12);%12
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = [""]; %#ok<*NBRAK>
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
        
        set(gca,'xticklabel','');
        
        
        %ylabel(sprintf('Std of \n norm distribution'), 'FontSize', 11, 'FontWeight', 'bold');
        ylabel(sprintf('Std of norm'), 'FontSize', 9);
        temp_title = title(sprintf('History weight'),'fontsize',9);
        temp_title.Interpreter = 'none';
        
    end
end

%% Part 5, history weight
weight_trialHistory_summary = allPart_summary_mulFOV.part5_summary_mulFOV.weight_trialHistory_summary;
AUROC_optimal_trialHistory_summary = allPart_summary_mulFOV.part5_summary_mulFOV.AUROC_optimal_trialHistory_summary;


if false
    weight_trialHistory_summary_D = allPart_summary_mulFOV.part5_summary_mulFOV.weight_trialHistory_summary;
    AUROC_optimal_trialHistory_summary_D = allPart_summary_mulFOV.part5_summary_mulFOV.AUROC_optimal_trialHistory_summary;
    
    temp1 = weight_trialHistory_summary_D;
    temp1 = rescale(weight_trialHistory_summary_D,0,1);
    weight_trialHistory_summary_D_n11n = temp1./sum(temp1,2);
    
    %weight_trialHistory_summary = weight_trialHistory_summary_D_n11n;
end

if false
    weight_trialHistory_summary_Z = allPart_summary_mulFOV.part5_summary_mulFOV.weight_trialHistory_summary;
    AUROC_optimal_trialHistory_summary_Z = allPart_summary_mulFOV.part5_summary_mulFOV.AUROC_optimal_trialHistory_summary;
    
    temp1 = weight_trialHistory_summary_Z;
    temp1 = rescale(weight_trialHistory_summary_Z,0,1);
    weight_trialHistory_summary_Z_n11n = temp1./sum(temp1,2);
    
    weight_trialHistory_summary = weight_trialHistory_summary_Z_n11n;
end

if false
    %weight_trialHistory_summary = [weight_trialHistory_summary_D;weight_trialHistory_summary_Z];
    weight_trialHistory_summary = [weight_trialHistory_summary_D_n11n;weight_trialHistory_summary_Z_n11n];
    AUROC_optimal_trialHistory_summary = [AUROC_optimal_trialHistory_summary_D,AUROC_optimal_trialHistory_summary_Z];
    
    %temptempBoolIndex = AUROC_optimal_trialHistory_summary >= prctile(AUROC_optimal_trialHistory_summary,0);
    
    
    [M,I] = sort(AUROC_optimal_trialHistory_summary_Z,'descend');
    temp1 = false(size(AUROC_optimal_trialHistory_summary_Z));
    temp1(I(1:length(AUROC_optimal_trialHistory_summary_D))) = true;
    temp2 = find(temp1==false) + length(AUROC_optimal_trialHistory_summary_D);
    temptempBoolIndex = true(size(AUROC_optimal_trialHistory_summary));
    temptempBoolIndex(temp2) = false;
    
    
    weight_trialHistory_summary = weight_trialHistory_summary(temptempBoolIndex,:);
    AUROC_optimal_trialHistory_summary = AUROC_optimal_trialHistory_summary(temptempBoolIndex);
end


if true
    weight_trialHistory_summary = weight_trialHistory_summary(:,1:10);
    weight_trialHistory_summary = weight_trialHistory_summary./sum(weight_trialHistory_summary,2);
end

if true
    temp_1 = weight_trialHistory_summary;
    %    temptempBoolIndex1 = ~isoutlier(temp_1(:,1),'median',1)';
    temptempBoolIndex1 = ~isoutlier(temp_1(:,1),'mean',1)';
    %    temptempBoolIndex1 = ~isoutlier(temp_1(:,1),'grubbs',1)';
    %    temptempBoolIndex1 = ~isoutlier(temp_1(:,1),'gesd',1)';
    %    temptempBoolIndex1 = ~isoutlier(temp_1(:,1),'quartiles',1)';
    
    temptempBoolIndex1;
    
    weight_trialHistory_summary = weight_trialHistory_summary(temptempBoolIndex1,:);
    AUROC_optimal_trialHistory_summary = AUROC_optimal_trialHistory_summary(temptempBoolIndex1);
end

weight_trialHistory_summary_mean = mean(weight_trialHistory_summary,1);
% weight_trialHistory_summary_mean = median(weight_trialHistory_summary,1);
weight_trialHistory_summary_sem = std(weight_trialHistory_summary,1,1)./sqrt(size(weight_trialHistory_summary,1));

% AUROC_optimal_trialHistory_summary_median = median(AUROC_optimal_trialHistory_summary);
AUROC_optimal_trialHistory_summary_median = mean(AUROC_optimal_trialHistory_summary);

temp_chanceLevel = 1/size(weight_trialHistory_summary,2);%1/20


% [~,temp_p] = ttest(weight_trialHistory_summary,0);
[~,temp_p] = ttest(weight_trialHistory_summary,temp_chanceLevel,'tail','right');

if false
    [~,temp_p] = ttest(weight_trialHistory_summary(:,1),0);
end

if false
    temp_prc = prctile(weight_trialHistory_summary,5);
    temptempBoolIndex = temp_prc > 0;
end

p_threshold = 0.05;%0.01,0.05

if if_plot == 1 && if_plot_part5 == 1
    %% Plot history weight
    tempPlot_historyIndicator = 1;%big(0);small(1)
    fig = figure('Name','asd','NumberTitle','off');
    %set(gcf,'Position',[1110 400 240*0.8*0.85*1.85*0.7*1.5*0.95 143*1.02*0.88*1.2*1.2]);
    %set(gcf,'Position',[1110 400 240*0.8*0.85*1.85*0.7*1.5*0.95 143*1.02*0.88*1.2*1.2*1.05]);
    %set(gcf,'Position',[1110 400 301.2*0.48*1.04 194.1*0.92*1.04]);
    %set(gcf,'Position',[1110 400 301.2*0.48*1.04*1.15 194.1*0.92*1.04]);
    %set(gcf,'Position',[1110 400 301.2*0.48*1.04*1.15*0.95*0.97 194.1*0.92*1.04*0.82*1.06]);
    %set(gcf,'Position',[1110 400 301.2*0.48*1.04*1.15*0.95*0.97*1.06 194.1*0.92*1.04*0.82*1.06]);
    %set(gcf,'Position',[1110 400 301.2*0.48*1.04*1.15*0.95*0.97*1.06*0.9 194.1*0.92*1.04*0.82*1.06]);
    if tempPlot_historyIndicator == 0
    set(gcf,'Position',[1110 400 301.2*0.48*1.04*1.15*0.95*0.97*1.06 194.1*0.92*1.04*0.82*1.06]);
    elseif tempPlot_historyIndicator == 1
    set(gcf,'Position',[1110 400 301.2*0.48*1.04*1.15*0.95*0.97*1.06 194.1*0.92*1.04*0.82*1.06*0.75*0.95]);
    end
    
    t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
    
    
    nexttile
    
    x = 1:length(weight_trialHistory_summary_mean);
    y = weight_trialHistory_summary_mean;
    
    %     for tempi=1:size(weight_trialHistory_summary,1)
    %         scatter(x,weight_trialHistory_summary(tempi,:),3,'filled','MarkerFaceColor',[1 1 1]*0.5,...
    %             'MarkerFaceAlpha',0.5);
    %         %scatter(x,weight_trialHistory_summary(tempi,:),2.5,'filled','MarkerFaceColor',[1 1 1]*0.5);
    %         hold on
    %     end
    
    
    %plot(x,y,'color',[0.25 0.25 0.25],'linewidth',2);%2
    plot(x,y,'color',[0 0 0],'linewidth',2);%2
    hold on
    
    y_sem = weight_trialHistory_summary_sem;
    %patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0.25 0.25 0.25],'FaceAlpha',0.05,'EdgeColor',[0.62 0.62 0.62]);
    %hold on
    
    
    for tempi=1:size(weight_trialHistory_summary,1)
        scatter(x,weight_trialHistory_summary(tempi,:),9,'filled','MarkerFaceColor',[1 1 1]*0.35,...
            'MarkerFaceAlpha',0.3);
        hold on
    end
    
    
    %[y_min,y_max] = bounds([y-y_sem y+y_sem]);%old
    [y_min,y_max] = bounds(weight_trialHistory_summary(:));
    
    y_min = min([y_min,0]);
    
    %y_max = y_max +
    
    scatter(x(temp_p<p_threshold),y_max+(y_max-y_min)*0.15,8,[0 0 0],'*');
    
    
    %plot([x(1) x(end)],temp_chanceLevel*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1);
    %hold on
    
    %xticks([1:8,10,x(end)]);
    
    %xticks([1:2:7,10,x(end)]);
    %xticks([1:7,10,x(end)]);
    xticks([1:7,10]);
    
    xtickangle(0);
    
    set(gca,'linewidth',1.5)
    xlim([x(1)-1 x(end)+1]);
    ylim([y_min-(y_max-y_min)*0.1 y_max+(y_max-y_min)*0.20]);%0.1
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    
    %yticks([0 0.1]);
    %yticks([0 temp_chanceLevel floor(y_max*10)/10]);
    if if_monkey_D0_Z1 == 0
        yticks([0 floor(y_max*10)/10]);
    elseif if_monkey_D0_Z1 == 1
        yticks([0 0.1]);
    end
    
    xlabel('Trial history', 'FontSize', 9, 'FontWeight', 'normal');
    ylabel('Weight', 'FontSize', 9, 'FontWeight', 'normal');
    title(sprintf('History weight'),'FontSize',9,'FontWeight','normal');
    
    if tempPlot_historyIndicator == 0
    %text(14,y_max+(y_max-y_min)*0.10,...
    text(7,y_max+(y_max-y_min)*0.00,...
        sprintf('AUROC = %.3f',AUROC_optimal_trialHistory_summary_median),'FontSize',7,...
        'HorizontalAlignment','center');
    elseif tempPlot_historyIndicator == 1
    text(7,y_max+(y_max-y_min)*0.00,...
        sprintf('AUROC = %.3f',AUROC_optimal_trialHistory_summary_median),'FontSize',7,...
        'HorizontalAlignment','center');        
    end
    
    %subtitle(sprintf('AUROC = %.3f',AUROC_optimal_trialHistory_summary_median),'FontSize',7,'FontWeight','normal');
    
end

%% Part 5, history weight
% weight_trialHistory_summary_raw = allPart_summary_mulFOV.part5_summary_mulFOV.weight_trialHistory_summary;
%
% % weight_trialHistory_summary = allPart_summary_mulFOV.part5_summary_mulFOV.weight_trialHistory_summary(:,1:20);
% weight_trialHistory_summary = allPart_summary_mulFOV.part5_summary_mulFOV.weight_trialHistory_summary;

weight_trialHistory_summary_predict = nan(size(weight_trialHistory_summary));

for tempi=1:size(weight_trialHistory_summary,1)
    temp1 = weight_trialHistory_summary(tempi,:);
    
    x = 1:length(temp1);
    y = temp1;
    % cftool(x,y);
    
    [xData, yData] = prepareCurveData(x,y);
    ft = fittype( 'exp1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0 0];
    [fitresult, gof] = fit(xData,yData,ft,opts);
    
    temp_r2 = gof.rsquare;
    
    y_predict = feval(fitresult,x)';
    
    if false
        close all
        plot(x,y);
        hold on
        plot(x,y_predict);
        hold on
    end
    
    weight_trialHistory_summary_predict(tempi,:) = y_predict;
    
end

weight_trialHistory_summary_predict;
if false
    close all
    plot(x,mean(weight_trialHistory_summary,1));
    hold on
    plot(x,mean(weight_trialHistory_summary_predict,1));
    hold on
end

p_weight_trialHistory_summary_predict = nan(1,size(weight_trialHistory_summary,2));
for tempi=1:size(weight_trialHistory_summary,2)
    [~,p_weight_trialHistory_summary_predict(tempi)] = ttest(weight_trialHistory_summary_predict(:,tempi),0);
    %[~,p_weight_trialHistory_summary_predict(tempi)] = ttest2(weight_trialHistory_summary_predict(:,tempi),weight_trialHistory_summary_raw(:));
end

%% Part 5F, Left: Trial history of baseline
historyRewardMean_baselineHigh_summary = allPart_summary_mulFOV.part5_summary_mulFOV.historyRewardMean_baselineHigh_summary;
historyRewardMean_baselineLow_summary = allPart_summary_mulFOV.part5_summary_mulFOV.historyRewardMean_baselineLow_summary;

temp_1 = historyRewardMean_baselineHigh_summary';
temp_2 = historyRewardMean_baselineLow_summary';

[~,temp_p,~,~] = ttest(temp_1,temp_2);

if false
    [~,temp_p,~,~] = ttest(temp_1,temp_2,'tail','right');
end

if if_plot == 1 && if_plot_part5 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*0.5*0.9 200*1.8*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    %xtl = ["High-meta"; "Low-meta"];
    xtl = ["High-baseline"; "Low-baseline"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %     if if_monkey_D0_Z1 == 0
    %         %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.315;%0.56,0.4
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.31;%0.56,0.4
    %     elseif if_monkey_D0_Z1 == 1
    %         %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.25;%0.56,0.4
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.305;%0.56,0.4
    %     end
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.20;
    
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%25
    
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Weighted reward', 'FontSize', 9);
    %title(sprintf('Baseline'),'fontsize',10);
end

%% Part 5F, Right: Trial history of meta
historyRewardMean_metaHigh_summary = allPart_summary_mulFOV.part5_summary_mulFOV.historyRewardMean_metaHigh_summary;
historyRewardMean_metaLow_summary = allPart_summary_mulFOV.part5_summary_mulFOV.historyRewardMean_metaLow_summary;

temp_1 = historyRewardMean_metaHigh_summary';
temp_2 = historyRewardMean_metaLow_summary';

[~,temp_p,~,~] = ttest(temp_1,temp_2);
% [~,temp_p,~,~] = ttest(temp_1,temp_2,'tail','right');

if if_plot == 1 && if_plot_part5 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*0.5*0.9 200*1.8*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    
    
    temp_y_min = min([temp_1;temp_2]);
    temp_y_max = max([temp_1;temp_2]);
    
    if if_plot_violinplot0_pairline1 == 0
        violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 2.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xticks([1 2]);
    xtl = ["High-meta"; "Low-meta"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %     if if_monkey_D0_Z1 == 0
    %         %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.38;%0.56,0.4
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.335;%0.56,0.4
    %     elseif if_monkey_D0_Z1 == 1
    %         %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.28;%0.56,0.4
    %         xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.30;%0.56,0.4
    %     end
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.20;
    
    %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%25
    
    set(gca,'xticklabel','');
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Weighted reward', 'FontSize', 9);
    %title(sprintf('Meta-memory'),'fontsize',10);
end


%% Part 5F, Right: Trial history of mismatch
historyRewardMean_OverMismatch_summary = allPart_summary_mulFOV.part5_summary_mulFOV.historyRewardMean_OverMismatch_summary;
historyRewardMean_HighMatch_summary = allPart_summary_mulFOV.part5_summary_mulFOV.historyRewardMean_HighMatch_summary;
historyRewardMean_LowMatch_summary = allPart_summary_mulFOV.part5_summary_mulFOV.historyRewardMean_LowMatch_summary;
historyRewardMean_UnderMismatch_summary = allPart_summary_mulFOV.part5_summary_mulFOV.historyRewardMean_UnderMismatch_summary;

temp_1 = historyRewardMean_HighMatch_summary';
temp_2 = historyRewardMean_UnderMismatch_summary';

temp_3 = historyRewardMean_LowMatch_summary';
temp_4 = historyRewardMean_OverMismatch_summary';

[~,temp_p12,~,~] = ttest(temp_1,temp_2);
[~,temp_p34,~,~] = ttest(temp_3,temp_4);

if if_plot == 1 && if_plot_part5 == 1
    fig = figure('Name','','NumberTitle','off');
    %set(gcf,'Position',[10 450 355*0.5*1.8 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 450 355*0.5*1.8*0.9*0.9 200*1.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 355*0.5*1.8*0.9*0.9*0.9 200*1.8*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    
    nexttile
    
    temp_y_min = min([temp_1;temp_2;temp_3;temp_4]);
    temp_y_max = max([temp_1;temp_2;temp_3;temp_4]);
    
    if if_plot_violinplot0_pairline1 == 0
        temp_data = [temp_1';temp_2';temp_3';temp_4'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        g4 = repmat({'D'},length(temp_4),1);
        
        temp_label = [g1;g2;g3;g4];
        
        temptemp_color2 = [...
            color_choiceMemoryHigh;
            color_choiceOffloadHigh;
            color_choiceOffloadLow;
            color_choiceMemoryLow
            ];
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.5,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
        h(1).ViolinPlot.FaceAlpha = 0.3;
        h(2).ViolinPlot.FaceAlpha = 0.3;
        h(3).ViolinPlot.FaceAlpha = 0.3;
        h(4).ViolinPlot.FaceAlpha = 0.3;
        
        
    elseif if_plot_violinplot0_pairline1 == 1
        plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([3-0.35 3+0.35],[1 1]*mean(temp_3),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        plot([4-0.35 4+0.35],[1 1]*mean(temp_4),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
        hold on
        
        
        temptemp_color2 = [...
            color_choiceMemoryHigh;
            color_choiceOffloadHigh;
            color_choiceOffloadLow;
            color_choiceMemoryLow
            ];
        
        for tempi=1:length(temp_1)
            plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);%[1 1 1]*0.4
            hold on
            
            plot([3 4],[temp_3(tempi) temp_4(tempi)],'Color',[1 1 1]*0.4);
            hold on
        end
        
        scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',temptemp_color2(1,:),... %[1 1 1]*0.05
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7,'MarkerEdgeColor',[1 1 1]*0.05);
        hold on
        scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',temptemp_color2(2,:),...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7,'MarkerEdgeColor',[1 1 1]*0.05);
        hold on
        scatter(3*ones(1,length(temp_3)),temp_3,15,'filled','MarkerFaceColor',temptemp_color2(3,:),...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7,'MarkerEdgeColor',[1 1 1]*0.05);
        hold on
        scatter(4*ones(1,length(temp_4)),temp_4,15,'filled','MarkerFaceColor',temptemp_color2(4,:),...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7,'MarkerEdgeColor',[1 1 1]*0.05);
        hold on
    end
    
    tempTxt = sprintf('');
    if temp_p12 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p12 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p12 < 0.05
        tempTxt = sprintf('*');
    end
    text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    tempTxt = sprintf('');
    if temp_p34 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p34 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p34 < 0.05
        tempTxt = sprintf('*');
    end
    text(3.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    set(gca,'linewidth',1.5)
    xlim([0.15 4.85])
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
    set(gca, 'FontSize', 8)
    
    xtl = ["High-match"; "Under-mismatch";"Low-match"; "Over-mismatch"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    if if_monkey_D0_Z1 == 0
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.445;%0.56,0.4
    elseif if_monkey_D0_Z1 == 1
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.325;%0.56,0.4
    end
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',8);%25
    set(gca,'xticklabel','');
    
    
    set(gca,'box','off');% 取消右、上边框
    ylabel('Weighted reward', 'FontSize', 9);
    %temp_title = title(sprintf(''),'FontSize',12,'FontWeight','bold');
    %temp_title.Interpreter = 'none';
end



% % overMismatch-->highMatch-->lowMatch-->underMismatch
% x = [1*ones(size(temp_1,2),1);2*ones(size(temp_2,2),1);3*ones(size(temp_3,2),1);4*ones(size(temp_4,2),1)];
% y = [temp_1,temp_2,temp_3,temp_4]';
%
% temp_mdl = fitglm(x,y,'linear');
% p_linear = temp_mdl.Coefficients.pValue(2);
%
%
% % [~,temp_p,~,~] = ttest(temp_1,temp_4,'Tail','Right');
% [~,temp_p,~,~] = ttest(temp_1,temp_4); %#ok<*ASGLU>
% % [~,temp_p,~,~] = ttest2(temp_1,temp_4); %#ok<*ASGLU>
%
% [~,temp_p,~,~] = ttest([temp_1 temp_2],[temp_3 temp_4]);
%
%
% color_choiceMemoryLow = [1,133,113]/255; %[5,113,176]/255
% color_choiceMemoryHigh = [128,205,193]/255; %[146,197,222]/255
% color_choiceOffloadLow = [223,194,125]/255; %[244,165,130]/255
% color_choiceOffloadHigh = [166,97,26]/255; %[202,0,32]/255
%
% if if_plot == 1 && if_plot_part5 == 1
%     fig = figure('Name','locDistri','NumberTitle','off');
%     set(gcf,'Position',[50+0 50+0 240*0.80 336*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
%
%     t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
%     nexttile
%
%     temp_y_min = min([temp_1 temp_2 temp_3 temp_4]);
%     temp_y_max = max([temp_1 temp_2 temp_3 temp_4]);
%
%     temp_data = [temp_1';temp_2';temp_3';temp_4'];
%
%     g1 = repmat({'A'},length(temp_1),1);
%     g2 = repmat({'B'},length(temp_2),1);
%     g3 = repmat({'C'},length(temp_3),1);
%     g4 = repmat({'D'},length(temp_4),1);
%
%     temp_label = [g1;g2;g3;g4];
%
%     temptemp_color2 = ...
%         [color_choiceMemoryLow;
%         color_choiceMemoryHigh;
%         color_choiceOffloadLow;
%         color_choiceOffloadHigh];
%
%     h = violinplot(temp_data,temp_label,'ViolinAlpha',0.5,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
%         'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
%     h(1).ViolinPlot.FaceAlpha = 0.3;
%     h(2).ViolinPlot.FaceAlpha = 0.3;
%     h(3).ViolinPlot.FaceAlpha = 0.3;
%     h(4).ViolinPlot.FaceAlpha = 0.3;
%
%
%     set(gca,'linewidth',1.5)
%
%     xlim([0 5]);
%     ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.2]);%0.3
%     set(gca, 'FontSize', 12) %14
%     %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
%     set(gca,'box','off');% 取消右、上边框
%
%
%     %xtl = ["OverMismatch", "HighMatch", "LowMatch", "UnderMismatch"];
%     xtl = ["", "", "", ""];
%     xt=get(gca,'XTick');
%     yt=get(gca,'YTick');
%     xtext_xp=xt;
%     xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
%     text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
%
%     set(gca,'xticklabel','');
%
%
%     ylabel('Mean reward', 'FontSize', 12, 'FontWeight', 'bold');
%     temp_title = title(sprintf('Delay1'),'fontsize',12);
%     temp_title.Interpreter = 'none';
%
% end



%% Part 6: Neuron tuning structure
if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
    tempStr = 'trial-level';
elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
    tempStr = 'seq-level';
elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
    tempStr = 'mix-level';
end

%% Single neuron tuning property of memoryPrecision, choiceMemory, baseline (trial-level tuning)
if false
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
        %set(gcf,'Position',[510 450 240*3 240*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[510 450 720 290]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[510 450 240*3*0.95 240*1.05*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[510 450 240*3*0.67 240*1.05*0.67*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        
        %t.Title.String = sprintf('Multi-FOV (%s)\n Single neuron tuning property, n=%d',tempStr,length(beta_memoryPrecision_summary));
        t.Title.String = sprintf('Multi-FOV (%s)\n Regression coefficients of neurons, n=%d',tempStr,length(beta_memoryPrecision_summary));
        t.Title.FontSize = 9;%11
        t.Title.Interpreter = 'none';
        
        
        tempSize = 4;%10
        
        nexttile
        
        %x = allPart_summary_mulFOV.part6_summary_mulFOV.beta_memoryPrecision_summary;
        %y = allPart_summary_mulFOV.part6_summary_mulFOV.beta_choiceMemory_summary;
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_memoryPrecision_summary;
            y = beta_choiceMemory_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_seqPrecision_summary;
            y = beta_gProb_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_seqPrecision_summary;
            y = beta_gProb_summary;
            %y = beta_choiceMemory_summary;
        end
        
        [temp_r,temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        xy_min = -1*xy_max;
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        tempTxt = 'p>0.05';
        if temp_p < 0.05
            tempTxt = 'p<0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p<0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p<0.001';
        end
        %     text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.18,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
        %     text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.05,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        
        
        %axis equal
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Memory strength', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r=%.3f, %s',temp_r,tempTxt), 'FontSize', 9, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
        
        nexttile
        
        %x = allPart_summary_mulFOV.part6_summary_mulFOV.beta_memoryPrecision_summary;
        %y = allPart_summary_mulFOV.part6_summary_mulFOV.beta_choiceMemory_baseline_summary;
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_memoryPrecision_summary;
            y = beta_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_seqPrecision_summary;
            y = beta_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_seqPrecision_summary;
            y = beta_choiceMemory_baseline_summary;
        end
        
        [temp_r,temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        xy_min = -1*xy_max;
        
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        
        tempTxt = 'p>0.05';
        if temp_p < 0.05
            tempTxt = 'p<0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p<0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p<0.001';
        end
        %     if if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        %         text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.98,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
        %         text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.85,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        %     elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        %         text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.18,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
        %         text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.05,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        %     end
        
        
        %axis equal
        
        set(gca,'linewidth',1.5)
        
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Memory strength', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Baseline', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r=%.3f, %s',temp_r,tempTxt), 'FontSize', 9, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
        nexttile
        
        %x = allPart_summary_mulFOV.part6_summary_mulFOV.beta_choiceMemory_summary;
        %y = allPart_summary_mulFOV.part6_summary_mulFOV.beta_choiceMemory_baseline_summary;
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_choiceMemory_summary;
            y = beta_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_gProb_summary;
            y = beta_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_gProb_summary;
            %x = beta_choiceMemory_summary;
            y = beta_choiceMemory_baseline_summary;
        end
        
        [temp_r,temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        xy_min = -1*xy_max;
        
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        
        tempTxt = 'p>0.05';
        if temp_p < 0.05
            tempTxt = 'p<0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p<0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p<0.001';
        end
        %     if if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        %         text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.98,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
        %         text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.85,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        %     elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        %         text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.18,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
        %         text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.05,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        %     end
        
        
        
        %axis equal
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Baseline', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r=%.3f, %s',temp_r,tempTxt), 'FontSize', 9, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
    end
end




%% Single neuron tuning property of memory (location), choiceMemory, baseline
if true
    if if_plot == 1 && if_plot_part6 == 1
        
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            temp_1 = beta_6loc_summary_peak;
            temp_2 = beta_choiceMemory_summary;
            temp_3 = beta_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            temp_1 = beta_6loc_summary_peak;
            temp_2 = beta_gProb_summary;
            temp_3 = beta_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            temp_1 = beta_6loc_summary_peak;
            temp_2 = beta_gProb_summary;
            temp_3 = beta_choiceMemory_baseline_summary;
        end
        
        [temp_123_min_raw,temp_123_max_raw] = bounds([temp_1;temp_2;temp_3]);
        
        temp_123_max = max([abs(temp_123_min_raw),abs(temp_123_max_raw)]);
        temp_123_min = -1*temp_123_max;
        
        
        fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
        %set(gcf,'Position',[510 450 240*3*0.67 240*1.05*0.67*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[510 450 240*3*0.67*0.81 240*1.05*0.67*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        t = tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
        
        t.Title.String = sprintf('Multi-FOV (%s)\n Regression coefficients of neurons, n=%d',tempStr,length(beta_memoryPrecision_summary(selectiveNeuronBoolIndex)));
        t.Title.FontSize = 9;%11
        t.Title.Interpreter = 'none';
        
        
        tempSize = 4;%10
        
        nexttile
        
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_6loc_summary_peak;
            y = beta_choiceMemory_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_6loc_summary_peak;
            y = beta_gProb_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_6loc_summary_peak;
            y = beta_gProb_summary;
        end
        x = x(selectiveNeuronBoolIndex);
        y = y(selectiveNeuronBoolIndex);
        
        [temp_r,temp_p] = corr(x,y);
        
        %[x_min,x_max] = bounds(x);
        %[y_min,y_max] = bounds(y);
        %xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        %xy_min = -1*xy_max;
        
        xy_max = temp_123_max;
        xy_min = temp_123_min;
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        
        set(gca,'linewidth',1.5)
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('WM', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Meta-WM', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r = %.3f, %s',temp_r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
        
        nexttile
        
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_6loc_summary_peak;
            y = beta_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_6loc_summary_peak;
            y = beta_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_6loc_summary_peak;
            y = beta_choiceMemory_baseline_summary;
        end
        x = x(selectiveNeuronBoolIndex);
        y = y(selectiveNeuronBoolIndex);
        
        
        [temp_r,temp_p] = corr(x,y);
        
        %[x_min,x_max] = bounds(x);
        %[y_min,y_max] = bounds(y);
        %xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        %xy_min = -1*xy_max;
        
        xy_max = temp_123_max;
        xy_min = temp_123_min;
        
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        
        set(gca,'linewidth',1.5)
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('WM', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Baseline', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r = %.3f, %s',temp_r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
        nexttile
        
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_choiceMemory_summary;
            y = beta_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_gProb_summary;
            y = beta_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_gProb_summary;
            y = beta_choiceMemory_baseline_summary;
        end
        x = x(selectiveNeuronBoolIndex);
        y = y(selectiveNeuronBoolIndex);
        
        
        [temp_r,temp_p] = corr(x,y);
        
        %[x_min,x_max] = bounds(x);
        %[y_min,y_max] = bounds(y);
        %xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        %xy_min = -1*xy_max;
        
        xy_max = temp_123_max;
        xy_min = temp_123_min;
        
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        
        set(gca,'linewidth',1.5)
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-WM', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Baseline', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r = %.3f, %s',temp_r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
    end
end

%% Single neuron tuning property of memory (location), choiceMemory, baseline
if true
    if if_plot == 1 && if_plot_part6 == 1
        
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            temp_1 = r2_6loc_summary;
            temp_2 = r2_choiceMemory_summary;
            temp_3 = r2_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            temp_1 = r2_6loc_summary;
            temp_2 = r2_gProb_summary;
            temp_3 = r2_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            temp_1 = r2_6loc_summary;
            temp_2 = r2_gProb_summary;
            temp_3 = r2_choiceMemory_baseline_summary;
        end
        
        [temp_123_min_raw,temp_123_max_raw] = bounds([temp_1;temp_2;temp_3]);
        
        temp_123_max = max([temp_123_min_raw,temp_123_max_raw]);
        %temp_123_min = -1*temp_123_max;
        temp_123_min = min([temp_123_min_raw,temp_123_max_raw]);
        
        
        fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
        %set(gcf,'Position',[510 450 240*3*0.67 240*1.05*0.67*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[510 450 240*3*0.67*0.81 240*1.05*0.67*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        t = tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
        
        t.Title.String = sprintf('Multi-FOV (%s)\n Regression explained variance of neurons, n=%d',tempStr,length(beta_memoryPrecision_summary(selectiveNeuronBoolIndex)));
        t.Title.FontSize = 9;%11
        t.Title.Interpreter = 'none';
        
        
        tempSize = 4;%10
        
        nexttile
        
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = r2_6loc_summary;
            y = r2_choiceMemory_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = r2_6loc_summary;
            y = r2_gProb_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = r2_6loc_summary;
            y = r2_gProb_summary;
        end
        x = x(selectiveNeuronBoolIndex);
        y = y(selectiveNeuronBoolIndex);
        
        [temp_r,temp_p] = corr(x,y);
        
        %[x_min,x_max] = bounds(x);
        %[y_min,y_max] = bounds(y);
        %xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        %xy_min = -1*xy_max;
        
        xy_max = temp_123_max;
        xy_min = temp_123_min;
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        %         plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        %         plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        
        set(gca,'linewidth',1.5)
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('WM', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Meta-WM', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r = %.3f, %s',temp_r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
        
        nexttile
        
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = r2_6loc_summary;
            y = r2_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = r2_6loc_summary;
            y = r2_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = r2_6loc_summary;
            y = r2_choiceMemory_baseline_summary;
        end
        x = x(selectiveNeuronBoolIndex);
        y = y(selectiveNeuronBoolIndex);
        
        
        [temp_r,temp_p] = corr(x,y);
        
        %[x_min,x_max] = bounds(x);
        %[y_min,y_max] = bounds(y);
        %xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        %xy_min = -1*xy_max;
        
        xy_max = temp_123_max;
        xy_min = temp_123_min;
        
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        %         plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        %         plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        
        set(gca,'linewidth',1.5)
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('WM', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Baseline', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r = %.3f, %s',temp_r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
        nexttile
        
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = r2_choiceMemory_summary;
            y = r2_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = r2_gProb_summary;
            y = r2_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = r2_gProb_summary;
            y = r2_choiceMemory_baseline_summary;
        end
        x = x(selectiveNeuronBoolIndex);
        y = y(selectiveNeuronBoolIndex);
        
        
        [temp_r,temp_p] = corr(x,y);
        
        %[x_min,x_max] = bounds(x);
        %[y_min,y_max] = bounds(y);
        %xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        %xy_min = -1*xy_max;
        
        xy_max = temp_123_max;
        xy_min = temp_123_min;
        
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        %         plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        %         plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        
        set(gca,'linewidth',1.5)
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-WM', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Baseline', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r = %.3f, %s',temp_r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
    end
end


%% Single neuron tuning property of memory (location), choiceMemory, baseline
if true
    if if_plot == 1 && if_plot_part6 == 1
        
        temp_1 = beta_6loc_summary_norm;
        temp_2 = beta_gProb_summary_norm;
        temp_3 = beta_choiceMemory_baseline_summary_norm;
        
        
        [temp_123_min_raw,temp_123_max_raw] = bounds([temp_1;temp_2;temp_3]);
        
        temp_123_max = max([temp_123_min_raw,temp_123_max_raw]);
        %temp_123_min = -1*temp_123_max;
        temp_123_min = min([temp_123_min_raw,temp_123_max_raw]);
        
        
        fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
        %set(gcf,'Position',[510 450 240*3*0.67 240*1.05*0.67*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[510 450 240*3*0.67*0.81 240*1.05*0.67*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        t = tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
        
        t.Title.String = sprintf('Multi-FOV (%s)\n Regression coefficients of neurons, n=%d',tempStr,length(beta_memoryPrecision_summary(selectiveNeuronBoolIndex)));
        t.Title.FontSize = 9;%11
        t.Title.Interpreter = 'none';
        
        
        tempSize = 4;%10
        
        nexttile
        
        x = beta_6loc_summary_norm;
        y = beta_gProb_summary_norm;
        
        x = x(selectiveNeuronBoolIndex);
        y = y(selectiveNeuronBoolIndex);
        
        [temp_r,temp_p] = corr(x,y);
        
        %[x_min,x_max] = bounds(x);
        %[y_min,y_max] = bounds(y);
        %xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        %xy_min = -1*xy_max;
        
        xy_max = temp_123_max;
        xy_min = temp_123_min;
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        %         plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        %         plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        
        set(gca,'linewidth',1.5)
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('WM', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Meta-WM', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r = %.3f, %s',temp_r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
        
        nexttile
        
        x = beta_6loc_summary_norm;
        y = beta_choiceMemory_baseline_summary_norm;
        
        x = x(selectiveNeuronBoolIndex);
        y = y(selectiveNeuronBoolIndex);
        
        
        [temp_r,temp_p] = corr(x,y);
        
        %[x_min,x_max] = bounds(x);
        %[y_min,y_max] = bounds(y);
        %xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        %xy_min = -1*xy_max;
        
        xy_max = temp_123_max;
        xy_min = temp_123_min;
        
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        %         plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        %         plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        
        set(gca,'linewidth',1.5)
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('WM', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Baseline', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r = %.3f, %s',temp_r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
        nexttile
        
        x = beta_gProb_summary_norm;
        y = beta_choiceMemory_baseline_summary_norm;
        
        x = x(selectiveNeuronBoolIndex);
        y = y(selectiveNeuronBoolIndex);
        
        
        [temp_r,temp_p] = corr(x,y);
        
        %[x_min,x_max] = bounds(x);
        %[y_min,y_max] = bounds(y);
        %xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        %xy_min = -1*xy_max;
        
        xy_max = temp_123_max;
        xy_min = temp_123_min;
        
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        %         plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        %         plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        %         hold on
        
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        
        set(gca,'linewidth',1.5)
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-WM', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Baseline', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('r = %.3f, %s',temp_r,tempTxt), 'FontSize', 8, 'FontWeight', 'normal');
        temp_title.Interpreter = 'none';
        
    end
end


%% ePRIRS
if if_ePRIRS == 1
    % x = beta_memoryPrecision_summary;
    % y = beta_choiceMemory_summary;
    % z = beta_choiceMemory_baseline_summary;
    %     if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
    %         x = beta_memoryPrecision_summary;
    %         y = beta_choiceMemory_summary;
    %         z = beta_choiceMemory_baseline_summary;
    %     elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
    %         x = beta_seqPrecision_summary;
    %         y = beta_gProb_summary;
    %         z = beta_gProb_baseline_summary;
    %     elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
    %         x = beta_seqPrecision_summary;
    %         y = beta_gProb_summary;
    %         z = beta_choiceMemory_baseline_summary;
    %     end
    
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        x = beta_6loc_summary_peak;
        y = beta_choiceMemory_summary;
        z = beta_choiceMemory_baseline_summary;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        x = beta_6loc_summary_peak;
        y = beta_gProb_summary;
        z = beta_gProb_baseline_summary;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        x = beta_6loc_summary_peak;
        y = beta_gProb_summary;
        z = beta_choiceMemory_baseline_summary;
    end
    
    temp_vectors = [x,y,z];
    
    if if_compute_part6 == 1
        %if false
        [clusteriness, p, dists] = pairsClusterTest_elliptical(temp_vectors);
    end
    
    temp_angles = dists.data;
    temp_angles_bootstrap = dists.bootstrap;
    % temp_angles_bootstrap_median = median(temp_angles_bootstrap,1)';
    % temp_angles_bootstrap_median = median(temp_angles_bootstrap,2);
    % temp_angles_bootstrap_median = temp_angles_bootstrap(1,:)';
    temp_angles_bootstrap_median = reshape(temp_angles_bootstrap,[],1);
    
    %% Part 6: Neuron tuning structure, ePRIRS
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name','asd','NumberTitle','off');
        %set(gcf,'Position',[10 450 355/2*2 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        % set(gcf,'Position',[10 450 177*1.3 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        % set(gcf,'Position',[10 450 177*2 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 450 412-10 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 450 402*0.95 200*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 450 402*0.43*1.15*1.05*0.98 200*0.53*1.08*1.47*0.98]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
        
        nexttile
        
        x1 = temp_angles;
        %x1 = temp_angles_bootstrap_median;
        x2 = temp_angles_bootstrap_median;
        
        [~,temp_p] = ttest2(x1,x2);
        
        
        h_NumBins = 50;%10-->100-->400
        
        x = x1;
        h1 = histogram(x,'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.5);
        hold on
        %h1.NumBins = h_NumBins;
        h1.EdgeColor = [0.4660 0.6740 0.1880];%[0.25 0.75 0.25]
        
        x = x2;
        h2 = histogram(x,'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.5);
        hold on
        %h2.NumBins = h_NumBins;
        h2.EdgeColor = [0 0 0];
        
        
        x_min = min([x1;x2]);
        %x_max = max([x1;x2]);
        x_max = prctile([x1;x2],99);%97.5
        
        h1.BinLimits = [x_min,x_max];
        h1.NumBins = h_NumBins;
        
        h2.BinLimits = [x_min,x_max];
        h2.NumBins = h_NumBins;
        
        y1 = h1.Values;
        y2 = h2.Values;
        
        
        [y1_min,y1_max] = bounds(y1);
        [y2_min,y2_max] = bounds(y2);
        y_min = min([y1_min,y2_min]);
        y_max = max([y1_max,y2_max]);
        
        y_min = max([y_min,1]);
        
        if if_plot_additionalSmooth == 1
            n=1000;
            n=2^ceil(log2(n)); % round up n to the next power of 2;
            
            [pdf1,xmesh1,~] = ksdensity(x1,'NumPoints',n,'Function','pdf');
            
            %plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',1.5,'color',[0 0 0]);
            plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',2,'color',[0.4660 0.6740 0.1880]*0.7);
            hold on
            
            %[pdf2,xmesh2,~] = ksdensity(x2,'NumPoints',n,'Function','pdf');
            %plot(xmesh2,pdf2*sum(y2)*h2.BinWidth,':','LineWidth',1.5,'color',[0.5 0.5 0.5]*0.5);
            %hold on
        end
        
        le = legend('Observed','Bootstrap',...
            'Location','northeast','fontsize',7);
        le.ItemTokenSize = ones(1,2)*8;
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text((x_min+x_max)*0.3,y_max+(y_max-y_min)*0.00,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        
        %set(gca,'YScale','log');
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
        xlim([x_min-(x_max-x_min)*0.02 x_max+(x_max-x_min)*0.02]);
        ylim([y_min y_max+(y_max-y_min)*0.1]);%0.1
        %xticks([0 1]);
        set(gca, 'FontSize', 10)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Nearest-neighbour angle', 'FontSize', 10, 'FontWeight', 'normal');
        %ylabel('Frequency', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Pdf', 'FontSize', 10, 'FontWeight', 'normal');
        
        temp_title = title(sprintf('Distribution tests: ePAIRS'),'FontSize',10,'FontWeight','normal');
        temp_title.Interpreter = 'none';
        
    end
    
end

%% Tuning clustering
% if if_ePRIRS == 1
if false
    % x = beta_memoryPrecision_summary;
    % y = beta_choiceMemory_summary;
    % z = beta_choiceMemory_baseline_summary;
    
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        x = beta_memoryPrecision_summary;
        y = beta_choiceMemory_summary;
        z = beta_choiceMemory_baseline_summary;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        x = beta_seqPrecision_summary;
        y = beta_gProb_summary;
        z = beta_gProb_baseline_summary;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        x = beta_seqPrecision_summary;
        y = beta_gProb_summary;
        z = beta_choiceMemory_baseline_summary;
    end
    
    temp_vectors = [x,y,z];
    
    
    KList = 1:6;%1:6
    % criterion = 'CalinskiHarabasz';
    % criterion = 'DaviesBouldin';
    criterion = 'gap';
    % criterion = 'silhouette';
    
    if if_compute_part6 == 1
        temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','sqEuclidean');
        %temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','cosine');
        %temp_eva = evalclusters(temp_vectors,'gmdistribution',criterion,'KList',KList,'Distance','sqEuclidean');
    end
    
    CriterionValues_tuning = temp_eva.CriterionValues;
    clusterID_tuning = temp_eva.OptimalY;
    clusterNum_tuning = temp_eva.OptimalK;
    
    
    multi_rgbColor = ...
        [228,26,28;
        55,126,184;
        77,175,74;
        152,78,163;
        255,127,0;
        [255,255,51].*0.8]/255;
    
    
    
    %% Part 6: Neuron tuning structure, clustering
    if if_plot == 1 && if_plot_part6 == 1
        %% Plot clusterNum vs CriterionValues
        fig = figure('Name',' ','NumberTitle','off'); %#ok<*NASGU>
        %set(gcf,'Position',[600 520 264 260*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[600 520 308 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[600 520 308*0.95 200*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        h_line = [];
        
        x = KList;
        y = CriterionValues_tuning;
        h_line = [h_line plot(x,y,'LineWidth',1,'color',[0.25 0.25 0.25])];
        hold on
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        
        plot([1 1].*clusterNum_tuning,[y_min-(y_max-y_min)*0.8 CriterionValues_tuning(clusterNum_tuning)],':','color',[0.25 0.25 0.25],'linewidth',1);
        hold on
        
        plot([x_min-(x_max-x_min)*0.08 clusterNum_tuning],[1 1].*CriterionValues_tuning(clusterNum_tuning),':','color',[0.25 0.25 0.25],'linewidth',1);
        hold on
        
        xticks(KList);
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        set(gca, 'FontSize', 11)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Cluster number', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Clustering scores', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('Optimal cluster number'),'FontSize',11);
        temp_title.Interpreter = 'none';
        
        
        %% Single neuron tuning property of memoryPrecision and choiceMemory (trial-level tuning)
        fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
        set(gcf,'Position',[510 450 240*3 240*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        
        %t.Title.String = sprintf('Multi-FOV (trial-level)\n Single neuron tuning property, n=%d',length(x));
        t.Title.String = sprintf('Multi-FOV (%s)\n Single neuron tuning property, n=%d',tempStr,length(beta_memoryPrecision_summary));
        t.Title.FontSize = 11;
        t.Title.Interpreter = 'none';
        
        tempSize = 4;%10
        
        nexttile
        
        %x = beta_memoryPrecision_summary;
        %y = beta_choiceMemory_summary;
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_memoryPrecision_summary;
            y = beta_choiceMemory_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_seqPrecision_summary;
            y = beta_gProb_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_seqPrecision_summary;
            y = beta_gProb_summary;
        end
        
        [temp_r,temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        xy_min = -1*xy_max;
        
        %scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        %hold on
        
        for tempi=1:clusterNum_tuning
            temptempBoolIndex = clusterID_tuning==tempi;
            temptempColor = multi_rgbColor(tempi,:);
            scatter(x(temptempBoolIndex),y(temptempBoolIndex),tempSize,...
                'filled','MarkerFaceColor',temptempColor,'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
        
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        
        %axis equal
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Memory strength', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
        temp_title.Interpreter = 'none';
        
        
        nexttile
        %x = beta_memoryPrecision_summary;
        %y = beta_choiceMemory_baseline_summary;
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_memoryPrecision_summary;
            y = beta_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_seqPrecision_summary;
            y = beta_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_seqPrecision_summary;
            y = beta_choiceMemory_baseline_summary;
        end
        
        [temp_r,temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        xy_min = -1*xy_max;
        
        
        %scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        %hold on
        
        for tempi=1:clusterNum_tuning
            temptempBoolIndex = clusterID_tuning==tempi;
            temptempColor = multi_rgbColor(tempi,:);
            scatter(x(temptempBoolIndex),y(temptempBoolIndex),tempSize,...
                'filled','MarkerFaceColor',temptempColor,'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        %axis equal
        
        set(gca,'linewidth',1.5)
        
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Memory strength', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Baseline', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
        temp_title.Interpreter = 'none';
        
        nexttile
        
        %x = beta_choiceMemory_summary;
        %y = beta_choiceMemory_baseline_summary;
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_choiceMemory_summary;
            y = beta_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_gProb_summary;
            y = beta_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_gProb_summary;
            y = beta_choiceMemory_baseline_summary;
        end
        
        [temp_r,temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        xy_min = -1*xy_max;
        
        
        %scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        %hold on
        
        for tempi=1:clusterNum_tuning
            temptempBoolIndex = clusterID_tuning==tempi;
            temptempColor = multi_rgbColor(tempi,:);
            scatter(x(temptempBoolIndex),y(temptempBoolIndex),tempSize,...
                'filled','MarkerFaceColor',temptempColor,'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        %axis equal
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Baseline', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
        temp_title.Interpreter = 'none';
        
    end
end

%% Part 6E: Neuron spatial organization, clustering index
numFOV = size(allPart_summary_mulFOV.part6_summary_mulFOV.clusterIndex_memoryPrecision_disBin_mean_summary,2);

tempBoolIndex_validFOV = true(1,numFOV); %#ok<*PREALL>

temptempBoolIndex_A = r_mean_correct_summary>prctile(r_mean_correct_summary,25);
temptempBoolIndex_B = AUROC_meta_delay1_summary>prctile(AUROC_meta_delay1_summary,25);
temptempBoolIndex_C = AUROC_meta_baseline_summary>prctile(AUROC_meta_baseline_summary,25);

% tempBoolIndex_validFOV = temptempBoolIndex_A & temptempBoolIndex_B & temptempBoolIndex_C;
% tempBoolIndex_validFOV = temptempBoolIndex_A | temptempBoolIndex_B | temptempBoolIndex_C;

% tempBoolIndex_validFOV = temptempBoolIndex_A;
% tempBoolIndex_validFOV = temptempBoolIndex_C;


% sum(tempBoolIndex_validFOV)

clusterIndex_memoryPrecision_disBin_mean_summary = allPart_summary_mulFOV.part6_summary_mulFOV.clusterIndex_memoryPrecision_disBin_mean_summary(:,tempBoolIndex_validFOV);
clusterIndex_choiceMemory_disBin_mean_summary = allPart_summary_mulFOV.part6_summary_mulFOV.clusterIndex_choiceMemory_disBin_mean_summary(:,tempBoolIndex_validFOV);
clusterIndex_choiceMemory_baseline_disBin_mean_summary = allPart_summary_mulFOV.part6_summary_mulFOV.clusterIndex_choiceMemory_baseline_disBin_mean_summary(:,tempBoolIndex_validFOV);

temp_1 = clusterIndex_memoryPrecision_disBin_mean_summary;
temp_2 = clusterIndex_choiceMemory_disBin_mean_summary;
temp_3 = clusterIndex_choiceMemory_baseline_disBin_mean_summary;

if false
    a1 = isoutlier(temp_1,'median',2);
    %a1 = isoutlier(temp_1,'mean',2);
    %a1 = isoutlier(temp_1,'quartiles',2);
    %a1 = isoutlier(temp_1,'grubbs',2);
    %a1 = isoutlier(temp_1,'gesd',2);
    
    a2 = isoutlier(temp_2,'median',2);
    a3 = isoutlier(temp_3,'median',2);
    
end


temp_1_mean = mean(temp_1,2)';
temp_2_mean = mean(temp_2,2)';
temp_3_mean = mean(temp_3,2)';

temp_1_sem = std(temp_1,1,2)'./sqrt(size(temp_1,2));
temp_2_sem = std(temp_2,1,2)'./sqrt(size(temp_2,2));
temp_3_sem = std(temp_3,1,2)'./sqrt(size(temp_3,2));

% x = [];
% y = [];
% for tempi=1:size(temp_1,1)
%    x = [x; tempi*ones(size(temp_1,2),1)]; %#ok<*AGROW>
%    y = [y; temp_1(tempi,:)'];
% end
% temp_mdl = fitglm(x,y,'linear');
% p_linear_1 = temp_mdl.Coefficients.pValue(2);
%
% x = [];
% y = [];
% for tempi=1:size(temp_2,1)
%    x = [x; tempi*ones(size(temp_2,2),1)]; %#ok<*AGROW>
%    y = [y; temp_2(tempi,:)'];
% end
% temp_mdl = fitglm(x,y,'linear');
% p_linear_2 = temp_mdl.Coefficients.pValue(2);
%
% x = [];
% y = [];
% for tempi=1:size(temp_3,1)
%    x = [x; tempi*ones(size(temp_3,2),1)]; %#ok<*AGROW>
%    y = [y; temp_3(tempi,:)'];
% end
% temp_mdl = fitglm(x,y,'linear');
% p_linear_3 = temp_mdl.Coefficients.pValue(2);


x = 1:length(temp_1_mean);
y = temp_1_mean;
temp_mdl = fitglm(x,y,'linear');
p_linear_1_mean = temp_mdl.Coefficients.pValue(2);

% [r,p]=corr(x',y');

x = 1:length(temp_2_mean);
y = temp_2_mean;
temp_mdl = fitglm(x,y,'linear');
p_linear_2_mean = temp_mdl.Coefficients.pValue(2);

x = 1:length(temp_3_mean);
y = temp_3_mean;
temp_mdl = fitglm(x,y,'linear');
p_linear_3_mean = temp_mdl.Coefficients.pValue(2);


% prctile_high_memoryPrecision_shuffled_summary = allPart_summary_mulFOV.part6_summary_mulFOV.prctile_high_memoryPrecision_shuffled_summary;
% prctile_high_choiceMemory_shuffled_summary = allPart_summary_mulFOV.part6_summary_mulFOV.prctile_high_choiceMemory_shuffled_summary;
% prctile_high_choiceMemory_baseline_shuffled_summary = allPart_summary_mulFOV.part6_summary_mulFOV.prctile_high_choiceMemory_baseline_shuffled_summary;

clusterIndex_memoryPrecision_shuffled_meanB_summary = allPart_summary_mulFOV.part6_summary_mulFOV.clusterIndex_memoryPrecision_shuffled_meanB_summary(:,tempBoolIndex_validFOV);
clusterIndex_choiceMemory_shuffled_meanB_summary = allPart_summary_mulFOV.part6_summary_mulFOV.clusterIndex_choiceMemory_shuffled_meanB_summary(:,tempBoolIndex_validFOV);
clusterIndex_choiceMemory_baseline_shuffled_meanB_summary = allPart_summary_mulFOV.part6_summary_mulFOV.clusterIndex_choiceMemory_baseline_shuffled_meanB_summary(:,tempBoolIndex_validFOV);

% temp_1B = prctile_high_memoryPrecision_shuffled_summary;
% temp_2B = prctile_high_choiceMemory_shuffled_summary;
% temp_3B = prctile_high_choiceMemory_baseline_shuffled_summary;

temp_1B = clusterIndex_memoryPrecision_shuffled_meanB_summary;
temp_2B = clusterIndex_choiceMemory_shuffled_meanB_summary;
temp_3B = clusterIndex_choiceMemory_baseline_shuffled_meanB_summary;

if false
    a1 = isoutlier(temp_1B,'median',2);
    a1 = isoutlier(temp_1B,'mean',2);
    a1 = isoutlier(temp_1B,'quartiles',2);
    a1 = isoutlier(temp_1B,'grubbs',2);
    a1 = isoutlier(temp_1B,'gesd',2);
    
    temp_1B = clusterIndex_memoryPrecision_shuffled_meanB_summary;
    temp_2B = clusterIndex_choiceMemory_shuffled_meanB_summary;
    temp_3B = clusterIndex_choiceMemory_baseline_shuffled_meanB_summary;
    
end

temp_1B_mean = mean(temp_1B,2)';
temp_2B_mean = mean(temp_2B,2)';
temp_3B_mean = mean(temp_3B,2)';

temp_1B_sem = std(temp_1B,1,2)'./sqrt(size(temp_1B,2));
temp_2B_sem = std(temp_2B,1,2)'./sqrt(size(temp_2B,2));
temp_3B_sem = std(temp_3B,1,2)'./sqrt(size(temp_3B,2));

temptempBoolIndex1 = temp_1>temp_1B;
temptempBoolIndex2 = temp_2>temp_2B;
temptempBoolIndex3 = temp_3>temp_3B;

temptempBoolIndex1_mean = temp_1_mean>temp_1B_mean;
temptempBoolIndex2_mean = temp_2_mean>temp_2B_mean;
temptempBoolIndex3_mean = temp_3_mean>temp_3B_mean;

[~,temp_p1] = ttest(temp_1',temp_1B');
[~,temp_p2] = ttest(temp_2',temp_2B');
[~,temp_p3] = ttest(temp_3',temp_3B');

if true
    [~,temp_p1] = ttest(temp_1',temp_1B','tail','right');
    [~,temp_p2] = ttest(temp_2',temp_2B','tail','right');
    [~,temp_p3] = ttest(temp_3',temp_3B','tail','right');
end

dis_range = ...
    [0 200;
    150 350;
    300 500;
    450 650;
    600 800];

if if_plot == 1 && if_plot_part6 == 1
    fig = figure('Name',' ','NumberTitle','off');
    %set(gcf,'Position',[50 520 720 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50 520 720*0.9 200*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    
    
    %% memory precision
    nexttile
    
    h_line = [];
    
    x1 = 1:size(dis_range,1);
    y1 = temp_1_mean;
    y1_sem = temp_1_sem;
    h_line = [h_line plot(x1,y1,'LineWidth',1,'color',[0.25 0.25 0.25])];
    hold on
    
    x2 = 1:size(dis_range,1);
    y2 = temp_1B_mean;
    y2_sem = temp_1B_sem;
    h_line = [h_line plot(x2,y2,'LineWidth',1,'color',[0.75 0.75 0.75])];
    hold on
    
    x_min = min([x1,x2]);
    x_max = max([x1,x2]);
    y_min = min([y1-y1_sem,y2-y2_sem]);
    y_max = max([y1+y1_sem,y2+y2_sem]);
    
    scatter(x1(temp_p1<0.05),y_max+(y_max-y_min)*0.08,20,[0 0 0],'*');
    hold on
    
    le = legend(h_line,'Data','Shuffled','Location','southwest','fontsize',8);%10
    le.ItemTokenSize = ones(1,2)*10;
    legend('boxoff');
    
    
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
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Spatial size (um)','FontSize',9);
    ylabel('Clustering index', 'FontSize', 9);
    temp_title = title(sprintf('Memory'),'FontSize',9);
    temp_title.Interpreter = 'none';
    
    
    
    %% meta-memory
    nexttile
    
    h_line = [];
    
    x1 = 1:size(dis_range,1);
    y1 = temp_2_mean;
    y1_sem = temp_2_sem;
    h_line = [h_line plot(x1,y1,'LineWidth',1,'color',[0.25 0.25 0.25])];
    hold on
    
    x2 = 1:size(dis_range,1);
    y2 = temp_2B_mean;
    y2_sem = temp_2B_sem;
    h_line = [h_line plot(x2,y2,'LineWidth',1,'color',[0.75 0.75 0.75])];
    hold on
    
    x_min = min([x1,x2]);
    x_max = max([x1,x2]);
    y_min = min([y1-y1_sem,y2-y2_sem]);
    y_max = max([y1+y1_sem,y2+y2_sem]);
    
    scatter(x1(temp_p2<0.05),y_max+(y_max-y_min)*0.08,20,[0 0 0],'*');
    hold on
    
    %le = legend(h_line,'Data','Shuffled','Location','southwest','fontsize',8);%10
    %le.ItemTokenSize = ones(1,2)*10;
    
    
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
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Spatial size (um)','FontSize',9);
    ylabel('Clustering index', 'FontSize',9);
    temp_title = title(sprintf('Meta-memory'),'FontSize',9);
    temp_title.Interpreter = 'none';
    
    
    
    %% baseline meta-memory
    nexttile
    
    h_line = [];
    
    x1 = 1:size(dis_range,1);
    y1 = temp_3_mean;
    y1_sem = temp_3_sem;
    h_line = [h_line plot(x1,y1,'LineWidth',1,'color',[0.25 0.25 0.25])];
    hold on
    
    x2 = 1:size(dis_range,1);
    y2 = temp_3B_mean;
    y2_sem = temp_3B_sem;
    h_line = [h_line plot(x2,y2,'LineWidth',1,'color',[0.75 0.75 0.75])];
    hold on
    
    x_min = min([x1,x2]);
    x_max = max([x1,x2]);
    y_min = min([y1-y1_sem,y2-y2_sem]);
    y_max = max([y1+y1_sem,y2+y2_sem]);
    
    scatter(x1(temp_p3<0.05),y_max+(y_max-y_min)*0.08,20,[0 0 0],'*');
    hold on
    
    %le = legend(h_line,'Data','Shuffled','Location','southwest','fontsize',8);%10
    %le.ItemTokenSize = ones(1,2)*10;
    
    
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
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Spatial size (um)', 'FontSize',9);
    ylabel('Clustering index', 'FontSize',9);
    temp_title = title(sprintf('Baseline Meta'),'FontSize',9);
    temp_title.Interpreter = 'none';
    
end


%% Part 6F: Neuron spatial organization, centroid distance
tempBoolIndex_validFOV = true(1,numFOV); %#ok<*PREALL>

temptempBoolIndex_A = r_mean_correct_summary>prctile(r_mean_correct_summary,25);
temptempBoolIndex_B = AUROC_meta_delay1_summary>prctile(AUROC_meta_delay1_summary,25);
temptempBoolIndex_C = AUROC_meta_baseline_summary>prctile(AUROC_meta_baseline_summary,25);

% tempBoolIndex_validFOV = temptempBoolIndex_A & temptempBoolIndex_B & temptempBoolIndex_C;
% tempBoolIndex_validFOV = temptempBoolIndex_C;

% sum(tempBoolIndex_validFOV)

centriodDis_AB_all_summary = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_summary(tempBoolIndex_validFOV);
centriodDis_AC_all_summary = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_AC_all_summary(tempBoolIndex_validFOV);
centriodDis_BC_all_summary = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_BC_all_summary(tempBoolIndex_validFOV);

temp_1 = centriodDis_AB_all_summary;
temp_2 = centriodDis_AC_all_summary;
temp_3 = centriodDis_BC_all_summary;

% temp_1B = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_shuffled_prctileB_summary;
% temp_2B = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_AC_shuffled_prctileB_summary;
% temp_3B = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_BC_shuffled_prctileB_summary;
% temp_1B = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_shuffled_mean_summary;
% temp_2B = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_AC_shuffled_mean_summary;
% temp_3B = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_BC_shuffled_mean_summary;
temp_1B = (allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_shuffled_prctileA_summary(tempBoolIndex_validFOV)+...
    allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_shuffled_prctileB_summary(tempBoolIndex_validFOV))/2;
temp_2B = (allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_AC_shuffled_prctileA_summary(tempBoolIndex_validFOV)+...
    allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_AC_shuffled_prctileB_summary(tempBoolIndex_validFOV))/2;
temp_3B = (allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_BC_shuffled_prctileA_summary(tempBoolIndex_validFOV)+...
    allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_BC_shuffled_prctileB_summary(tempBoolIndex_validFOV))/2;

temp_1B_mean = mean(temp_1B);
temp_2B_mean = mean(temp_2B);
temp_3B_mean = mean(temp_3B);


% [~,temp_p_1] = ttest(temp_1,temp_1B);
% [~,temp_p_2] = ttest(temp_2,temp_2B);
% [~,temp_p_3] = ttest(temp_3,temp_3B);
[~,temp_p_1] = ttest(temp_1,temp_1B_mean);
[~,temp_p_2] = ttest(temp_2,temp_2B_mean);
[~,temp_p_3] = ttest(temp_3,temp_3B_mean);
% [~,temp_p_1] = ttest2(temp_1,temp_1B);
% [~,temp_p_2] = ttest2(temp_2,temp_2B);
% [~,temp_p_3] = ttest2(temp_3,temp_3B);


temp_1_mean = mean(temp_1);
temp_2_mean = mean(temp_2);
temp_3_mean = mean(temp_3);

% [~,temp_p_12] = ttest(temp_1,temp_2_mean);
% [~,temp_p_13] = ttest(temp_1,temp_3_mean);

[~,temp_p_12] = ttest(temp_1,temp_2);
[~,temp_p_13] = ttest(temp_1,temp_3);
[~,temp_p_23] = ttest(temp_2,temp_3);

if false
    %[~,temp_p_12] = ttest(temp_1,temp_2,'Tail','left');
    %[~,temp_p_13] = ttest(temp_1,temp_3,'Tail','left');
    
    temp_1B = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_shuffled_prctileB_summary;
    temp_2B = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_AC_shuffled_prctileB_summary;
    temp_3B = allPart_summary_mulFOV.part6_summary_mulFOV.centriodDis_all_BC_shuffled_prctileB_summary;
    
    temp1_boolIndex = temp_1<temp_1B;
    temp2_boolIndex = temp_2>temp_2B;
    temp3_boolIndex = temp_3>temp_3B;
    
    temp123_boolIndex = temp1_boolIndex & temp2_boolIndex & temp3_boolIndex;
    
    temp_123 = [temp_1;temp_2;temp_3];
    temp_123B = [temp_1B;temp_2B;temp_3B];
    
    %temp_1_T = temp_1';
    
    temp_123_round = round(temp_123*10)/10;
    temp_123B_round = round(temp_123B*10)/10;
    
end


% if if_plot == 1 && if_plot_part6 == 1
%     fig = figure('Name','','NumberTitle','off');
%     set(gcf,'Position',[50+0 50+0 240*0.80*1.1*1.1 336*1.11*0.5]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
%
%     t = tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
%     nexttile
%
%
%     temp_y_min = min([temp_1 temp_2_mean temp_3_mean]);
%     temp_y_max = max([temp_1 temp_2_mean temp_3_mean]);
%
%
%     temp_data = temp_1';
%     temp_label = repmat({'A'},length(temp_1),1);
%
%     temptemp_color1 = [1 1 1]*0.5;
%
%
%
%     violinplot(temp_data,temp_label,'ViolinColor',temptemp_color1,'BoxColor',[1 1 1]*0.2,...
%         'GroupOrder',{'A'});
%     hold on
%
%     h_line = [];
%
%     %h_line = [h_line plot([0 1.5],temp_1_mean*[1 1],'--','color',[0.4 0.4 0.4],'linewidth',1)];
%     h_line = [h_line plot([3 4],temp_1_mean*[1 1],'--','color',[0.4 0.4 0.4],'linewidth',1)];
%     hold on
%
%     h_line = [h_line plot([0 1.5],temp_2_mean*[1 1],'--','color',[0 0.4470 0.7410],'linewidth',1)];
%     hold on
%
%     h_line = [h_line plot([0 1.5],temp_3_mean*[1 1],'--','color',[0.8500 0.3250 0.0980],'linewidth',1)];
%     hold on
%
%
%
%     tempTxt = sprintf('');
%     if temp_p_12 < 0.001
%         tempTxt = sprintf('***');
%     elseif temp_p_12 < 0.01
%         tempTxt = sprintf('**');
%     elseif temp_p_12 < 0.05
%         tempTxt = sprintf('*');
%     end
%     text(1.75,temp_2_mean+2,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
%         'HorizontalAlignment','center');
%
%     tempTxt = sprintf('');
%     if temp_p_13 < 0.001
%         tempTxt = sprintf('***');
%     elseif temp_p_13 < 0.01
%         tempTxt = sprintf('**');
%     elseif temp_p_13 < 0.05
%         tempTxt = sprintf('*');
%     end
%     text(1.75,temp_3_mean-2,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
%         'HorizontalAlignment','center');
%
%     le = legend(h_line,'Precision & Meta','Precision & Baseline','Meta & Baseline','Location','northeast','fontsize',8);%10
%     le.ItemTokenSize = ones(1,2)*10;
%
%     set(gca,'linewidth',1.5)
%
%     xlim([0 2]);
%     ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*1.15]);%0.3
%     set(gca, 'FontSize', 12) %14
%     %set(gca,'YTick',0:0.1:1,'FontSize', 12);%12
%     set(gca,'box','off');% 取消右、上边框
%
%
%     xtl = [""]; %#ok<*NBRAK>
%     xt=get(gca,'XTick');
%     yt=get(gca,'YTick');
%     xtext_xp=xt;
%     xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.33;
%     text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',11);%9-->12
%
%     set(gca,'xticklabel','');
%
%
%     ylabel('Centroid distance', 'FontSize', 12, 'FontWeight', 'bold');
%     temp_title = title(sprintf('Spatial organization'),'fontsize',11);
%     temp_title.Interpreter = 'none';
%
% end


if if_plot == 1 && if_plot_part6 == 1
    color_precision = [252,141,89]/255;
    color_meta = [145,191,219]/255;
    color_baseline = 0.7*[255,255,191]/255;
    
    
    fig = figure('Name','Centroid distance','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80 336*1.11]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[50+0 50+0 240*0.80 336*1.11*0.67*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50+0 50+0 240*0.80 336*1.11*0.67*0.9*1.05*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_y_min = min([temp_1 temp_2 temp_3]);
    temp_y_max = max([temp_1 temp_2 temp_3]);
    
    plot([1-0.45 1+0.45],temp_1B_mean*[1 1],'color',[0.6350 0.0780 0.1840],'linewidth',4);
    hold on
    
    plot([2-0.45 2+0.45],temp_2B_mean*[1 1],'color',[0.6350 0.0780 0.1840],'linewidth',4);
    hold on
    
    plot([3-0.45 3+0.45],temp_3B_mean*[1 1],'color',[0.6350 0.0780 0.1840],'linewidth',4);
    hold on
    
    
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
    
    
    tempTxt = sprintf('');
    if temp_p_1 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p_1 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p_1 < 0.05
        tempTxt = sprintf('*');
    end
    text(1,temp_y_min-(temp_y_max-temp_y_min)*0.1,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    tempTxt = sprintf('');
    if temp_p_2 < 0.001
        tempTxt = sprintf('***');
    elseif temp_p_2 < 0.01
        tempTxt = sprintf('**');
    elseif temp_p_2 < 0.05
        tempTxt = sprintf('*');
    end
    text(2,temp_y_min-(temp_y_max-temp_y_min)*0.1,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
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
    text(3,temp_y_min-(temp_y_max-temp_y_min)*0.1,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    %plot([3-0.2 3+0.2],temp_y_max+(temp_y_max-temp_y_min)*0.05*[1 1],'color',[0.25 0.25 0.25],'linewidth',1.5);
    %hold on
    
    
    set(gca,'linewidth',1.5)
    
    xlim([0 4]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.1]);%0.3
    set(gca, 'FontSize', 10) %14
    %set(gca,'YTick',0:0.2:1,'FontSize', 12);%12
    set(gca,'box','off');% 取消右、上边框
    
    
    %xtl = ["Precision & Meta", "Precision & Baseline", "Meta & Baseline"];
    
    if if_tuning_location0_precision1 == 1
        temp_label1 = ['\color[rgb]',sprintf('{%f %f %f}',color_precision),'Strength',...
            '\color{black} & ',...
            '\color[rgb]',sprintf('{%f %f %f}',color_meta),'Meta'];
        
        temp_label2 = ['\color[rgb]',sprintf('{%f %f %f}',color_precision),'Strength',...
            '\color{black} & ',...
            '\color[rgb]',sprintf('{%f %f %f}',color_baseline),'Baseline'];
    elseif if_tuning_location0_precision1 == 0
        temp_label1 = ['\color[rgb]',sprintf('{%f %f %f}',color_precision),'Memory',...
            '\color{black} & ',...
            '\color[rgb]',sprintf('{%f %f %f}',color_meta),'Meta'];
        
        temp_label2 = ['\color[rgb]',sprintf('{%f %f %f}',color_precision),'Memory',...
            '\color{black} & ',...
            '\color[rgb]',sprintf('{%f %f %f}',color_baseline),'Baseline'];
    end
    
    temp_label3 = ['\color[rgb]',sprintf('{%f %f %f}',color_meta),'Meta',...
        '\color{black} & ',...
        '\color[rgb]',sprintf('{%f %f %f}',color_baseline),'Baseline'];
    
    xtl = [string(temp_label1),string(temp_label2),string(temp_label3)];
    
    
    
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.5;
    %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.61;
    %handle_txt = text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',8);%9-->12
    
    xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.31;
    handle_txt = text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Centroid distance', 'FontSize', 10, 'FontWeight', 'normal');
    temp_title = title(sprintf('Spatial organization'),'fontsize',10);
    temp_title.Interpreter = 'none';
    
    a = 1;
    
end


%% Spatial distance (new)
allPart_summary_mulFOV;

spatialDisNew_1to1_mean_summary = allPart_summary_mulFOV.part6_summary_mulFOV.spatialDisNew_1to1_mean_summary;
spatialDisNew_2to2_mean_summary = allPart_summary_mulFOV.part6_summary_mulFOV.spatialDisNew_2to2_mean_summary;
spatialDisNew_3to3_mean_summary = allPart_summary_mulFOV.part6_summary_mulFOV.spatialDisNew_3to3_mean_summary;

spatialDisNew_1toOthers_mean_summary = allPart_summary_mulFOV.part6_summary_mulFOV.spatialDisNew_1toOthers_mean_summary;
spatialDisNew_2toOthers_mean_summary = allPart_summary_mulFOV.part6_summary_mulFOV.spatialDisNew_2toOthers_mean_summary;
spatialDisNew_3toOthers_mean_summary = allPart_summary_mulFOV.part6_summary_mulFOV.spatialDisNew_3toOthers_mean_summary;

spatialCentriodDis_1toOthers_summary = allPart_summary_mulFOV.part6_summary_mulFOV.spatialCentriodDis_1toOthers_summary;
spatialCentriodDis_2toOthers_summary = allPart_summary_mulFOV.part6_summary_mulFOV.spatialCentriodDis_2toOthers_summary;
spatialCentriodDis_3toOthers_summary = allPart_summary_mulFOV.part6_summary_mulFOV.spatialCentriodDis_3toOthers_summary;

spatialDisNew_1to1_mean_summary_1d = nan(1,numFOV);
spatialDisNew_2to2_mean_summary_1d = nan(1,numFOV);
spatialDisNew_3to3_mean_summary_1d = nan(1,numFOV);

spatialDisNew_1toOthers_mean_summary_1d = nan(1,numFOV);
spatialDisNew_2toOthers_mean_summary_1d = nan(1,numFOV);
spatialDisNew_3toOthers_mean_summary_1d = nan(1,numFOV);


for tempi=1:numFOV
    %     spatialDisNew_1to1_mean_summary_1d(tempi) = mean(spatialDisNew_1to1_mean_summary{tempi});
    %     spatialDisNew_2to2_mean_summary_1d(tempi) = mean(spatialDisNew_2to2_mean_summary{tempi});
    %     spatialDisNew_3to3_mean_summary_1d(tempi) = mean(spatialDisNew_3to3_mean_summary{tempi});
    %
    %     spatialDisNew_1toOthers_mean_summary_1d(tempi) = mean(spatialDisNew_1toOthers_mean_summary{tempi});
    %     spatialDisNew_2toOthers_mean_summary_1d(tempi) = mean(spatialDisNew_2toOthers_mean_summary{tempi});
    %     spatialDisNew_3toOthers_mean_summary_1d(tempi) = mean(spatialDisNew_3toOthers_mean_summary{tempi});
    
    spatialDisNew_1to1_mean_summary_1d(tempi) = median(spatialDisNew_1to1_mean_summary{tempi});
    spatialDisNew_2to2_mean_summary_1d(tempi) = median(spatialDisNew_2to2_mean_summary{tempi});
    spatialDisNew_3to3_mean_summary_1d(tempi) = median(spatialDisNew_3to3_mean_summary{tempi});
    
    spatialDisNew_1toOthers_mean_summary_1d(tempi) = median(spatialDisNew_1toOthers_mean_summary{tempi});
    spatialDisNew_2toOthers_mean_summary_1d(tempi) = median(spatialDisNew_2toOthers_mean_summary{tempi});
    spatialDisNew_3toOthers_mean_summary_1d(tempi) = median(spatialDisNew_3toOthers_mean_summary{tempi});
    
end


temp_1 = spatialDisNew_1to1_mean_summary_1d;
temp_2 = spatialDisNew_2to2_mean_summary_1d;
temp_3 = spatialDisNew_3to3_mean_summary_1d;

if if_plot_distance_pairwise0_centroid1 == 0
    temp_x1 = spatialDisNew_1toOthers_mean_summary_1d;
    temp_x2 = spatialDisNew_2toOthers_mean_summary_1d;
    temp_x3 = spatialDisNew_3toOthers_mean_summary_1d;
elseif if_plot_distance_pairwise0_centroid1 == 1
    temp_x1 = spatialCentriodDis_1toOthers_summary;
    temp_x2 = spatialCentriodDis_2toOthers_summary;
    temp_x3 = spatialCentriodDis_3toOthers_summary;
end

temp_1_mean = median(temp_1);
temp_2_mean = median(temp_2);
temp_3_mean = median(temp_3);

temp_x1_mean = median(temp_x1);
temp_x2_mean = median(temp_x2);
temp_x3_mean = median(temp_x3);


% [~,temp_p_1] = ttest(temp_1,temp_x1);
% [~,temp_p_2] = ttest(temp_2,temp_x2);
% [~,temp_p_3] = ttest(temp_3,temp_x3);

[~,temp_p_1] = ttest2(temp_1,temp_x1);
[~,temp_p_2] = ttest2(temp_2,temp_x2);
[~,temp_p_3] = ttest2(temp_3,temp_x3);

% [~,temp_p_1] = ttest(temp_1,temp_x1,'tail','left');
% [~,temp_p_2] = ttest(temp_2,temp_x2,'tail','left');
% [~,temp_p_3] = ttest(temp_3,temp_x3,'tail','left');

if false
    %% Test
    spatialDisNew_1to1_mean_summary_collapsed = [];
    spatialDisNew_2to2_mean_summary_collapsed = [];
    spatialDisNew_3to3_mean_summary_collapsed = [];
    
    spatialDisNew_1toOthers_mean_summary_collapsed = [];
    spatialDisNew_2toOthers_mean_summary_collapsed = [];
    spatialDisNew_3toOthers_mean_summary_collapsed = [];
    
    %for tempi=1:numFOV
    if true
        tempi = 4;% Ding with 2 features match: FOV(2,5,11,12; 4)
        %tempi = 7;% Zelku with 2 features match: FOV(3,6,7)
        spatialDisNew_1to1_mean_summary_collapsed = [spatialDisNew_1to1_mean_summary_collapsed; spatialDisNew_1to1_mean_summary{tempi}];
        spatialDisNew_2to2_mean_summary_collapsed = [spatialDisNew_2to2_mean_summary_collapsed; spatialDisNew_2to2_mean_summary{tempi}];
        spatialDisNew_3to3_mean_summary_collapsed = [spatialDisNew_3to3_mean_summary_collapsed; spatialDisNew_3to3_mean_summary{tempi}];
        
        spatialDisNew_1toOthers_mean_summary_collapsed = [spatialDisNew_1toOthers_mean_summary_collapsed; spatialDisNew_1toOthers_mean_summary{tempi}];
        spatialDisNew_2toOthers_mean_summary_collapsed = [spatialDisNew_2toOthers_mean_summary_collapsed; spatialDisNew_2toOthers_mean_summary{tempi}];
        spatialDisNew_3toOthers_mean_summary_collapsed = [spatialDisNew_3toOthers_mean_summary_collapsed; spatialDisNew_3toOthers_mean_summary{tempi}];
    end
    
    %[~,temp_p_1] = ttest(spatialDisNew_1to1_mean_summary_collapsed,spatialDisNew_1toOthers_mean_summary_collapsed);
    %[~,temp_p_2] = ttest(spatialDisNew_2to2_mean_summary_collapsed,spatialDisNew_2toOthers_mean_summary_collapsed);
    %[~,temp_p_3] = ttest(spatialDisNew_3to3_mean_summary_collapsed,spatialDisNew_3toOthers_mean_summary_collapsed);
    
    [~,temp_p_1] = ttest2(spatialDisNew_1to1_mean_summary_collapsed,spatialDisNew_1toOthers_mean_summary_collapsed);
    [~,temp_p_2] = ttest2(spatialDisNew_2to2_mean_summary_collapsed,spatialDisNew_2toOthers_mean_summary_collapsed);
    [~,temp_p_3] = ttest2(spatialDisNew_3to3_mean_summary_collapsed,spatialDisNew_3toOthers_mean_summary_collapsed);
    
    fprintf('dis_1to1 = %.1f, dis_1toOthers = %.1f, p = %.3f.\n',...
        mean(spatialDisNew_1to1_mean_summary_collapsed),mean(spatialDisNew_1toOthers_mean_summary_collapsed),temp_p_1);
    
    fprintf('dis_2to2 = %.1f, dis_2toOthers = %.1f, p = %.3f.\n',...
        mean(spatialDisNew_2to2_mean_summary_collapsed),mean(spatialDisNew_2toOthers_mean_summary_collapsed),temp_p_2);
    
    fprintf('dis_3to3 = %.1f, dis_3toOthers = %.1f, p = %.3f.\n',...
        mean(spatialDisNew_3to3_mean_summary_collapsed),mean(spatialDisNew_3toOthers_mean_summary_collapsed),temp_p_3);
    
    temp_1 = spatialDisNew_1to1_mean_summary_collapsed';
    temp_2 = spatialDisNew_2to2_mean_summary_collapsed';
    temp_3 = spatialDisNew_3to3_mean_summary_collapsed';
    
    temp_x1 = spatialDisNew_1toOthers_mean_summary_collapsed';
    temp_x2 = spatialDisNew_2toOthers_mean_summary_collapsed';
    temp_x3 = spatialDisNew_3toOthers_mean_summary_collapsed';
    
    temp_1_mean = median(temp_1);
    temp_2_mean = median(temp_2);
    temp_3_mean = median(temp_3);
    
    temp_x1_mean = median(temp_x1);
    temp_x2_mean = median(temp_x2);
    temp_x3_mean = median(temp_x3);
    
end

%% Plot Within groups VS. Between groups' pc_distance distribution
if if_plot == 1 && if_plot_part6 == 1
    
    fprintf('dis_3to3 = %.1f, dis_3toOthers = %.1f, p = %.3f.\n',...
        mean(temp_3),mean(temp_x3),temp_p_3);
    
    fprintf('dis_1to1 = %.1f, dis_1toOthers = %.1f, p = %.3f.\n',...
        mean(temp_1),mean(temp_x1),temp_p_1);
    
    fprintf('dis_2to2 = %.1f, dis_2toOthers = %.1f, p = %.3f.\n',...
        mean(temp_2),mean(temp_x2),temp_p_2);
    
    
    
    fig = figure('Name','Spatial distance (new)','NumberTitle','off');
    %set(gcf,'Position',[50+0 50+0 240*0.80 336*1.11*0.67*0.9*1.05*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[50+0 50+0 240*0.80*0.9 336*1.11*0.67*0.9*1.05*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    
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
    
    
    %xtl = ["WM", "Meta", "BSL"];
    xtl = ["BSL", "WM", "Meta"];    
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.33;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',8);%9-->12
    
    set(gca,'xticklabel','');
    
    %yticks([0 40 80]);
    
    ylabel('Distance (um)', 'FontSize', 9, 'FontWeight', 'normal');
    %temp_title = title(sprintf('Centroid distance'),'fontsize',9);
    %temp_title.Interpreter = 'none';
    
end

%% Single neuron tuning property of memoryPrecision and choiceMemory (trial-level tuning)
if if_plot_3d_beta == 1 && if_plot_part6 == 1
    fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
    %set(gcf,'Position',[10 150 600 600]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 150 400 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    set(gcf,'Position',[10 150 150 150]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为pointt
    t = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');
    
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        x = beta_memoryPrecision_summary;
        y = beta_choiceMemory_summary;
        z = beta_choiceMemory_baseline_summary;
        tempStr = 'trial-level';
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        x = beta_seqPrecision_summary;
        y = beta_gProb_summary;
        z = beta_gProb_baseline_summary;
        tempStr = 'seq-level';
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        %x = beta_seqPrecision_summary;
        %y = beta_gProb_summary;
        %z = beta_choiceMemory_baseline_summary;
        
        x = beta_6loc_summary_norm;
        y = beta_gProb_summary_norm;
        z = beta_choiceMemory_baseline_summary_norm;
        
        tempStr = 'mix-level';
    end
    
    %         x = r2_seqPrecision_summary;
    %         %x = r2_6loc_summary;
    %         y = r2_gProb_summary;
    %         z = r2_choiceMemory_baseline_summary;
    
    %     x = abs(x);
    %     y = abs(y);
    %     z = abs(z);
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    [z_min,z_max] = bounds(z);
    
    tempSize = 3;%10
    scatter3(x,y,z,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    %     for tempi=1:clusterNum_tuning
    %         temptempBoolIndex = clusterID_tuning==tempi;
    %         temptempColor = multi_rgbColor(tempi,:);
    %         scatter3(x(temptempBoolIndex),y(temptempBoolIndex),z(temptempBoolIndex),tempSize,...
    %             'filled','MarkerFaceColor',temptempColor,'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    %         hold on
    %     end
    
    axis equal
    
    xticks([]);
    yticks([]);
    zticks([]);
    
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 10)
    set(gca,'box','off');% 取消右、上边框
    %     xlabel('Memory precision', 'FontSize', 10, 'FontWeight', 'normal');
    %     ylabel('Meta-memory', 'FontSize', 10, 'FontWeight', 'normal');
    %     zlabel('Baseline', 'FontSize', 10, 'FontWeight', 'normal');
end


% if false
%     % if true
%     temp_range_AB_multi = cell(numFOV_AB_summary,1);
%
%     for tempIndexFOV_AB=1:numFOV_AB_summary
%         %for tempIndexFOV_AB=1:2
%         currentSessionIndex_AB_summary = tempIndexFOV_AB;
%         currentSessionIndex_AB = currentSessionIndex_AB_summary;
%
%         currentSessionIndex_A = multiFOV_matrix_summary(currentSessionIndex_AB,1);
%         currentSessionIndex_B = multiFOV_matrix_summary(currentSessionIndex_AB,2);
%
%
%         currentSession = currentABSession_multi{currentSessionIndex_AB};
%
%         fprintf('currentSession = %s.\n',currentSession);
%         output_path = 'D:\twoPhotonData_motionCorrected\twoSessions';
%
%         temp_load = load([output_path '\' currentSession]);
%         decodingDataSimplified = temp_load.decodingDataSimplified_AB;
%
%
%         temp_range_B = FOVAllCellRange_multiFOV(currentSessionIndex_B,1):FOVAllCellRange_multiFOV(currentSessionIndex_B,2);
%         %temp_range_AB_raw = decodingDataSimplified.extraForMerged.tempMappingCellIndex(:,2)' + temp_range_B(1) - 1;
%
%         cellIndex_suite2p_B_raw = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range_B);
%
%         cellIndex_suite2p_B = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end);
%
%         temp_range_AB = nan(1,length(cellIndex_suite2p_B));
%         for tempi=1:length(cellIndex_suite2p_B)
%             temp_range_AB(tempi) = find(cellIndex_suite2p_B_raw==cellIndex_suite2p_B(tempi));
%         end
%         temp_range_AB = temp_range_AB + temp_range_B(1) - 1;
%
%         cellIndex_suite2p_B2 = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range_AB);
%
%         temp_range_AB_multi{tempIndexFOV_AB} = temp_range_AB;
%     end
%
%     allPart_summary_mulFOV.others = struct;
%     allPart_summary_mulFOV.others.temp_range_AB_multi = temp_range_AB_multi;
%
%     %beta_6loc_raw = glm_beta_6locMean_delay1Bin_multiFOV(temp_range_AB,:);
% end

if false
    %     temp_range_AB_multi = allPart_summary_mulFOV.others.temp_range_AB_multi;
    %
    %     temp_range_AB_multi_collapsed = [];
    %     for tempi=1:length(temp_range_AB_multi)
    %         temp_range_AB_multi_collapsed = [temp_range_AB_multi_collapsed temp_range_AB_multi{tempi}];
    %     end
    %     %temp_range_AB_multi_collapsed = temp_range_AB_multi{8};
    %
    %     glm_beta_6locMean_delay1Bin_multiFOV;
    %
    %     beta_6loc_raw = glm_beta_6locMean_delay1Bin_multiFOV(temp_range_AB_multi_collapsed,:);
    %     std_beta_6loc_raw = std(beta_6loc_raw,1,2);
    %     %beta_6loc;
    %
    %     beta_6loc = rescale(beta_6loc_raw,-1,1);
    %
    % %     std_beta_6loc = std_beta_6loc_raw-mean(std_beta_6loc_raw);
    % %     if abs(max(std_beta_6loc)) > abs(min(std_beta_6loc))
    % %         std_beta_6loc = std_beta_6loc./abs(max(std_beta_6loc));
    % %     else
    % %         std_beta_6loc = std_beta_6loc./abs(min(std_beta_6loc));
    % %     end
    %
    %     %std_beta_6loc = rescale(std_beta_6loc,-1,1);
    %     %std_beta_6loc = std_beta_6loc_raw;
    %
    %     std_beta_6loc = rescale(std_beta_6loc_raw,0,1);
    
    %if if_ePAIRS == 1
    if true
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_memoryPrecision_summary;
            y = beta_choiceMemory_summary;
            z = beta_choiceMemory_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_seqPrecision_summary;
            y = beta_gProb_summary;
            z = beta_gProb_baseline_summary;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_seqPrecision_summary;
            y = beta_gProb_summary;
            z = beta_choiceMemory_baseline_summary;
        end
        
        %         x = r2_seqPrecision_summary;
        %         y = r2_gProb_summary;
        %         z = r2_choiceMemory_baseline_summary;
        %         l = r2_6loc_summary;
        
        %x = r2_seqPrecision_summary_n11n;
        %y = r2_gProb_summary_n11n;
        %z = r2_choiceMemory_baseline_summary_n11n;
        %l = r2_6loc_summary_n11n;
        
        
        %l = std_beta_6loc;
        %l = beta_6loc_summary;
        l = beta_6loc_summary_peak;
        
        
        temp_vectors = [x,y,z,l];
        [clusteriness, p, dists] = pairsClusterTest_elliptical(temp_vectors);
        
        temp_angles = dists.data;
        temp_angles_bootstrap = dists.bootstrap;
        % temp_angles_bootstrap_median = median(temp_angles_bootstrap,1)';
        %temp_angles_bootstrap_median = median(temp_angles_bootstrap,2);
        % temp_angles_bootstrap_median = temp_angles_bootstrap(1,:)';
        temp_angles_bootstrap_median = reshape(temp_angles_bootstrap,[],1);
        
        
        %temp_vectors;
        %temp_vectors = [x,y,z,l];
        temp_vectors = [x,y,z];
        
        KList = 1:6;%1:6
        % criterion = 'CalinskiHarabasz';
        % criterion = 'DaviesBouldin';
        criterion = 'gap';
        % criterion = 'silhouette';
        
        % temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList);
        temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','sqEuclidean');
        %temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','cosine');
        % temp_eva = evalclusters(temp_vectors,'gmdistribution',criterion,'KList',KList,'Distance','sqEuclidean');
        
        %temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','sqEuclidean',...
        %   'ReferenceDistribution','uniform');
        %temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','sqEuclidean','SearchMethod','firstMaxSE');
        
        CriterionValues_tuning = temp_eva.CriterionValues;
        clusterID_tuning = temp_eva.OptimalY;
        clusterNum_tuning = temp_eva.OptimalK;
        
        
        %temp_vectors = [x,y,z];
        %temp_vectors = [x,y,z,l];
        
        
        temp_vectors = [x,y,z];
        [coeff3,score3,latent3,tsquared3,explained3,mu3] = pca(temp_vectors);
        %temp_eva3 = evalclusters(score3(:,1:2),'kmeans',criterion,'KList',KList,'Distance','sqEuclidean');
        %temp_eva3 = evalclusters(score3(:,1:2),'kmeans',criterion,'KList',KList,'Distance','cosine');
        %temp_eva3 = evalclusters(score3(:,1:2),'kmeans',criterion,'KList',KList,'Distance','correlation');
        temp_eva3 = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','sqEuclidean');
        
        clusterID_tuning_temp_eva3 = temp_eva3.OptimalY;
        clusterNum_tuning_temp_eva3 = temp_eva3.OptimalK;
        
        
        temp_vectors = [x,y,z,l];
        %temp_vectors = [y,z,l];
        [coeff4,score4,latent4,tsquared4,explained4,mu4] = pca(temp_vectors);
        %temp_eva4 = evalclusters(score4(:,1:3),'kmeans',criterion,'KList',KList,'Distance','sqEuclidean');
        %temp_eva4 = evalclusters(score4(:,1:3),'kmeans',criterion,'KList',KList,'Distance','cosine');
        %temp_eva4 = evalclusters(score4(:,1:3),'kmeans',criterion,'KList',KList,'Distance','correlation');
        temp_eva4 = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','sqEuclidean');
        
        clusterID_tuning_temp_eva4 = temp_eva4.OptimalY;
        clusterNum_tuning_temp_eva4 = temp_eva4.OptimalK;
        
        
        
        %% Plot score3
        fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
        %set(gcf,'Position',[510 450 240*3 240*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[510 450 720 290]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[710 450 240*0.95 240*1.05*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        temp1 = score3(:,1);
        temp2 = score3(:,2);
        
        tempSize = 4;%10
        
        [temp_r,temp_p] = corr(temp1,temp2);
        
        [x_min,x_max] = bounds(temp1);
        [y_min,y_max] = bounds(temp2);
        xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        xy_min = -1*xy_max;
        
        %scatter(temp1,temp2,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        %hold on
        
        for tempi=1:clusterNum_tuning_temp_eva3
            temptempBoolIndex = clusterID_tuning_temp_eva3==tempi;
            temptempColor = multi_rgbColor(tempi,:);
            scatter(temp1(temptempBoolIndex),temp2(temptempBoolIndex),tempSize,...
                'filled','MarkerFaceColor',temptempColor,'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        set(gca,'linewidth',1.5)
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('temp1', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('temp2', 'FontSize', 12, 'FontWeight', 'bold');
        temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
        temp_title.Interpreter = 'none';
        
        
        
        %% Plot score4
        fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
        set(gcf,'Position',[10 150 600 600]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');
        
        temp1 = score4(:,1);
        temp2 = score4(:,2);
        temp3 = score4(:,3);
        
        [x_min,x_max] = bounds(temp1);
        [y_min,y_max] = bounds(temp2);
        [z_min,z_max] = bounds(temp3);
        
        tempSize = 10;%10
        %scatter3(temp1,temp2,temp3,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        %hold on
        for tempi=1:clusterNum_tuning_temp_eva4
            temptempBoolIndex = clusterID_tuning_temp_eva4==tempi;
            temptempColor = multi_rgbColor(tempi,:);
            scatter3(temp1(temptempBoolIndex),temp2(temptempBoolIndex),temp3(temptempBoolIndex),tempSize,...
                'filled','MarkerFaceColor',temptempColor,'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
            hold on
        end
        
        
        
        axis equal
        
        set(gca,'linewidth',1.5)
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('temp1', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('temp2', 'FontSize', 12, 'FontWeight', 'bold');
        zlabel('temp3', 'FontSize', 12, 'FontWeight', 'bold');
        
        
    end
end

if false
    
    
    fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
    %set(gcf,'Position',[510 450 240*3 240*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[510 450 720 290]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[510 450 240*3*0.95 240*1.05*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    
    t.Title.String = sprintf('Multi-FOV (%s)\n Single neuron tuning property, n=%d',tempStr,length(beta_memoryPrecision_summary));
    t.Title.FontSize = 11;
    t.Title.Interpreter = 'none';
    
    tempSize = 4;%10
    
    nexttile
    
    x = beta_6loc_summary_peak;
    y = beta_seqPrecision_summary;
    
    [temp_r,temp_p] = corr(x,y);
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
    xy_min = -1*xy_max;
    
    scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    
    
    %axis equal
    
    set(gca,'linewidth',1.5)
    
    xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    
    set(gca, 'FontSize', 12)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Location', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Memory strength', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
    temp_title.Interpreter = 'none';
    
    
    nexttile
    
    x = beta_6loc_summary_peak;
    y = beta_choiceMemory_summary;
    
    
    [temp_r,temp_p] = corr(x,y);
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
    xy_min = -1*xy_max;
    
    
    scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    
    %axis equal
    
    set(gca,'linewidth',1.5)
    
    xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    
    set(gca, 'FontSize', 12)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Location', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
    temp_title.Interpreter = 'none';
    
    nexttile
    
    x = beta_6loc_summary_peak;
    y = beta_choiceMemory_baseline_summary;
    
    [temp_r,temp_p] = corr(x,y);
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
    xy_min = -1*xy_max;
    
    
    scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    
    %axis equal
    
    set(gca,'linewidth',1.5)
    
    xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    
    set(gca, 'FontSize', 12)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Location', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Baseline', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
    temp_title.Interpreter = 'none';
    
end


%% PCA along feature dimension
if if_PCA_feature0_neuron1 == 0
    % if false
    % if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
    %     x = beta_memoryPrecision_summary;
    %     y = beta_choiceMemory_summary;
    %     z = beta_choiceMemory_baseline_summary;
    % elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
    %     x = beta_seqPrecision_summary;
    %     y = beta_gProb_summary;
    %     z = beta_gProb_baseline_summary;
    % elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
    %     x = beta_seqPrecision_summary;
    %     y = beta_gProb_summary;
    %     %y = beta_choiceMemory_summary;
    %     z = beta_choiceMemory_baseline_summary;
    % end
    
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        x = beta_6loc_summary_peak;
        y = beta_choiceMemory_summary;
        z = beta_choiceMemory_baseline_summary;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        x = beta_6loc_summary_peak;
        y = beta_gProb_summary;
        z = beta_gProb_baseline_summary;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        x = beta_6loc_summary_peak;
        y = beta_gProb_summary;
        z = beta_choiceMemory_baseline_summary;
    end
    
    
    temp_vectors = [x,y,z];
    [coeff3,score3,latent3,tsquared3,explained3,mu3] = pca(temp_vectors);
    
    temp1 = nan(1,length(explained3));
    for tempi=1:length(temp1)
        temp1(tempi) = sum(explained3(1:tempi));
    end
    x_optimal = 2;
    y_optimal = temp1(x_optimal);
    
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name',' ','NumberTitle','off'); %#ok<*NASGU>
        set(gcf,'Position',[600 520 308*0.95*0.6*0.8*1.1*0.97 200*0.95*0.6*1.2*1.2*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        h_line = [];
        
        x = 1:length(temp1);
        y = temp1;
        h_line = [h_line plot(x,y,'LineWidth',1,'color',[0.25 0.25 0.25])];
        hold on
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        
        %     plot([1 1].*x_optimal,[y_min-(y_max-y_min)*0.8 temp1(x_optimal)],':','color',[0.25 0.25 0.25],'linewidth',1);
        %     hold on
        %
        %     plot([x_min-(x_max-x_min)*0.08 x_optimal],[1 1].*temp1(x_optimal),':','color',[0.25 0.25 0.25],'linewidth',1);
        %     hold on
        
        xticks(x);
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Number of PC', 'FontSize', 9, 'FontWeight', 'normal');
        %ylabel('Cumulative explained variance', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Explained variance', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('Optimal cluster number'),'FontSize',11);
        %temp_title.Interpreter = 'none';
    end
    
    
    temp1 = coeff3(:,1)';
    temp2 = coeff3(:,2)';
    temp3 = coeff3(:,3)';
    
    % [y_min,y_max] = bounds([temp1 temp2]);
    [y_min,y_max] = bounds([temp1 temp2 temp3]);
    
    
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name',' ','NumberTitle','off'); %#ok<*NASGU>
        %set(gcf,'Position',[600 520 308*0.95*0.6*0.8*1.1*0.97 200*0.95*0.6*1.2*1.2*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[600 520 150 153*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[600 520 150 153*0.8*0.63]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        h_line = [];
        
        x = 1:length(temp1);
        y = temp1;
        h_line = [h_line plot(x,y,'LineWidth',1,'color',[0.25 0.25 0.25])];
        hold on
        
        [x_min,x_max] = bounds(x);
        [temp_y_min,temp_y_max] = bounds(y);
        
        
        xticks(x);
        
        
        %xtl = ["Precision", "Meta", "BSL"];
        xtl = ["WM", "Meta", "BSL"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %     if if_monkey_D0_Z1 == 0
        %         %xtext_yp=(yt(1))*ones(1,length(xt))-(y_max-y_min)*0.70;
        %         xtext_yp=(yt(1))*ones(1,length(xt))-(y_max-y_min)*0.4;
        %     elseif if_monkey_D0_Z1 == 1
        %         xtext_yp=(yt(1))*ones(1,length(xt))-(y_max-y_min)*1.05;
        %     end
        xtext_yp=(y_min)*ones(1,length(xt))-(y_max-y_min)*0.60;
        
        %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%9-->12
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        yticks([0 1]);
        
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        xlim([0.5 3.5]);
        ylim([y_min-(y_max-y_min)*0.15 y_max+(y_max-y_min)*0.15]);
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        %xlabel('Number of PC', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Coefficients', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('PC1'),'FontSize',9);
        temp_title.Interpreter = 'none';
    end
    
    
    % temp1 = coeff3(:,2)';
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name',' ','NumberTitle','off'); %#ok<*NASGU>
        %set(gcf,'Position',[800 520 150 153*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[600 520 150 153*0.8*0.63]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        
        t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        h_line = [];
        
        x = 1:length(temp2);
        y = temp2;
        h_line = [h_line plot(x,y,'LineWidth',1,'color',[0.25 0.25 0.25])];
        hold on
        
        [x_min,x_max] = bounds(x);
        %[y_min,y_max] = bounds(y);
        
        
        xticks(x);
        
        
        %xtl = ["Precision", "Meta", "BSL"];
        xtl = ["WM", "Meta", "BSL"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %     if if_monkey_D0_Z1 == 0
        %         %xtext_yp=(yt(1))*ones(1,length(xt))-(y_max-y_min)*1.25;
        %         xtext_yp=(yt(1))*ones(1,length(xt))-(y_max-y_min)*0.95;
        %     elseif if_monkey_D0_Z1 == 1
        %         xtext_yp=(yt(1))*ones(1,length(xt))-(y_max-y_min)*0.60;
        %     end
        xtext_yp=(y_min)*ones(1,length(xt))-(y_max-y_min)*0.60;
        
        %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%9-->12
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        yticks([0 1]);
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        xlim([0.5 3.5]);
        ylim([y_min-(y_max-y_min)*0.15 y_max+(y_max-y_min)*0.15]);
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        %xlabel('Number of PC', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Coefficients', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('PC2'),'FontSize',9);
        temp_title.Interpreter = 'none';
    end
    
    
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name',' ','NumberTitle','off'); %#ok<*NASGU>
        %set(gcf,'Position',[1000 520 150 153*0.8]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[600 520 150 153*0.8*0.63]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        
        t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        
        h_line = [];
        
        x = 1:length(temp3);
        y = temp3;
        h_line = [h_line plot(x,y,'LineWidth',1,'color',[0.25 0.25 0.25])];
        hold on
        
        [x_min,x_max] = bounds(x);
        
        xticks(x);
        
        
        xtl = ["WM", "Meta", "BSL"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %     if if_monkey_D0_Z1 == 0
        %         xtext_yp=(yt(1))*ones(1,length(xt))-(y_max-y_min)*0.95;
        %     elseif if_monkey_D0_Z1 == 1
        %         xtext_yp=(yt(1))*ones(1,length(xt))-(y_max-y_min)*0.60;
        %     end
        xtext_yp=(y_min)*ones(1,length(xt))-(y_max-y_min)*0.60;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9-->12
        
        set(gca,'xticklabel','');
        
        yticks([0 1]);
        
        set(gca,'linewidth',1.5)
        xlim([0.5 3.5]);
        ylim([y_min-(y_max-y_min)*0.15 y_max+(y_max-y_min)*0.15]);
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        ylabel('Coefficients', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('PC3'),'FontSize',9);
        temp_title.Interpreter = 'none';
    end
    
end



%% PCA along neuron dimension
if if_plot == 1 & if_plot_part6 == 1 %#ok<*AND2>
    if if_PCA_feature0_neuron1 == 1
        
        %if_tuning_r20_beta1 = 0;
        
        if if_tuning_r20_beta1 == 0
            x_raw = r2_6loc_summary;
            y_raw = r2_gProb_summary;
            z_raw = r2_choiceMemory_baseline_summary;
            
        elseif if_tuning_r20_beta1 == 1
            x_raw = beta_6loc_summary_norm;
            y_raw = beta_gProb_summary_norm;
            z_raw = beta_choiceMemory_baseline_summary_norm;
        end
        
        if false
            %if true
            
            %temp1 = rand([size(x_raw,1),3]);
            %x_raw = temp1(:,1);
            %y_raw = temp1(:,2);
            %z_raw = temp1(:,3);
            
            x_raw = rand(size(x_raw));
            y_raw = rand(size(y_raw));
            z_raw = rand(size(z_raw));
            
            %x_raw = rand(size(x_raw));
            %y_raw = x_raw;
            %z_raw = rand(size(z_raw));
            
            %x_raw = rand(size(x_raw));
            %y_raw = x_raw*2+0.5;
            %z_raw = rand(size(z_raw));
        end
        
        %x = (x_raw-mean(x_raw))./std(x_raw);
        %y = (y_raw-mean(y_raw))./std(y_raw);
        %z = (z_raw-mean(z_raw))./std(z_raw);
        
        x = x_raw;
        y = y_raw;
        z = z_raw;
        
        temp_singleFOVIndex = 4;%8,4,7
        if if_PCA_singleFOV0_multiFOV1 == 0
            temp_FOVIndex_mulFOV;
            temptempBoolIndex = temp_FOVIndex_mulFOV == temp_singleFOVIndex;
            x = x(temptempBoolIndex);
            y = y(temptempBoolIndex);
            z = z(temptempBoolIndex);
            
        end
        
        
        %     x = [ones(1000,1);ones(2000,1)];
        %     y = [ones(1000,1);ones(1000,1);ones(1000,1)];
        %     z = [ones(2000,1);ones(1000,1)];
        
        temp_vectors = [x,y,z];
        [coeff3,score3,latent3,tsquared3,explained3,mu3] = pca(temp_vectors');
        
        temp1 = nan(1,length(explained3));
        for tempi=1:length(temp1)
            temp1(tempi) = sum(explained3(1:tempi));
        end
        
        %if if_plot == 1 && if_plot_part6 == 1
        if false
            fig = figure('Name',' ','NumberTitle','off'); %#ok<*NASGU>
            set(gcf,'Position',[600 200 308*0.95*0.6*0.8*1.1*0.97 200*0.95*0.6*1.2*1.2*0.93]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
            
            nexttile
            
            h_line = [];
            
            x = 1:length(temp1);
            y = temp1;
            h_line = [h_line plot(x,y,'LineWidth',1,'color',[0.25 0.25 0.25])];
            hold on
            
            [x_min,x_max] = bounds(x);
            [y_min,y_max] = bounds(y);
            
            
            xticks(x);
            
            set(gca,'linewidth',1.5)
            xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
            ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
            set(gca, 'FontSize', 9)
            set(gca,'box','off');% 取消右、上边框
            xlabel('Number of PC', 'FontSize', 9, 'FontWeight', 'normal');
            ylabel('Explained variance', 'FontSize', 9, 'FontWeight', 'normal');
        end
        
        
        temp_vectors;
        
        %% Shuffle
        if true
            temp_shuffleNum = 10;%100
            temp_vectors_shuffled = nan([size(temp_vectors),temp_shuffleNum]);
            temp_roiNum = size(temp_vectors,1);
            parfor tempi=1:temp_shuffleNum
                temptemp_temp_vectors_shuffled = temp_vectors;
                for tempj=1:temp_roiNum
                    temp1 = randperm(3);
                    temptemp_temp_vectors_shuffled(tempj,:) = temp_vectors(tempj,temp1);
                    
                    %temp1 = ones(1,3)*randperm(3,1);
                    %temptemp_temp_vectors_shuffled(tempj,:) = temp_vectors(tempj,temp1);
                    
                    temp_vectors_shuffled(tempj,:,tempi) = temptemp_temp_vectors_shuffled(tempj,:);
                end
            end
            
            temp_vectors_shuffled_mean = mean(temp_vectors_shuffled,3);
            
            temp_score3_shuffled_mean = (temp_vectors_shuffled_mean'-mu3)*coeff3;
            
            temp_dis12 = norm(score3(1,:)-score3(2,:));
            temp_dis13 = norm(score3(1,:)-score3(3,:));
            temp_dis23 = norm(score3(2,:)-score3(3,:));
            
            temp_dis12_shuffled = nan(1,temp_shuffleNum);
            temp_dis13_shuffled = nan(1,temp_shuffleNum);
            temp_dis23_shuffled = nan(1,temp_shuffleNum);
            
            temp_score3_shuffled = nan([size(temp_score3_shuffled_mean),temp_shuffleNum]);
            for tempi=1:temp_shuffleNum
                temp_score3_shuffled(:,:,tempi) = (temp_vectors_shuffled(:,:,tempi)'-mu3)*coeff3;
                
                temp_dis12_shuffled(tempi) = norm(temp_score3_shuffled(1,:,tempi)-temp_score3_shuffled(2,:,tempi));
                temp_dis13_shuffled(tempi) = norm(temp_score3_shuffled(1,:,tempi)-temp_score3_shuffled(3,:,tempi));
                temp_dis23_shuffled(tempi) = norm(temp_score3_shuffled(2,:,tempi)-temp_score3_shuffled(3,:,tempi));
            end
            
            [~,temp_p_dis12] = ttest(temp_dis12_shuffled,temp_dis12);
            [~,temp_p_dis13] = ttest(temp_dis13_shuffled,temp_dis13);
            [~,temp_p_dis23] = ttest(temp_dis23_shuffled,temp_dis23);
            
            %fprintf('p_dis12 = %.3f, p_dis13 = %.3f, p_dis23 = %.3f.\n',temp_p_dis12,temp_p_dis13,temp_p_dis23);
            
            
            
            [coeff3_shuffled,score3_shuffled,latent3_shuffled,tsquared3_shuffled,explained3_shuffled,mu3_shuffled] = pca(temp_vectors_shuffled_mean');
            score3_projToShuffle = (temp_vectors'-mu3_shuffled)*coeff3_shuffled;
            
        end
        
        
        
        %% Resample
        if true
            temp_resampleNum_B = 1000;%100
            %temp_resampleNum_B = size(temp_vectors,1);
            
            temp_resample_lesionProb = 0.95;%0.1,0.5,0.95,0.99,0.95
            
            temp_vectors_resampledB = nan([size(temp_vectors),temp_resampleNum_B]);
            temp_roiNum = size(temp_vectors,1);
            parfor tempi=1:temp_resampleNum_B
                temptemp_temp_vectors_resampledB = temp_vectors;
                
                temptempBoolIndex = false(temp_roiNum,1);
                
                temp1 = sort(randperm(temp_roiNum,round(temp_roiNum*temp_resample_lesionProb)));
                %temp1 = sort(randperm(temp_roiNum,temp_roiNum-1));
                temptempBoolIndex(temp1) = true;
                
                temptemp_temp_vectors_resampledB = temp_vectors;
                
                temp2 = mean(temp_vectors,2);
                
                %temp2 = temp2 .* 0;
                %temp2 = rand(size(temp_vectors));
                
                temptemp_temp_vectors_resampledB(temptempBoolIndex,:) = temp2(temptempBoolIndex) * [1,1,1];
                
                temp_vectors_resampledB(:,:,tempi) = temptemp_temp_vectors_resampledB;
                
            end
            
            %         parfor tempi=1:temp_resampleNum_B
            %             temptemp_temp_vectors_resampledB = temp_vectors;
            %
            %             temptempBoolIndex = true(temp_roiNum,1);
            %
            %             %temp1 = sort(randperm(temp_roiNum,round(temp_roiNum*(1-temp_resampleProb))));
            %             %temptempBoolIndex(temp1) = true;
            %
            %             temptempBoolIndex(tempi) = false;
            %
            %             temptemp_temp_vectors_resampledB = temp_vectors;
            %
            %             temp2 = mean(temp_vectors,2);
            %
            %             temp2 = temp2 .* 0;
            %
            %             temptemp_temp_vectors_resampledB(temptempBoolIndex,:) = temp2(temptempBoolIndex) * [1,1,1];
            %
            %             temp_vectors_resampledB(:,:,tempi) = temptemp_temp_vectors_resampledB;
            %
            %         end
            
            temp_vectors_resampledB_mean = mean(temp_vectors_resampledB,3);
            
            temp_score3_resampledB_mean = (temp_vectors_resampledB_mean'-mu3)*coeff3;
            
            temp_score3_resampledB = nan([size(temp_score3_resampledB_mean),temp_resampleNum_B]);
            for tempi=1:temp_resampleNum_B
                temp_score3_resampledB(:,:,tempi) = (temp_vectors_resampledB(:,:,tempi)'-mu3)*coeff3;
            end
            
            
            %% PC_distance within groups & between groups
            % temp_pcDisNew_1to1_mean
            temp_coordinate = squeeze(temp_score3_resampledB(1,:,:))';
            pcDisNew  = nan(temp_resampleNum_B,temp_resampleNum_B);
            for tempi=1:temp_resampleNum_B
                temp_1 = temp_coordinate(tempi,:);
                for tempj=1:temp_resampleNum_B
                    temp_2 = temp_coordinate(tempj,:);
                    temp_dis = norm(temp_1-temp_2);
                    pcDisNew(tempi,tempj) = temp_dis;
                end
            end
            temp_pcDisNew_1to1 = pcDisNew;
            temp_pcDisNew_1to1_mean = mean(temp_pcDisNew_1to1,2);
            
            % temp_pcDisNew_2to2_mean
            temp_coordinate = squeeze(temp_score3_resampledB(2,:,:))';
            pcDisNew  = nan(temp_resampleNum_B,temp_resampleNum_B);
            for tempi=1:temp_resampleNum_B
                temp_1 = temp_coordinate(tempi,:);
                for tempj=1:temp_resampleNum_B
                    temp_2 = temp_coordinate(tempj,:);
                    temp_dis = norm(temp_1-temp_2);
                    pcDisNew(tempi,tempj) = temp_dis;
                end
            end
            temp_pcDisNew_2to2 = pcDisNew;
            temp_pcDisNew_2to2_mean = mean(temp_pcDisNew_2to2,2);
            
            % temp_pcDisNew_3to3_mean
            temp_coordinate = squeeze(temp_score3_resampledB(3,:,:))';
            pcDisNew  = nan(temp_resampleNum_B,temp_resampleNum_B);
            for tempi=1:temp_resampleNum_B
                temp_1 = temp_coordinate(tempi,:);
                for tempj=1:temp_resampleNum_B
                    temp_2 = temp_coordinate(tempj,:);
                    temp_dis = norm(temp_1-temp_2);
                    pcDisNew(tempi,tempj) = temp_dis;
                end
            end
            temp_pcDisNew_3to3 = pcDisNew;
            temp_pcDisNew_3to3_mean = mean(temp_pcDisNew_3to3,2);
            
            %
            temp_pcDisNew_1to1_mean;
            temp_pcDisNew_2to2_mean;
            temp_pcDisNew_3to3_mean;
            
            centriod_1 = temp_score3_resampledB_mean(1,:);
            centriod_2 = temp_score3_resampledB_mean(2,:);
            centriod_3 = temp_score3_resampledB_mean(3,:);
            
            pcCentriodDis_1toOthers = norm(centriod_1-(centriod_2+centriod_3)/2);
            pcCentriodDis_2toOthers = norm(centriod_2-(centriod_1+centriod_3)/2);
            pcCentriodDis_3toOthers = norm(centriod_3-(centriod_1+centriod_2)/2);
            
            [~,temp_p_11_1Others_pc] = ttest(temp_pcDisNew_1to1_mean,pcCentriodDis_1toOthers);
            [~,temp_p_22_2Others_pc] = ttest(temp_pcDisNew_2to2_mean,pcCentriodDis_2toOthers);
            [~,temp_p_33_3Others_pc] = ttest(temp_pcDisNew_3to3_mean,pcCentriodDis_3toOthers);
            
            %         fprintf('pcDis_1to1 = %.2f, pcCentriodDis_1toOthers = %.2f, p = %.3f.\n',...
            %             mean(temp_pcDisNew_1to1_mean),pcCentriodDis_1toOthers,temp_p_11_1Others_pc);
            %
            %         fprintf('pcDis_2to2 = %.2f, pcCentriodDis_2toOthers = %.2f, p = %.3f.\n',...
            %             mean(temp_pcDisNew_2to2_mean),pcCentriodDis_2toOthers,temp_p_22_2Others_pc);
            %
            %         fprintf('pcDis_3to3 = %.2f, pcCentriodDis_3toOthers = %.2f, p = %.3f.\n',...
            %             mean(temp_pcDisNew_3to3_mean),pcCentriodDis_3toOthers,temp_p_33_3Others_pc);
            
            
        end
        
        %% Across FOV
        if if_PCA_singleFOV0_multiFOV1 == 1
            
            temp_FOVNum = length(allPart_summary_mulFOV.temp_cellIndexMapping_AB_mulFOV);
            temp_FOVIndex_mulFOV;
            
            temp_vectors_eachFOV = nan([size(temp_vectors),temp_FOVNum]);
            temp_roiNum = size(temp_vectors,1);
            for tempi=1:temp_FOVNum
                temptemp_temp_vectors_eachFOV = temp_vectors;
                
                temptempBoolIndex = temp_FOVIndex_mulFOV == tempi;
                
                % each FOV
                temptempBoolIndex2 = ~temptempBoolIndex;
                
                % leave one FOV out
                %temptempBoolIndex2 = temptempBoolIndex;
                
                temp2 = mean(temp_vectors,2);
                
                %temp2 = temp2 .* 0;
                %temp2 = rand(size(temp_vectors));
                
                
                temptemp_temp_vectors_eachFOV(temptempBoolIndex2,:) = temp2(temptempBoolIndex2) * [1,1,1];
                
                temp_vectors_eachFOV(:,:,tempi) = temptemp_temp_vectors_eachFOV;
            end
            
            
            temp_vectors_eachFOV_mean = mean(temp_vectors_eachFOV,3);
            
            temp_score3_eachFOV_mean = (temp_vectors_eachFOV_mean'-mu3)*coeff3;
            
            temp_score3_eachFOV = nan([size(temp_score3_eachFOV_mean),temp_FOVNum]);
            for tempi=1:temp_FOVNum
                temp_score3_eachFOV(:,:,tempi) = (temp_vectors_eachFOV(:,:,tempi)'-mu3)*coeff3;
            end
            
            
            %% PC_distance_eachFOV within groups & between groups
            % temp_pcDisNewFOV_1to1_mean
            temp_coordinate = squeeze(temp_score3_eachFOV(1,:,:))';
            pcDisNewFOV  = nan(temp_FOVNum,temp_FOVNum);
            for tempi=1:temp_FOVNum
                temp_1 = temp_coordinate(tempi,:);
                for tempj=1:temp_FOVNum
                    temp_2 = temp_coordinate(tempj,:);
                    temp_dis = norm(temp_1-temp_2);
                    pcDisNewFOV(tempi,tempj) = temp_dis;
                end
            end
            temp_pcDisNewFOV_1to1 = pcDisNewFOV;
            temp_pcDisNewFOV_1to1_mean = mean(temp_pcDisNewFOV_1to1,2);
            
            % temp_pcDisNewFOV_2to2_mean
            temp_coordinate = squeeze(temp_score3_eachFOV(2,:,:))';
            pcDisNewFOV  = nan(temp_FOVNum,temp_FOVNum);
            for tempi=1:temp_FOVNum
                temp_1 = temp_coordinate(tempi,:);
                for tempj=1:temp_FOVNum
                    temp_2 = temp_coordinate(tempj,:);
                    temp_dis = norm(temp_1-temp_2);
                    pcDisNewFOV(tempi,tempj) = temp_dis;
                end
            end
            temp_pcDisNewFOV_2to2 = pcDisNewFOV;
            temp_pcDisNewFOV_2to2_mean = mean(temp_pcDisNewFOV_2to2,2);
            
            % temp_pcDisNewFOV_3to3_mean
            temp_coordinate = squeeze(temp_score3_eachFOV(3,:,:))';
            pcDisNewFOV  = nan(temp_FOVNum,temp_FOVNum);
            for tempi=1:temp_FOVNum
                temp_1 = temp_coordinate(tempi,:);
                for tempj=1:temp_FOVNum
                    temp_2 = temp_coordinate(tempj,:);
                    temp_dis = norm(temp_1-temp_2);
                    pcDisNewFOV(tempi,tempj) = temp_dis;
                end
            end
            temp_pcDisNewFOV_3to3 = pcDisNewFOV;
            temp_pcDisNewFOV_3to3_mean = mean(temp_pcDisNewFOV_3to3,2);
            
            %
            temp_pcDisNewFOV_1to1_mean;
            temp_pcDisNewFOV_2to2_mean;
            temp_pcDisNewFOV_3to3_mean;
            
            centriod_1 = temp_score3_eachFOV_mean(1,:);
            centriod_2 = temp_score3_eachFOV_mean(2,:);
            centriod_3 = temp_score3_eachFOV_mean(3,:);
            
            pcCentriodDisFOV_1toOthers = norm(centriod_1-(centriod_2+centriod_3)/2);
            pcCentriodDisFOV_2toOthers = norm(centriod_2-(centriod_1+centriod_3)/2);
            pcCentriodDisFOV_3toOthers = norm(centriod_3-(centriod_1+centriod_2)/2);
            
            [~,temp_p_11_1Others_pcFOV] = ttest(temp_pcDisNewFOV_1to1_mean,pcCentriodDisFOV_1toOthers,'tail','left');
            [~,temp_p_22_2Others_pcFOV] = ttest(temp_pcDisNewFOV_2to2_mean,pcCentriodDisFOV_2toOthers,'tail','left');
            [~,temp_p_33_3Others_pcFOV] = ttest(temp_pcDisNewFOV_3to3_mean,pcCentriodDisFOV_3toOthers,'tail','left');
            
            %         fprintf('pcDis_1to1 = %.2f, pcCentriodDisFOV_1toOthers = %.2f, p = %.3f.\n',...
            %             mean(temp_pcDisNewFOV_1to1_mean),pcCentriodDisFOV_1toOthers,temp_p_11_1Others_pcFOV);
            %
            %         fprintf('pcDis_2to2 = %.2f, pcCentriodDisFOV_2toOthers = %.2f, p = %.3f.\n',...
            %             mean(temp_pcDisNewFOV_2to2_mean),pcCentriodDisFOV_2toOthers,temp_p_22_2Others_pcFOV);
            %
            %         fprintf('pcDis_3to3 = %.2f, pcCentriodDisFOV_3toOthers = %.2f, p = %.3f.\n',...
            %             mean(temp_pcDisNewFOV_3to3_mean),pcCentriodDisFOV_3toOthers,temp_p_33_3Others_pcFOV);
            
        end
        
        
        
        if false
            a = [ones(1000,1);zeros(2000,1)];
            b = [zeros(1000,1);ones(1000,1);zeros(1000,1)];
            c = [zeros(2000,1);ones(1000,1)];
            
            temp_degree_ab = subspace(a,b)*180/pi;
            temp_degree_ac = subspace(a,c)*180/pi;
            temp_degree_bc = subspace(b,c)*180/pi;
            
            
            a = [1 0 1]';
            b = [0 1 0]';
            
            temp_degree_ab = subspace(a,b)*180/pi
            
        end
        
        
        %% Angles, resampling
        if true
            a = temp_vectors(:,1);
            b = temp_vectors(:,2);
            c = temp_vectors(:,3);
            
            temp_degree_ab = subspace(a,b)*180/pi;
            temp_degree_ac = subspace(a,c)*180/pi;
            temp_degree_bc = subspace(b,c)*180/pi;
            
            temp_resampleNum = 100;
            temp_vectors_resampled = nan([size(temp_vectors),temp_resampleNum]);
            for tempi=1:temp_resampleNum
                temptempIndex = sort(randi(size(temp_vectors,1),[1,size(temp_vectors,1)]));
                
                temp_vectors_resampled(:,:,tempi) = temp_vectors(temptempIndex,:);
            end
            
            
            temp_degree_ab_resampled = nan(1,temp_resampleNum);
            temp_degree_ac_resampled = nan(1,temp_resampleNum);
            temp_degree_bc_resampled = nan(1,temp_resampleNum);
            temp_degree_aa_resampled = nan(1,temp_resampleNum);
            temp_degree_bb_resampled = nan(1,temp_resampleNum);
            temp_degree_cc_resampled = nan(1,temp_resampleNum);
            
            %a0 = temp_vectors(:,1);
            %b0 = temp_vectors(:,2);
            %c0 = temp_vectors(:,3);
            
            for tempi=1:temp_resampleNum
                
                a = temp_vectors_resampled(:,1,tempi);
                b = temp_vectors_resampled(:,2,tempi);
                c = temp_vectors_resampled(:,3,tempi);
                
                temp_degree_ab_resampled(tempi) = subspace(a,b)*180/pi;
                temp_degree_ac_resampled(tempi) = subspace(a,c)*180/pi;
                temp_degree_bc_resampled(tempi) = subspace(b,c)*180/pi;
                temp_degree_aa_resampled(tempi) = subspace(a,a)*180/pi;
                temp_degree_bb_resampled(tempi) = subspace(b,b)*180/pi;
                temp_degree_cc_resampled(tempi) = subspace(c,c)*180/pi;
                
            end
            temp_degree_ab_resampled_mean = mean(temp_degree_ab_resampled);
            temp_degree_ac_resampled_mean = mean(temp_degree_ac_resampled);
            temp_degree_bc_resampled_mean = mean(temp_degree_bc_resampled);
            temp_degree_aa_resampled_mean = mean(temp_degree_aa_resampled);
            temp_degree_bb_resampled_mean = mean(temp_degree_bb_resampled);
            temp_degree_cc_resampled_mean = mean(temp_degree_cc_resampled);
            
            [~,temp_p_ab_resampled] = ttest(temp_degree_ab_resampled,0);
            [~,temp_p_ac_resampled] = ttest(temp_degree_ac_resampled,0);
            [~,temp_p_bc_resampled] = ttest(temp_degree_bc_resampled,0);
            
            %fprintf('p_degree_ab = %.3f, p_degree_ac = %.3f, p_degree_bc = %.3f.\n',...
            %    temp_p_ac_resampled,temp_p_bc_resampled,temp_p_dis23);
            
        end
        
        
        
        %% Angles, each FOV
        if if_PCA_singleFOV0_multiFOV1 == 1
            a = temp_vectors(:,1);
            b = temp_vectors(:,2);
            c = temp_vectors(:,3);
            
            temp_degree_ab = subspace(a,b)*180/pi;
            temp_degree_ac = subspace(a,c)*180/pi;
            temp_degree_bc = subspace(b,c)*180/pi;
            
            temp_FOVNum = length(allPart_summary_mulFOV.temp_cellIndexMapping_AB_mulFOV);
            temp_FOVIndex_mulFOV;
            
            %temp_vectors_eachFOV = nan([size(temp_vectors),temp_FOVNum]);
            temp_vectors_eachFOV = cell(temp_FOVNum,1);
            for tempi=1:temp_FOVNum
                %temptempIndex = sort(randi(size(temp_vectors,1),[1,size(temp_vectors,1)]));
                
                temptempBoolIndex = temp_FOVIndex_mulFOV == tempi;
                
                % each FOV
                temp_vectors_eachFOV{tempi} = temp_vectors(temptempBoolIndex,:);
                
                % leave one FOV out
                %temp_vectors_eachFOV{tempi} = temp_vectors(~temptempBoolIndex,:);
            end
            
            
            temp_degree_ab_eachFOV = nan(1,temp_FOVNum);
            temp_degree_ac_eachFOV = nan(1,temp_FOVNum);
            temp_degree_bc_eachFOV = nan(1,temp_FOVNum);
            temp_degree_aa_eachFOV = nan(1,temp_FOVNum);
            temp_degree_bb_eachFOV = nan(1,temp_FOVNum);
            temp_degree_cc_eachFOV = nan(1,temp_FOVNum);
            
            for tempi=1:temp_FOVNum
                
                a = temp_vectors_eachFOV{tempi}(:,1);
                b = temp_vectors_eachFOV{tempi}(:,2);
                c = temp_vectors_eachFOV{tempi}(:,3);
                
                temp_degree_ab_eachFOV(tempi) = subspace(a,b)*180/pi;
                temp_degree_ac_eachFOV(tempi) = subspace(a,c)*180/pi;
                temp_degree_bc_eachFOV(tempi) = subspace(b,c)*180/pi;
                temp_degree_aa_eachFOV(tempi) = subspace(a,a)*180/pi;
                temp_degree_bb_eachFOV(tempi) = subspace(b,b)*180/pi;
                temp_degree_cc_eachFOV(tempi) = subspace(c,c)*180/pi;
                
            end
            temp_degree_ab_eachFOV_mean = mean(temp_degree_ab_eachFOV);
            temp_degree_ac_eachFOV_mean = mean(temp_degree_ac_eachFOV);
            temp_degree_bc_eachFOV_mean = mean(temp_degree_bc_eachFOV);
            temp_degree_aa_eachFOV_mean = mean(temp_degree_aa_eachFOV);
            temp_degree_bb_eachFOV_mean = mean(temp_degree_bb_eachFOV);
            temp_degree_cc_eachFOV_mean = mean(temp_degree_cc_eachFOV);
            
            [~,temp_p_ab_eachFOV] = ttest(temp_degree_ab_eachFOV,0);
            [~,temp_p_ac_eachFOV] = ttest(temp_degree_ac_eachFOV,0);
            [~,temp_p_bc_eachFOV] = ttest(temp_degree_bc_eachFOV,0);
            
            %[~,temp_p_ab_eachFOV_90] = ttest(temp_degree_ab_eachFOV,90);
            %[~,temp_p_ac_eachFOV_90] = ttest(temp_degree_ac_eachFOV,90);
            %[~,temp_p_bc_eachFOV_90] = ttest(temp_degree_bc_eachFOV,90);
            
            %fprintf('p_degree_ab = %.3f, p_degree_ac = %.3f, p_degree_bc = %.3f.\n',...
            %   temp_p_ac_eachFOV,temp_p_bc_eachFOV,temp_p_dis23);
            
        end
        
        
        %% PCA along neuron dimension
        if if_plot == 1 && if_plot_part6 == 1
            fig = figure('Name','PCA along neuron dimension','NumberTitle','off');
            set(gcf,'Position',[300 200 240*0.85*0.95*0.95 240*1.05*0.67*1.05*0.95*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
            
            
            nexttile
            
            x = score3(:,1);
            y = score3(:,2);
            
            x_shuffled = temp_score3_shuffled_mean(:,1);
            y_shuffled = temp_score3_shuffled_mean(:,2);
            
            
            %                 x = score3_projToShuffle(:,1);
            %                 y = score3_projToShuffle(:,2);
            %
            %                 x_shuffled = score3_shuffled(:,1);
            %                 y_shuffled = score3_shuffled(:,2);
            
            
            [temp_r,temp_p] = corr(x,y);
            
            [x_min,x_max] = bounds([x;x_shuffled]);
            [y_min,y_max] = bounds([y;y_shuffled]);
            
            tempSize = 30;%10,4
            temp_MarkerFaceAlpha = 1;
            for tempi=1:3
                if tempi == 1
                    temptemp_color = color_memoryQuality;
                elseif tempi == 2
                    temptemp_color = color_meta;
                elseif tempi == 3
                    temptemp_color = color_prior;
                end
                
                scatter(x(tempi),y(tempi),tempSize,'filled','MarkerFaceColor',temptemp_color,...
                    'MarkerFaceAlpha',temp_MarkerFaceAlpha,'MarkerEdgeColor',[0 0 0]);
                hold on
                
                %scatter(x_shuffled(tempi),y_shuffled(tempi),tempSize/2,'filled','MarkerFaceColor',temptemp_color,...
                %    'MarkerFaceAlpha',0.5,'MarkerEdgeColor',[0 0 0]);
                %hold on
            end
            
            
            %axis equal
            
            set(gca,'linewidth',1.5)
            xlim([x_min-(x_max-x_min)*0.15 x_max+(x_max-x_min)*0.15]);
            ylim([y_min-(y_max-y_min)*0.15 y_max+(y_max-y_min)*0.15]);
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            %xlabel('PC1', 'FontSize', 9, 'FontWeight', 'normal');
            %ylabel('PC2', 'FontSize', 9, 'FontWeight', 'normal');
            
            xlabel(sprintf('PC1 (%.0f%%)',explained3(1)), 'FontSize', 9, 'FontWeight', 'normal');
            ylabel(sprintf('PC2 (%.0f%%)',explained3(2)), 'FontSize', 9, 'FontWeight', 'normal');
            
            %temp_title = title(sprintf('r=%.3f, %s',temp_r,tempTxt), 'FontSize', 9, 'FontWeight', 'normal');
            %temp_title.Interpreter = 'none';
        end
        
        %% PCA along neuron dimension (lesionProb)
        if if_plot == 1 && if_plot_part6 == 1
            fig = figure('Name','PCA along neuron dimension (lesionProb)','NumberTitle','off');
            %set(gcf,'Position',[500 200 240*0.85*0.95*0.95 240*1.05*0.67*1.05*0.95*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[500 200 240*0.85*0.95*0.95*0.93 240*1.05*0.67*1.05*0.95*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
            
            
            nexttile
            
            x = squeeze(temp_score3_resampledB(:,1,:));
            y = squeeze(temp_score3_resampledB(:,2,:));
            
            x_mean = temp_score3_resampledB_mean(:,1);
            y_mean = temp_score3_resampledB_mean(:,2);
            
            
            
            [x_min,x_max] = bounds(x(:));
            [y_min,y_max] = bounds(y(:));
            
            temp_h = [];
            
            tempSize = 8;%30
            temp_MarkerFaceAlpha = 0.4;%1
            for tempi=1:3
                if tempi == 1
                    temptemp_color = color_memoryQuality;
                elseif tempi == 2
                    temptemp_color = color_meta;
                elseif tempi == 3
                    temptemp_color = color_prior;
                end
                
                %scatter(x(tempi,:),y(tempi,:),tempSize,'filled','MarkerFaceColor',temptemp_color,...
                %   'MarkerFaceAlpha',temp_MarkerFaceAlpha,'MarkerEdgeColor',temptemp_color);%[0 0 0]
                
                scatter(x(tempi,:),y(tempi,:),tempSize,'filled','MarkerFaceColor',temptemp_color,...
                    'MarkerFaceAlpha',temp_MarkerFaceAlpha);%[0 0 0]
                hold on
                
                scatter(x_mean(tempi),y_mean(tempi),tempSize,'filled','MarkerFaceColor',temptemp_color,...
                    'MarkerFaceAlpha',temp_MarkerFaceAlpha,'MarkerEdgeColor',[0 0 0]);%[0 0 0]
                hold on
                
                
                temp_h = [temp_h scatter(100,100,tempSize,'filled','MarkerFaceColor',temptemp_color,...
                    'MarkerFaceAlpha',1)];
                hold on
            end
            
            le = legend(temp_h,'WM','Meta-WM','BSL','Location','southwest','fontsize',6.5,'NumColumns',3);%10
            le.ItemTokenSize = ones(1,3)*4;%5
            legend('boxoff');
            
            %axis equal
            
            set(gca,'linewidth',1.5)
            xlim([x_min-(x_max-x_min)*0.15 x_max+(x_max-x_min)*0.15]);
            ylim([y_min-(y_max-y_min)*0.25 y_max+(y_max-y_min)*0.15]);
            
            set(gca, 'FontSize', 8)
            set(gca,'box','off');% 取消右、上边框
            %xlabel('PC1', 'FontSize', 9, 'FontWeight', 'normal');
            %ylabel('PC2', 'FontSize', 9, 'FontWeight', 'normal');
            xlabel(sprintf('PC1 (%.0f%%)',explained3(1)), 'FontSize', 9, 'FontWeight', 'normal');
            ylabel(sprintf('PC2 (%.0f%%)',explained3(2)), 'FontSize', 9, 'FontWeight', 'normal');
            %temp_title = title(sprintf('r=%.3f, %s',temp_r,tempTxt), 'FontSize', 9, 'FontWeight', 'normal');
            %temp_title.Interpreter = 'none';
        end
        
        
        %% Plot Within groups VS. Between groups' pc_distance distribution (lesionProb)
        if true
            fig = figure('Name','Centroid distance (lesionProb)','NumberTitle','off');
            %set(gcf,'Position',[700 45 240*0.80 336*1.11*0.67*0.9*1.2*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[700 45 240*0.80*0.93 336*1.11*0.67*0.9*1.2*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
            nexttile
            
            temp_1 = temp_pcDisNew_1to1_mean';
            temp_2 = temp_pcDisNew_2to2_mean';
            temp_3 = temp_pcDisNew_3to3_mean';
            
            temp_x1_mean = pcCentriodDis_1toOthers;
            temp_x2_mean = pcCentriodDis_2toOthers;
            temp_x3_mean = pcCentriodDis_3toOthers;
            
            
            temp_p_1 = temp_p_11_1Others_pc;
            temp_p_2 = temp_p_22_2Others_pc;
            temp_p_3 = temp_p_33_3Others_pc;
            
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
            
            
            plot([1-0.45 1+0.45],temp_x1_mean*[1 1],'color', [0.25 0.25 0.25],'linewidth',4);%[0.6350 0.0780 0.1840]
            hold on
            
            plot([2-0.45 2+0.45],temp_x2_mean*[1 1],'color', [0.25 0.25 0.25],'linewidth',4);
            hold on
            
            plot([3-0.45 3+0.45],temp_x3_mean*[1 1],'color', [0.25 0.25 0.25],'linewidth',4);
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
            
            
            set(gca,'linewidth',1.5)
            
            xlim([0.3 3.7]);
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
            set(gca, 'FontSize', 8) %14
            set(gca,'box','off');% 取消右、上边框
            
            
            xtl = ["WM", "Meta", "BSL"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.33;
            text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',8);%9-->12
            
            set(gca,'xticklabel','');
            
            %yticks([0 40 80]);
            
            ylabel('Distance (PC)', 'FontSize', 9, 'FontWeight', 'normal');
            %temp_title = title(sprintf('Centroid distance'),'fontsize',9);
            %temp_title.Interpreter = 'none';
            
        end
        
        
        %% PCA along neuron dimension (eachFOV)
        if if_PCA_singleFOV0_multiFOV1 == 1
            if if_plot == 1 && if_plot_part6 == 1
                fig = figure('Name','PCA along neuron dimension (eachFOV)','NumberTitle','off');
                %set(gcf,'Position',[900 200 240*0.85*0.95 240*1.05*0.67*1.05*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                set(gcf,'Position',[900 200 240*0.85*0.95*0.95*0.93 240*1.05*0.67*1.05*0.95*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                
                t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
                
                
                nexttile
                
                x = squeeze(temp_score3_eachFOV(:,1,:));
                y = squeeze(temp_score3_eachFOV(:,2,:));
                
                x_mean = temp_score3_eachFOV_mean(:,1);
                y_mean = temp_score3_eachFOV_mean(:,2);
                
                
                
                [x_min,x_max] = bounds(x(:));
                [y_min,y_max] = bounds(y(:));
                
                temp_h = [];
                
                tempSize = 16;%30,8
                temp_MarkerFaceAlpha = 0.7;%1,0.4
                for tempi=1:3
                    if tempi == 1
                        temptemp_color = color_memoryQuality;
                    elseif tempi == 2
                        temptemp_color = color_meta;
                    elseif tempi == 3
                        temptemp_color = color_prior;
                    end
                    
                    %scatter(x(tempi,:),y(tempi,:),tempSize,'filled','MarkerFaceColor',temptemp_color,...
                    %   'MarkerFaceAlpha',temp_MarkerFaceAlpha,'MarkerEdgeColor',temptemp_color);%[0 0 0]
                    
                    scatter(x(tempi,:),y(tempi,:),tempSize,'filled','MarkerFaceColor',temptemp_color,...
                        'MarkerFaceAlpha',temp_MarkerFaceAlpha);%[0 0 0]
                    hold on
                    
                    scatter(x_mean(tempi),y_mean(tempi),tempSize,'filled','MarkerFaceColor',temptemp_color,...
                        'MarkerFaceAlpha',temp_MarkerFaceAlpha,'MarkerEdgeColor',[0 0 0]);%[0 0 0]
                    hold on
                    
                    temp_h = [temp_h scatter(100,100,tempSize,'filled','MarkerFaceColor',temptemp_color,...
                        'MarkerFaceAlpha',1)];
                    hold on
                    
                end
                
                le = legend(temp_h,'WM','Meta-WM','BSL','Location','southwest','fontsize',6.5,'NumColumns',3);%10
                le.ItemTokenSize = ones(1,3)*4;%5
                legend('boxoff');
                
                
                %axis equal
                
                set(gca,'linewidth',1.5)
                xlim([x_min-(x_max-x_min)*0.15 x_max+(x_max-x_min)*0.15]);
                ylim([y_min-(y_max-y_min)*0.30 y_max+(y_max-y_min)*0.15]);
                
                set(gca, 'FontSize', 8)
                set(gca,'box','off');% 取消右、上边框
                %xlabel('PC1', 'FontSize', 9, 'FontWeight', 'normal');
                %ylabel('PC2', 'FontSize', 9, 'FontWeight', 'normal');
                xlabel(sprintf('PC1 (%.0f%%)',explained3(1)), 'FontSize', 9, 'FontWeight', 'normal');
                ylabel(sprintf('PC2 (%.0f%%)',explained3(2)), 'FontSize', 9, 'FontWeight', 'normal');
                %temp_title = title(sprintf('r=%.3f, %s',temp_r,tempTxt), 'FontSize', 9, 'FontWeight', 'normal');
                %temp_title.Interpreter = 'none';
            end
        end
        
        
        %% Plot Within groups VS. Between groups' pc_distance distribution (eachFOV)
        if if_PCA_singleFOV0_multiFOV1 == 1
            if if_plot == 1 && if_plot_part6 == 1
                
                fig = figure('Name','Centroid distance (eachFOV)','NumberTitle','off');
                %set(gcf,'Position',[1100 45 240*0.80 336*1.11*0.67*0.9*1.2*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                set(gcf,'Position',[1100 45 240*0.80*0.93 336*1.11*0.67*0.9*1.2*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                
                t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
                nexttile
                
                temp_1 = temp_pcDisNewFOV_1to1_mean';
                temp_2 = temp_pcDisNewFOV_2to2_mean';
                temp_3 = temp_pcDisNewFOV_3to3_mean';
                
                temp_x1_mean = pcCentriodDisFOV_1toOthers;
                temp_x2_mean = pcCentriodDisFOV_2toOthers;
                temp_x3_mean = pcCentriodDisFOV_3toOthers;
                
                
                temp_p_1 = temp_p_11_1Others_pcFOV;
                temp_p_2 = temp_p_22_2Others_pcFOV;
                temp_p_3 = temp_p_33_3Others_pcFOV;
                
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
                
                
                plot([1-0.45 1+0.45],temp_x1_mean*[1 1],'color', [0.25 0.25 0.25],'linewidth',4);%[0.6350 0.0780 0.1840]
                hold on
                
                plot([2-0.45 2+0.45],temp_x2_mean*[1 1],'color', [0.25 0.25 0.25],'linewidth',4);
                hold on
                
                plot([3-0.45 3+0.45],temp_x3_mean*[1 1],'color', [0.25 0.25 0.25],'linewidth',4);
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
                
                
                set(gca,'linewidth',1.5)
                
                xlim([0.3 3.7]);
                ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
                set(gca, 'FontSize', 8) %14
                set(gca,'box','off');% 取消右、上边框
                
                
                xtl = ["WM", "Meta", "BSL"];
                xt=get(gca,'XTick');
                yt=get(gca,'YTick');
                xtext_xp=xt;
                xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.33;
                text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',0,'fontsize',8);%9-->12
                
                set(gca,'xticklabel','');
                
                %yticks([0 40 80]);
                
                ylabel('Distance (PC)', 'FontSize', 9, 'FontWeight', 'normal');
                %temp_title = title(sprintf('Centroid distance'),'fontsize',9);
                %temp_title.Interpreter = 'none';
                
            end
        end
        
        
        %% Angles between vecters of neuron tuning strength (resampling)
        if if_plot == 1 && if_plot_part6 == 1
            fig = figure('Name','Vecter angle (resampling)','NumberTitle','off');
            %set(gcf,'Position',[50+0 400+0 240*0.80 336*1.11*0.67*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[950+0 450+0 240*0.80*1.03 336*1.11*0.67*0.9*1.28*1.06]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
            nexttile
            
            temp_1 = temp_degree_ab_resampled;
            temp_2 = temp_degree_ac_resampled;
            temp_3 = temp_degree_bc_resampled;
            
            temp_chanceLevel = 0;
            
            temp_p_1 = temp_p_ab_resampled;
            temp_p_2 = temp_p_ac_resampled;
            temp_p_3 = temp_p_bc_resampled;
            
            if false
                temp_y_min = min([temp_1 temp_2 temp_3 temp_chanceLevel]);
                temp_y_max = max([temp_1 temp_2 temp_3 temp_chanceLevel]);
                
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
                
            end
            
            %temp1_SEM = std(temp_1)/sqrt(length(temp_1));
            %temp2_SEM = std(temp_2)/sqrt(length(temp_2));
            temp1_SEM = std(temp_1);
            temp2_SEM = std(temp_2);
            temp3_SEM = std(temp_3);
            
            temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM,temp_chanceLevel]);
            temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM,temp_chanceLevel]);
            
            %         bar([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)], ...
            %             'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',0.5);
            %         hold on
            plot([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],'color',[0 0 0],'linewidth',0.5);
            hold on
            errorbar([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM],...
                '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);
            hold on
            
            
            %x = [1 2 3];
            %y = [mean(temp_1) mean(temp_2) mean(temp_3)];
            %y_sem = [temp1_SEM temp2_SEM temp3_SEM];
            %patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0 0 0],'FaceAlpha',0.5,'EdgeColor',[0 0 0]);
            %hold on
            
            
            tempTxt = sprintf('');
            if temp_p_1 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p_1 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p_1 < 0.05
                tempTxt = sprintf('*');
            end
            text(1,temp_y_min+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            
            tempTxt = sprintf('');
            if temp_p_2 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p_2 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p_2 < 0.05
                tempTxt = sprintf('*');
            end
            text(2,temp_y_min+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            
            tempTxt = sprintf('');
            if temp_p_3 < 0.001
                tempTxt = sprintf('***');
            elseif temp_p_3 < 0.01
                tempTxt = sprintf('**');
            elseif temp_p_3 < 0.05
                tempTxt = sprintf('*');
            end
            text(3,temp_y_min+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            plot([0.3 3.7],[1 1].*temp_chanceLevel,'--','LineWidth',1,'color',[0.5 0.5 0.5]);
            hold on
            
            
            set(gca,'linewidth',1.5)
            
            xlim([0.3 3.7]);
            %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
            set(gca, 'FontSize', 8) %14
            set(gca,'box','off');% 取消右、上边框
            
            
            xtl = ["WM-Meta", "WM-BSL", "Meta-BSL"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.32;
            xtext_yp=temp_y_min*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.14;
            text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
            
            set(gca,'xticklabel','');
            
            %yticks([0 40 80]);
            
            ylabel('Angle (degree)', 'FontSize', 9, 'FontWeight', 'normal');
            temp_title = title(sprintf('Tuning strength vectors'),'fontsize',9);
            temp_title.Interpreter = 'none';
            
            
        end
        
        
        %% Angles between vecters of neuron tuning strength (eachFOV)
        
        if if_PCA_singleFOV0_multiFOV1 == 1
            if if_plot == 1 && if_plot_part6 == 1
                
                fig = figure('Name','Vecter angle (eachFOV)','NumberTitle','off');
                %set(gcf,'Position',[50+0 400+0 240*0.80 336*1.11*0.67*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                %set(gcf,'Position',[1175+0 450+0 240*0.80*1.03 336*1.11*0.67*0.9*1.28*1.06]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                set(gcf,'Position',[1175+0 450+0 240*0.80*1.03*0.95 336*1.11*0.67*0.9*1.28*1.06]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                
                t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
                nexttile
                
                temp_1 = temp_degree_ab_eachFOV;
                temp_2 = temp_degree_ac_eachFOV;
                temp_3 = temp_degree_bc_eachFOV;
                
                temp_chanceLevel = 0;
                
                temp_p_1 = temp_p_ab_eachFOV;
                temp_p_2 = temp_p_ac_eachFOV;
                temp_p_3 = temp_p_bc_eachFOV;
                
                if false
                    temp_y_min = min([temp_1 temp_2 temp_3 temp_chanceLevel]);
                    temp_y_max = max([temp_1 temp_2 temp_3 temp_chanceLevel]);
                    
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
                    
                end
                
                %         temp1_SEM = std(temp_1)/sqrt(length(temp_1));
                %         temp2_SEM = std(temp_2)/sqrt(length(temp_2));
                temp1_SEM = std(temp_1);
                temp2_SEM = std(temp_2);
                temp3_SEM = std(temp_3);
                
                temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM,temp_chanceLevel]);
                temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM,temp_chanceLevel]);
                
                %         bar([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)], ...
                %             'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',0.5);
                %         hold on
                plot([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],'color',[0 0 0],'linewidth',0.5);
                hold on
                errorbar([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM],...
                    '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);
                hold on
                
                
                %x = [1 2 3];
                %y = [mean(temp_1) mean(temp_2) mean(temp_3)];
                %y_sem = [temp1_SEM temp2_SEM temp3_SEM];
                %patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0 0 0],'FaceAlpha',0.5,'EdgeColor',[0 0 0]);
                %hold on
                
                
                tempTxt = sprintf('');
                if temp_p_1 < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p_1 < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p_1 < 0.05
                    tempTxt = sprintf('*');
                end
                text(1,temp_y_min+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                
                tempTxt = sprintf('');
                if temp_p_2 < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p_2 < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p_2 < 0.05
                    tempTxt = sprintf('*');
                end
                text(2,temp_y_min+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                
                tempTxt = sprintf('');
                if temp_p_3 < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p_3 < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p_3 < 0.05
                    tempTxt = sprintf('*');
                end
                text(3,temp_y_min+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                plot([0.3 3.7],[1 1].*temp_chanceLevel,'--','LineWidth',1,'color',[0.5 0.5 0.5]);
                hold on
                
                
                set(gca,'linewidth',1.5)
                
                xlim([0.3 3.7]);
                %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
                ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
                set(gca, 'FontSize', 8) %14
                set(gca,'box','off');% 取消右、上边框
                
                
                xtl = ["WM-Meta", "WM-BSL", "Meta-BSL"];
                xt=get(gca,'XTick');
                yt=get(gca,'YTick');
                xtext_xp=xt;
                %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.32;
                xtext_yp=temp_y_min*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.14;
                text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
                
                set(gca,'xticklabel','');
                
                %yticks([0 40 80]);
                
                ylabel('Angle (degree)', 'FontSize', 9, 'FontWeight', 'normal');
                temp_title = title(sprintf('Tuning strength vectors'),'fontsize',9);
                temp_title.Interpreter = 'none';
            end
            
        end
        
        
        
        %% Angles between vecters of neuron tuning strength (eachFOV)
        if if_PCA_singleFOV0_multiFOV1 == 1
            if if_plot == 1 && if_plot_part6 == 1
                
                fig = figure('Name','Vecter angle (eachFOV)','NumberTitle','off');
                %set(gcf,'Position',[50+0 400+0 240*0.80 336*1.11*0.67*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                %set(gcf,'Position',[1175+0 450+0 240*0.80*1.03 336*1.11*0.67*0.9*1.28*1.06]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                set(gcf,'Position',[1175+0 450+0 240*0.80*1.03*0.95 336*1.11*0.67*0.9*1.28*1.06]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
                
                t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
                nexttile
                
                temp_3 = temp_degree_ab_eachFOV;
                temp_1 = temp_degree_ac_eachFOV;
                temp_2 = temp_degree_bc_eachFOV;
                
                temp_chanceLevel = 0;
                
                temp_p_3 = temp_p_ab_eachFOV;
                temp_p_1 = temp_p_ac_eachFOV;
                temp_p_2 = temp_p_bc_eachFOV;
                
                if false
                    temp_y_min = min([temp_1 temp_2 temp_3 temp_chanceLevel]);
                    temp_y_max = max([temp_1 temp_2 temp_3 temp_chanceLevel]);
                    
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
                    
                end
                
                temp1_SEM = std(temp_1)/sqrt(length(temp_1));
                temp2_SEM = std(temp_2)/sqrt(length(temp_2));
                temp3_SEM = std(temp_3)/sqrt(length(temp_3));
                %temp1_SEM = std(temp_1);
                %temp2_SEM = std(temp_2);
                %temp3_SEM = std(temp_3);
                
                temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM,temp_chanceLevel]);
                temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM,temp_chanceLevel]);
                
                %         bar([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)], ...
                %             'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',0.5);
                %         hold on
                plot([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],'color',[0 0 0],'linewidth',0.5);
                hold on
                errorbar([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM],...
                    '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);
                hold on
                
                
                %x = [1 2 3];
                %y = [mean(temp_1) mean(temp_2) mean(temp_3)];
                %y_sem = [temp1_SEM temp2_SEM temp3_SEM];
                %patch([x(:); flipud(x(:))]', [y(:)+y_sem(:); flipud(y(:)-y_sem(:))]',[0 0 0],'FaceAlpha',0.5,'EdgeColor',[0 0 0]);
                %hold on
                
                
                tempTxt = sprintf('');
                if temp_p_1 < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p_1 < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p_1 < 0.05
                    tempTxt = sprintf('*');
                end
                text(1,temp_y_min+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                
                tempTxt = sprintf('');
                if temp_p_2 < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p_2 < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p_2 < 0.05
                    tempTxt = sprintf('*');
                end
                text(2,temp_y_min+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                
                tempTxt = sprintf('');
                if temp_p_3 < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p_3 < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p_3 < 0.05
                    tempTxt = sprintf('*');
                end
                text(3,temp_y_min+(temp_y_max-temp_y_min)*0.05,tempTxt,'Color',[0.25 0.25 0.25],'FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                plot([0.3 3.7],[1 1].*temp_chanceLevel,'--','LineWidth',1,'color',[0.5 0.5 0.5]);
                hold on
                
                
                set(gca,'linewidth',1.5)
                
                xlim([0.3 3.7]);
                %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
                ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
                set(gca, 'FontSize', 8) %14
                set(gca,'box','off');% 取消右、上边框
                
                
                %xtl = ["WM-Meta", "WM-BSL", "Meta-BSL"];
                xtl = ["BSL-WM", "BSL-Meta","WM-Meta"];
                xt=get(gca,'XTick');
                yt=get(gca,'YTick');
                xtext_xp=xt;
                %xtext_yp=(yt(1))*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.32;
                xtext_yp=temp_y_min*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.14;
                text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
                
                set(gca,'xticklabel','');
                
                %yticks([0 40 80]);
                
                ylabel('Angle (degree)', 'FontSize', 9, 'FontWeight', 'normal');
                temp_title = title(sprintf('Tuning strength vectors'),'fontsize',9);
                temp_title.Interpreter = 'none';
            end
            
        end
        
        
    end
end


%% Unit view of tuning strength
if if_plot == 1 & if_plot_part6 == 1 %#ok<*AND2>
% if true
    temp_degree_aa_prior_unit_summary = allPart_summary_mulFOV.part6_summary_mulFOV.betaAngle_choiceMemory_baseline_halfAB_norm_unitView_summary;
    temp_degree_bb_memory_unit_summary = allPart_summary_mulFOV.part6_summary_mulFOV.betaAngle_6loc_halfAB_norm_unitView_summary;
    temp_degree_cc_meta_unit_summary = allPart_summary_mulFOV.part6_summary_mulFOV.betaAngle_gProb_halfAB_norm_unitView_summary;
    
    temp_degree_ab_priorMemory_unit_summary = temp_degree_ac_eachFOV';
    temp_degree_ac_priorMeta_unit_summary = temp_degree_bc_eachFOV';
    temp_degree_bc_memoryMeta_unit_summary = temp_degree_ab_eachFOV';
    
    
    
    temp_degree_aa_prior_unit_summary_mean = nan(size(temp_degree_ab_priorMemory_unit_summary));
    temp_degree_bb_memory_unit_summary_mean = temp_degree_aa_prior_unit_summary_mean;
    temp_degree_cc_meta_unit_summary_mean = temp_degree_aa_prior_unit_summary_mean;
    
    temp_FOVNum = length(temp_degree_aa_prior_unit_summary_mean);
    
    temp_repeatNum = round(length(temp_degree_aa_prior_unit_summary)./temp_FOVNum);
    
    for tempi=1:temp_FOVNum
        temp_range = ((tempi-1)*temp_repeatNum+1):tempi*temp_repeatNum;
        
        temp_degree_aa_prior_unit_summary_mean(tempi) = mean(temp_degree_aa_prior_unit_summary(temp_range));
        temp_degree_bb_memory_unit_summary_mean(tempi) = mean(temp_degree_bb_memory_unit_summary(temp_range));
        temp_degree_cc_meta_unit_summary_mean(tempi) = mean(temp_degree_cc_meta_unit_summary(temp_range));
    end
    
    temp_degree_aa_prior_unit_summary_mean;
    temp_degree_bb_memory_unit_summary_mean;
    temp_degree_cc_meta_unit_summary_mean;
    
    temp_degree_ab_priorMemory_unit_summary;
    temp_degree_ac_priorMeta_unit_summary;
    temp_degree_bc_memoryMeta_unit_summary;
    
    
    %     [~,temp_p_aa_ab_unit] = ttest(temp_degree_aa_prior_unit_summary_mean,temp_degree_ab_priorMemory_unit_summary);
    %     [~,temp_p_aa_ac_unit] = ttest(temp_degree_aa_prior_unit_summary_mean,temp_degree_ac_priorMeta_unit_summary);
    %
    %     [~,temp_p_bb_ab_unit] = ttest(temp_degree_bb_memory_unit_summary_mean,temp_degree_ab_priorMemory_unit_summary);
    %     [~,temp_p_bb_bc_unit] = ttest(temp_degree_bb_memory_unit_summary_mean,temp_degree_bc_memoryMeta_unit_summary);
    %
    %     [~,temp_p_cc_ac_unit] = ttest(temp_degree_cc_meta_unit_summary_mean,temp_degree_ac_priorMeta_unit_summary);
    %     [~,temp_p_cc_bc_unit] = ttest(temp_degree_cc_meta_unit_summary_mean,temp_degree_bc_memoryMeta_unit_summary);
    
    
    [~,temp_p_aa_ab_unit] = ttest2(temp_degree_aa_prior_unit_summary,temp_degree_ab_priorMemory_unit_summary);
    [~,temp_p_aa_ac_unit] = ttest2(temp_degree_aa_prior_unit_summary,temp_degree_ac_priorMeta_unit_summary);
    
    [~,temp_p_bb_ab_unit] = ttest2(temp_degree_bb_memory_unit_summary,temp_degree_ab_priorMemory_unit_summary);
    [~,temp_p_bb_bc_unit] = ttest2(temp_degree_bb_memory_unit_summary,temp_degree_bc_memoryMeta_unit_summary);
    
    [~,temp_p_cc_ac_unit] = ttest2(temp_degree_cc_meta_unit_summary,temp_degree_ac_priorMeta_unit_summary);
    [~,temp_p_cc_bc_unit] = ttest2(temp_degree_cc_meta_unit_summary,temp_degree_bc_memoryMeta_unit_summary);
    
    fprintf('Degree between beta vectors (unit view): \nwithin (if_resample_replace=0) [%.1f, %.1f, %.1f], between [%.1f, %.1f, %.1f].\n',...
        mean(temp_degree_aa_prior_unit_summary_mean),mean(temp_degree_bb_memory_unit_summary_mean),mean(temp_degree_cc_meta_unit_summary_mean),...
        mean(temp_degree_ab_priorMemory_unit_summary),mean(temp_degree_ac_priorMeta_unit_summary),mean(temp_degree_bc_memoryMeta_unit_summary));
    
end

%% Plot unit view of tuning strength
if if_plot == 1 && if_plot_part6 == 1
    
    fig = figure('Name','Vecter angle (eachFOV)','NumberTitle','off');
    set(gcf,'Position',[1375+0 450+0 240*0.80*1.03*0.95*1.6*0.8 336*1.11*0.67*0.9*1.28*1.06]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_1 = temp_degree_ab_priorMemory_unit_summary;
    temp_2 = temp_degree_ac_priorMeta_unit_summary;
    temp_3 = temp_degree_bc_memoryMeta_unit_summary;
    
    %     temp_1B = temp_degree_aa_prior_unit_summary_mean;
    %     temp_2B = temp_degree_bb_memory_unit_summary_mean;
    %     temp_3B = temp_degree_cc_meta_unit_summary_mean;
    temp_1B = temp_degree_aa_prior_unit_summary;
    temp_2B = temp_degree_bb_memory_unit_summary;
    temp_3B = temp_degree_cc_meta_unit_summary;
    
    
    %     temp1_SEM = std(temp_1);
    %     temp2_SEM = std(temp_2);
    %     temp3_SEM = std(temp_3);
    %
    %     temp1B_SEM = std(temp_1B);
    %     temp2B_SEM = std(temp_2B);
    %     temp3B_SEM = std(temp_3B);
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    temp3_SEM = std(temp_3)/sqrt(length(temp_3));
    
    temp1B_SEM = std(temp_1B)/sqrt(length(temp_1B));
    temp2B_SEM = std(temp_2B)/sqrt(length(temp_2B));
    temp3B_SEM = std(temp_3B)/sqrt(length(temp_3B));
    
    
    
    temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM,...
        mean(temp_1B)-temp1B_SEM,mean(temp_2B)-temp2B_SEM,mean(temp_3B)-temp3B_SEM]);
    temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM,...
        mean(temp_1B)+temp1B_SEM,mean(temp_2B)+temp2B_SEM,mean(temp_3B)+temp3B_SEM]);
    
    plot([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],'color',[0 0 0],'linewidth',0.5);
    hold on
    errorbar([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM],...
        '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);
    hold on
    
    plot([4 5 6], [mean(temp_1B) mean(temp_2B) mean(temp_3B)],'color',[0 0 0],'linewidth',0.5);
    hold on
    errorbar([4 5 6], [mean(temp_1B) mean(temp_2B) mean(temp_3B)],[temp1B_SEM temp2B_SEM temp3B_SEM],...
        '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);
    hold on
    
    
    
    
    
    set(gca,'linewidth',1.5)
    
    %xlim([0.3 6.7]);
    xlim([0.5 6.5]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
    set(gca, 'FontSize', 8) %14
    set(gca,'box','off');% 取消右、上边框
    
    
    xtl = ["BSL-WM", "BSL-Meta","WM-Meta","BSL-BSL", "WM-WM","Meta-Meta"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=temp_y_min*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.14;
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.14;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('Angle (degree)', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('Tuning strength vectors'),'fontsize',9);
    temp_title.Interpreter = 'none';
end


%% Neuronal linear regression, using neuron activity to estimate task variables linearly.
% if if_plot == 1 & if_plot_part6 == 1 %#ok<*AND2>
if true
    beta_prior_linearAxis_summary = allPart_summary_mulFOV.part6_summary_mulFOV.beta_prior_linearAxis_summary;
    beta_memory_linearAxis_summary = allPart_summary_mulFOV.part6_summary_mulFOV.beta_memory_linearAxis_summary;
    beta_meta_linearAxis_summary = allPart_summary_mulFOV.part6_summary_mulFOV.beta_meta_linearAxis_summary;
    
    temp_degree_aa_prior_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_degree_aa_prior_summary;
    temp_degree_bb_memory_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_degree_bb_memory_summary;
    temp_degree_cc_meta_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_degree_cc_meta_summary;
    
    temp_degree_ab_priorMemory_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_degree_ab_priorMemory_summary;
    temp_degree_ac_priorMeta_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_degree_ac_priorMeta_summary;
    temp_degree_bc_memoryMeta_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_degree_bc_memoryMeta_summary;
    
    temp_vectors_linearAxix = [beta_prior_linearAxis_summary,beta_memory_linearAxis_summary,beta_meta_linearAxis_summary];
    
    
    temp_degree_aa_prior_summary_mean = nan(size(temp_degree_ab_priorMemory_summary));
    temp_degree_bb_memory_summary_mean = temp_degree_aa_prior_summary_mean;
    temp_degree_cc_meta_summary_mean = temp_degree_aa_prior_summary_mean;
    
    temp_FOVNum = length(temp_degree_aa_prior_summary_mean);
    
    temp_repeatNum = round(length(temp_degree_aa_prior_summary)./temp_FOVNum);
    
    for tempi=1:temp_FOVNum
        temp_range = ((tempi-1)*temp_repeatNum+1):tempi*temp_repeatNum;
        
        temp_degree_aa_prior_summary_mean(tempi) = mean(temp_degree_aa_prior_summary(temp_range));
        temp_degree_bb_memory_summary_mean(tempi) = mean(temp_degree_bb_memory_summary(temp_range));
        temp_degree_cc_meta_summary_mean(tempi) = mean(temp_degree_cc_meta_summary(temp_range));
    end
    
    temp_degree_aa_prior_summary_mean;
    temp_degree_bb_memory_summary_mean;
    temp_degree_cc_meta_summary_mean;
    
    temp_degree_ab_priorMemory_summary;
    temp_degree_ac_priorMeta_summary;
    temp_degree_bc_memoryMeta_summary;
    
    
    %     [~,temp_p_aa_ab] = ttest(temp_degree_aa_prior_summary_mean,temp_degree_ab_priorMemory_summary);
    %     [~,temp_p_aa_ac] = ttest(temp_degree_aa_prior_summary_mean,temp_degree_ac_priorMeta_summary);
    %
    %     [~,temp_p_bb_ab] = ttest(temp_degree_bb_memory_summary_mean,temp_degree_ab_priorMemory_summary);
    %     [~,temp_p_bb_bc] = ttest(temp_degree_bb_memory_summary_mean,temp_degree_bc_memoryMeta_summary);
    %
    %     [~,temp_p_cc_ac] = ttest(temp_degree_cc_meta_summary_mean,temp_degree_ac_priorMeta_summary);
    %     [~,temp_p_cc_bc] = ttest(temp_degree_cc_meta_summary_mean,temp_degree_bc_memoryMeta_summary);
    
    
    [~,temp_p_aa_ab] = ttest2(temp_degree_aa_prior_summary,temp_degree_ab_priorMemory_summary);
    [~,temp_p_aa_ac] = ttest2(temp_degree_aa_prior_summary,temp_degree_ac_priorMeta_summary);
    
    [~,temp_p_bb_ab] = ttest2(temp_degree_bb_memory_summary,temp_degree_ab_priorMemory_summary);
    [~,temp_p_bb_bc] = ttest2(temp_degree_bb_memory_summary,temp_degree_bc_memoryMeta_summary);
    
    [~,temp_p_cc_ac] = ttest2(temp_degree_cc_meta_summary,temp_degree_ac_priorMeta_summary);
    [~,temp_p_cc_bc] = ttest2(temp_degree_cc_meta_summary,temp_degree_bc_memoryMeta_summary);
    
    fprintf('Degree between beta vectors: \nwithin (if_resample_replace=0) [%.1f, %.1f, %.1f], between [%.1f, %.1f, %.1f].\n',...
        mean(temp_degree_aa_prior_summary_mean),mean(temp_degree_bb_memory_summary_mean),mean(temp_degree_cc_meta_summary_mean),...
        mean(temp_degree_ab_priorMemory_summary),mean(temp_degree_ac_priorMeta_summary),mean(temp_degree_bc_memoryMeta_summary));
    
    
    
    
    temp_VAF_ratio_aa_prior_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_VAF_ratio_priorPrior_summary;
    temp_VAF_ratio_bb_memory_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_VAF_ratio_memoryMemory_summary;
    temp_VAF_ratio_cc_meta_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_VAF_ratio_metaMeta_summary;
    
    temp_VAF_ratio_ab_priorMemory_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_VAF_ratio_priorMemory_summary;
    temp_VAF_ratio_ac_priorMeta_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_VAF_ratio_priorMeta_summary;
    temp_VAF_ratio_bc_memoryMeta_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_VAF_ratio_memoryMeta_summary;
    
    
    temp_VAF_ratio_aa_prior_summary_mean = nan(size(temp_VAF_ratio_ab_priorMemory_summary));
    temp_VAF_ratio_bb_memory_summary_mean = temp_VAF_ratio_aa_prior_summary_mean;
    temp_VAF_ratio_cc_meta_summary_mean = temp_VAF_ratio_aa_prior_summary_mean;
    
    temp_FOVNum = length(temp_VAF_ratio_aa_prior_summary_mean);
    
    temp_repeatNum = round(length(temp_VAF_ratio_aa_prior_summary)./temp_FOVNum);
    
    for tempi=1:temp_FOVNum
        temp_range = ((tempi-1)*temp_repeatNum+1):tempi*temp_repeatNum;
        
        temp_VAF_ratio_aa_prior_summary_mean(tempi) = mean(temp_VAF_ratio_aa_prior_summary(temp_range));
        temp_VAF_ratio_bb_memory_summary_mean(tempi) = mean(temp_VAF_ratio_bb_memory_summary(temp_range));
        temp_VAF_ratio_cc_meta_summary_mean(tempi) = mean(temp_VAF_ratio_cc_meta_summary(temp_range));
    end
    
    
    %     [~,temp_p_aa_ab_VAF] = ttest2(temp_VAF_ratio_aa_prior_summary,temp_VAF_ratio_ab_priorMemory_summary);
    %     [~,temp_p_aa_ac_VAF] = ttest2(temp_VAF_ratio_aa_prior_summary,temp_VAF_ratio_ac_priorMeta_summary);
    %
    %     [~,temp_p_bb_ab_VAF] = ttest2(temp_VAF_ratio_bb_memory_summary,temp_VAF_ratio_ab_priorMemory_summary);
    %     [~,temp_p_bb_bc_VAF] = ttest2(temp_VAF_ratio_bb_memory_summary,temp_VAF_ratio_bc_memoryMeta_summary);
    %
    %     [~,temp_p_cc_ac_VAF] = ttest2(temp_VAF_ratio_cc_meta_summary,temp_VAF_ratio_ac_priorMeta_summary);
    %     [~,temp_p_cc_bc_VAF] = ttest2(temp_VAF_ratio_cc_meta_summary,temp_VAF_ratio_bc_memoryMeta_summary);
    
    [~,temp_p_aa_ab_VAF] = ttest(temp_VAF_ratio_aa_prior_summary_mean,temp_VAF_ratio_ab_priorMemory_summary);
    [~,temp_p_aa_ac_VAF] = ttest(temp_VAF_ratio_aa_prior_summary_mean,temp_VAF_ratio_ac_priorMeta_summary);
    
    [~,temp_p_bb_ab_VAF] = ttest(temp_VAF_ratio_bb_memory_summary_mean,temp_VAF_ratio_ab_priorMemory_summary);
    [~,temp_p_bb_bc_VAF] = ttest(temp_VAF_ratio_bb_memory_summary_mean,temp_VAF_ratio_bc_memoryMeta_summary);
    
    [~,temp_p_cc_ac_VAF] = ttest(temp_VAF_ratio_cc_meta_summary_mean,temp_VAF_ratio_ac_priorMeta_summary);
    [~,temp_p_cc_bc_VAF] = ttest(temp_VAF_ratio_cc_meta_summary_mean,temp_VAF_ratio_bc_memoryMeta_summary);
    
    fprintf('VAF ratio between beta vectors (subspace view): within [%.2f, %.2f, %.2f], between [%.2f, %.2f, %.2f].\n\n',...
        mean(temp_VAF_ratio_aa_prior_summary),mean(temp_VAF_ratio_bb_memory_summary),mean(temp_VAF_ratio_cc_meta_summary),...
        mean(temp_VAF_ratio_ab_priorMemory_summary),mean(temp_VAF_ratio_ac_priorMeta_summary),mean(temp_VAF_ratio_bc_memoryMeta_summary));
    
end



%% (Angle) Plot Neuronal linear regression, using neuron activity to estimate task variables linearly.
if if_plot == 1 && if_plot_part6 == 1
    
    fig = figure('Name','Vecter angle (eachFOV)','NumberTitle','off');
    %set(gcf,'Position',[1375+0 450+0 240*0.80*1.03*0.95*1.6*0.8 336*1.11*0.67*0.9*1.28*1.06]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[1375+0 450+0 240.5*0.8*0.95 195.3]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_1 = temp_degree_ab_priorMemory_summary;
    temp_2 = temp_degree_ac_priorMeta_summary;
    temp_3 = temp_degree_bc_memoryMeta_summary;
    
    %     temp_1B = temp_degree_aa_prior_summary_mean;
    %     temp_2B = temp_degree_bb_memory_summary_mean;
    %     temp_3B = temp_degree_cc_meta_summary_mean;
    temp_1B = temp_degree_aa_prior_summary;
    temp_2B = temp_degree_bb_memory_summary;
    temp_3B = temp_degree_cc_meta_summary;
    
    
    %     temp1_SEM = std(temp_1);
    %     temp2_SEM = std(temp_2);
    %     temp3_SEM = std(temp_3);
    %
    %     temp1B_SEM = std(temp_1B);
    %     temp2B_SEM = std(temp_2B);
    %     temp3B_SEM = std(temp_3B);
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    temp3_SEM = std(temp_3)/sqrt(length(temp_3));
    
    temp1B_SEM = std(temp_1B)/sqrt(length(temp_1B));
    temp2B_SEM = std(temp_2B)/sqrt(length(temp_2B));
    temp3B_SEM = std(temp_3B)/sqrt(length(temp_3B));
    
    
    
    temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM,...
        mean(temp_1B)-temp1B_SEM,mean(temp_2B)-temp2B_SEM,mean(temp_3B)-temp3B_SEM]);
    temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM,...
        mean(temp_1B)+temp1B_SEM,mean(temp_2B)+temp2B_SEM,mean(temp_3B)+temp3B_SEM]);
    
    plot([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],'color',[0 0 0],'linewidth',0.5);
    hold on
    errorbar([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM],...
        '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);
    hold on
    
    plot([4 5 6], [mean(temp_1B) mean(temp_2B) mean(temp_3B)],'color',[0 0 0],'linewidth',0.5);
    hold on
    errorbar([4 5 6], [mean(temp_1B) mean(temp_2B) mean(temp_3B)],[temp1B_SEM temp2B_SEM temp3B_SEM],...
        '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);
    hold on
    
    
    
    
    
    set(gca,'linewidth',1.5)
    
    %xlim([0.3 6.7]);
    xlim([0.5 6.5]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
    set(gca, 'FontSize', 8) %14
    set(gca,'box','off');% 取消右、上边框
    
    xticks([1 2 3 4 5 6]);
    
    xtl = ["BSL-WM", "BSL-Meta","WM-Meta","BSL-BSL", "WM-WM","Meta-Meta"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=temp_y_min*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.14;
    %xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.14;
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.18;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',6.5);%9-->12
    
    set(gca,'xticklabel','');
    
    %ylabel('Angle (degree)', 'FontSize', 9, 'FontWeight', 'normal');
    ylabel('Angle (deg)', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('Angles between subspaces'),'fontsize',9);
    temp_title.Interpreter = 'none';
end


%% Two monkeys. (Angle) Plot Neuronal linear regression, using neuron activity to estimate task variables linearly.
if if_plotTwoMonkey_subspace == 1
    if if_plot == 1 && if_plot_part6 == 1
        
        fig = figure('Name','Two monkeys. Vecter angle (eachFOV)','NumberTitle','off');
        set(gcf,'Position',[1375+0 450+0 240.5*0.8*0.95 195.3]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
                
        temp_1D = temp_degree_ab_priorMemory_summary_Ding;
        temp_2D = temp_degree_ac_priorMeta_summary_Ding;
        temp_3D = temp_degree_bc_memoryMeta_summary_Ding;
        temp_4D = temp_degree_aa_prior_summary_Ding;
        temp_5D = temp_degree_bb_memory_summary_Ding;
        temp_6D = temp_degree_cc_meta_summary_Ding;
        
        temp_1Z = temp_degree_ab_priorMemory_summary;
        temp_2Z = temp_degree_ac_priorMeta_summary;
        temp_3Z = temp_degree_bc_memoryMeta_summary;
        temp_4Z = temp_degree_aa_prior_summary;
        temp_5Z = temp_degree_bb_memory_summary;
        temp_6Z = temp_degree_cc_meta_summary;
        
        if true
            [~,temp_p_13D] = ttest(temp_1D,temp_3D);
            [~,temp_p_23D] = ttest(temp_2D,temp_3D);
            
            [~,temp_p_13Z] = ttest(temp_1Z,temp_3Z);
            [~,temp_p_23Z] = ttest(temp_2Z,temp_3Z);
            
            [~,temp_p_12with3D] = ttest2([temp_1D;temp_2D],temp_3D);
            [~,temp_p_12with3Z] = ttest2([temp_1Z;temp_2Z],temp_3Z);
        end
        
        temp1D_SEM = std(temp_1D)/sqrt(length(temp_1D));
        temp2D_SEM = std(temp_2D)/sqrt(length(temp_2D));
        temp3D_SEM = std(temp_3D)/sqrt(length(temp_3D));
        temp4D_SEM = std(temp_4D)/sqrt(length(temp_4D));
        temp5D_SEM = std(temp_5D)/sqrt(length(temp_5D));
        temp6D_SEM = std(temp_6D)/sqrt(length(temp_6D));
        
        temp1Z_SEM = std(temp_1Z)/sqrt(length(temp_1Z));
        temp2Z_SEM = std(temp_2Z)/sqrt(length(temp_2Z));
        temp3Z_SEM = std(temp_3Z)/sqrt(length(temp_3Z));
        temp4Z_SEM = std(temp_4Z)/sqrt(length(temp_4Z));
        temp5Z_SEM = std(temp_5Z)/sqrt(length(temp_5Z));
        temp6Z_SEM = std(temp_6Z)/sqrt(length(temp_6Z));
        
        temp_y_min123 = min([mean(temp_1D)-temp1D_SEM,mean(temp_2D)-temp2D_SEM,mean(temp_3D)-temp3D_SEM,...
            mean(temp_1Z)-temp1Z_SEM,mean(temp_2Z)-temp2Z_SEM,mean(temp_3Z)-temp3Z_SEM]);
        temp_y_max123 = max([mean(temp_1D)+temp1D_SEM,mean(temp_2D)+temp2D_SEM,mean(temp_3D)+temp3D_SEM,...
            mean(temp_1Z)+temp1Z_SEM,mean(temp_2Z)+temp2Z_SEM,mean(temp_3Z)+temp3Z_SEM]);
        
        temp_y_min456 = min([mean(temp_4D)-temp4D_SEM,mean(temp_5D)-temp5D_SEM,mean(temp_6D)-temp6D_SEM,...
            mean(temp_4Z)-temp4Z_SEM,mean(temp_5Z)-temp5Z_SEM,mean(temp_6Z)-temp6Z_SEM]);
        temp_y_max456 = max([mean(temp_4D)+temp4D_SEM,mean(temp_5D)+temp5D_SEM,mean(temp_6D)+temp6D_SEM,...
            mean(temp_4Z)+temp4Z_SEM,mean(temp_5Z)+temp5Z_SEM,mean(temp_6Z)+temp6Z_SEM]);
        
        temp_y_min = min([temp_y_min123,temp_y_min456]);
        temp_y_max = max([temp_y_max123,temp_y_max456]);
        
        temp_lineWidth_errorbar = 1;%2
        temp_capSize = 7;%10
        temp_color_D = [1 1 1]*0.5;
        temp_color_Z = [1 1 1]*0;
        
        h_line = [];
        plot([1 2 3], [mean(temp_1D) mean(temp_2D) mean(temp_3D)],'color',temp_color_D,'linewidth',0.5);
        hold on
        h_line = [h_line errorbar([1 2 3], [mean(temp_1D) mean(temp_2D) mean(temp_3D)],[temp1D_SEM temp2D_SEM temp3D_SEM],...
            '.', 'Color', temp_color_D, 'MarkerFaceColor', temp_color_D, 'LineWidth', temp_lineWidth_errorbar, 'CapSize', temp_capSize)];
        hold on
        
        plot([4 5 6], [mean(temp_4D) mean(temp_5D) mean(temp_6D)],'color',temp_color_D,'linewidth',0.5);
        hold on
        errorbar([4 5 6], [mean(temp_4D) mean(temp_5D) mean(temp_6D)],[temp4D_SEM temp5D_SEM temp6D_SEM],...
            '.', 'Color', temp_color_D, 'MarkerFaceColor', temp_color_D, 'LineWidth', temp_lineWidth_errorbar, 'CapSize', temp_capSize);
        hold on
        
        
        plot([1 2 3], [mean(temp_1Z) mean(temp_2Z) mean(temp_3Z)],'color',temp_color_Z,'linewidth',0.5);
        hold on
        h_line = [h_line errorbar([1 2 3], [mean(temp_1Z) mean(temp_2Z) mean(temp_3Z)],[temp1Z_SEM temp2Z_SEM temp3Z_SEM],...
            '.', 'Color', temp_color_Z, 'MarkerFaceColor', temp_color_Z, 'LineWidth', temp_lineWidth_errorbar, 'CapSize', temp_capSize)];
        hold on
        
        plot([4 5 6], [mean(temp_4Z) mean(temp_5Z) mean(temp_6Z)],'color',temp_color_Z,'linewidth',0.5);
        hold on
        errorbar([4 5 6], [mean(temp_4Z) mean(temp_5Z) mean(temp_6Z)],[temp4Z_SEM temp5Z_SEM temp6Z_SEM],...
            '.', 'Color', temp_color_Z, 'MarkerFaceColor', temp_color_Z, 'LineWidth', temp_lineWidth_errorbar, 'CapSize', temp_capSize);
        hold on
        
        le = legend(h_line,'Monkey D','Monkey Z','Location','northeast','fontsize',6.5);%10
        le.ItemTokenSize = ones(1,2)*8;
        legend('boxoff');
        
        set(gca,'linewidth',1.5)
        
        %xlim([0.3 6.7]);
        xlim([0.5 6.5]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.05]);%0.3
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        xticks([1 2 3 4 5 6]);
        
        xtl = ["BSL-WM", "BSL-Meta","WM-Meta","BSL-BSL", "WM-WM","Meta-Meta"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.18;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',6.5);%9-->12
        
        set(gca,'xticklabel','');
        
        %ylabel('Angle (degree)', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Angle (deg)', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('Angles between subspaces'),'fontsize',9);
        temp_title = title(sprintf('Between subspaces'),'fontsize',9);
        
    end
end


%% (VAF) Plot Neuronal linear regression, using neuron activity to estimate task variables linearly.
if if_plot == 1 && if_plot_part6 == 1
    
    fig = figure('Name','Vecter angle (eachFOV)','NumberTitle','off');
    %set(gcf,'Position',[1375+0 450+0 240*0.80*1.03*0.95*1.6*0.8 336*1.11*0.67*0.9*1.28*1.06]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[1375+0 450+0 240.5*0.8*0.95 195.3]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile
    
    temp_1 = temp_VAF_ratio_ab_priorMemory_summary;
    temp_2 = temp_VAF_ratio_ac_priorMeta_summary;
    temp_3 = temp_VAF_ratio_bc_memoryMeta_summary;
    
    %     temp_1B = temp_VAF_ratio_aa_prior_summary_mean;
    %     temp_2B = temp_VAF_ratio_bb_memory_summary_mean;
    %     temp_3B = temp_VAF_ratio_cc_meta_summary_mean;
    temp_1B = temp_VAF_ratio_aa_prior_summary;
    temp_2B = temp_VAF_ratio_bb_memory_summary;
    temp_3B = temp_VAF_ratio_cc_meta_summary;
    
    
    %     temp1_SEM = std(temp_1);
    %     temp2_SEM = std(temp_2);
    %     temp3_SEM = std(temp_3);
    %
    %     temp1B_SEM = std(temp_1B);
    %     temp2B_SEM = std(temp_2B);
    %     temp3B_SEM = std(temp_3B);
    
    temp1_SEM = std(temp_1)/sqrt(length(temp_1));
    temp2_SEM = std(temp_2)/sqrt(length(temp_2));
    temp3_SEM = std(temp_3)/sqrt(length(temp_3));
    
    temp1B_SEM = std(temp_1B)/sqrt(length(temp_1B));
    temp2B_SEM = std(temp_2B)/sqrt(length(temp_2B));
    temp3B_SEM = std(temp_3B)/sqrt(length(temp_3B));
    
    
    
    temp_y_min = min([mean(temp_1)-temp1_SEM,mean(temp_2)-temp2_SEM,mean(temp_3)-temp3_SEM,...
        mean(temp_1B)-temp1B_SEM,mean(temp_2B)-temp2B_SEM,mean(temp_3B)-temp3B_SEM]);
    temp_y_max = max([mean(temp_1)+temp1_SEM,mean(temp_2)+temp2_SEM,mean(temp_3)+temp3_SEM,...
        mean(temp_1B)+temp1B_SEM,mean(temp_2B)+temp2B_SEM,mean(temp_3B)+temp3B_SEM]);
    
    plot([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],'color',[0 0 0],'linewidth',0.5);
    hold on
    errorbar([1 2 3], [mean(temp_1) mean(temp_2) mean(temp_3)],[temp1_SEM temp2_SEM temp3_SEM],...
        '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);
    hold on
    
    plot([4 5 6], [mean(temp_1B) mean(temp_2B) mean(temp_3B)],'color',[0 0 0],'linewidth',0.5);
    hold on
    errorbar([4 5 6], [mean(temp_1B) mean(temp_2B) mean(temp_3B)],[temp1B_SEM temp2B_SEM temp3B_SEM],...
        '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 10);
    hold on
    
    
    
    
    
    set(gca,'linewidth',1.5)
    
    %xlim([0.3 6.7]);
    xlim([0.5 6.5]);
    ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
    set(gca, 'FontSize', 8) %14
    set(gca,'box','off');% 取消右、上边框
    
    xticks([1 2 3 4 5 6]);
    
    xtl = ["BSL-WM", "BSL-Meta","WM-Meta","BSL-BSL", "WM-WM","Meta-Meta"];
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    xtext_xp=xt;
    %xtext_yp=temp_y_min*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.14;
    xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.18;
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',6.5);%9-->12
    
    set(gca,'xticklabel','');
    
    
    ylabel('VAF ratio', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('VAF between subspaces'),'fontsize',9);
    temp_title.Interpreter = 'none';
end




%% Two monkeys. (VAF) Plot Neuronal linear regression, using neuron activity to estimate task variables linearly.
if if_plotTwoMonkey_subspace == 1
    if if_plot == 1 && if_plot_part6 == 1
        
        fig = figure('Name','Two monkeys. VAF ratio (eachFOV)','NumberTitle','off');
        set(gcf,'Position',[1375+0 450+0 240.5*0.8*0.95 195.3]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_1D = temp_VAF_ratio_ab_priorMemory_summary_Ding;
        temp_2D = temp_VAF_ratio_ac_priorMeta_summary_Ding;
        temp_3D = temp_VAF_ratio_bc_memoryMeta_summary_Ding;
        temp_4D = temp_VAF_ratio_aa_prior_summary_Ding;
        temp_5D = temp_VAF_ratio_bb_memory_summary_Ding;
        temp_6D = temp_VAF_ratio_cc_meta_summary_Ding;
        
        temp_1Z = temp_VAF_ratio_ab_priorMemory_summary;
        temp_2Z = temp_VAF_ratio_ac_priorMeta_summary;
        temp_3Z = temp_VAF_ratio_bc_memoryMeta_summary;
        temp_4Z = temp_VAF_ratio_aa_prior_summary;
        temp_5Z = temp_VAF_ratio_bb_memory_summary;
        temp_6Z = temp_VAF_ratio_cc_meta_summary;
        
        if true
            [~,temp_p_13D] = ttest(temp_1D,temp_3D);
            [~,temp_p_23D] = ttest(temp_2D,temp_3D);
            
            [~,temp_p_13Z] = ttest(temp_1Z,temp_3Z);
            [~,temp_p_23Z] = ttest(temp_2Z,temp_3Z);
            
            [~,temp_p_12with3D] = ttest2([temp_1D;temp_2D],temp_3D);
            [~,temp_p_12with3Z] = ttest2([temp_1Z;temp_2Z],temp_3Z);
        end
        
        temp1D_SEM = std(temp_1D)/sqrt(length(temp_1D));
        temp2D_SEM = std(temp_2D)/sqrt(length(temp_2D));
        temp3D_SEM = std(temp_3D)/sqrt(length(temp_3D));
        temp4D_SEM = std(temp_4D)/sqrt(length(temp_4D));
        temp5D_SEM = std(temp_5D)/sqrt(length(temp_5D));
        temp6D_SEM = std(temp_6D)/sqrt(length(temp_6D));
        
        temp1Z_SEM = std(temp_1Z)/sqrt(length(temp_1Z));
        temp2Z_SEM = std(temp_2Z)/sqrt(length(temp_2Z));
        temp3Z_SEM = std(temp_3Z)/sqrt(length(temp_3Z));
        temp4Z_SEM = std(temp_4Z)/sqrt(length(temp_4Z));
        temp5Z_SEM = std(temp_5Z)/sqrt(length(temp_5Z));
        temp6Z_SEM = std(temp_6Z)/sqrt(length(temp_6Z));
        
        temp_y_min123 = min([mean(temp_1D)-temp1D_SEM,mean(temp_2D)-temp2D_SEM,mean(temp_3D)-temp3D_SEM,...
            mean(temp_1Z)-temp1Z_SEM,mean(temp_2Z)-temp2Z_SEM,mean(temp_3Z)-temp3Z_SEM]);
        temp_y_max123 = max([mean(temp_1D)+temp1D_SEM,mean(temp_2D)+temp2D_SEM,mean(temp_3D)+temp3D_SEM,...
            mean(temp_1Z)+temp1Z_SEM,mean(temp_2Z)+temp2Z_SEM,mean(temp_3Z)+temp3Z_SEM]);
        
        temp_y_min456 = min([mean(temp_4D)-temp4D_SEM,mean(temp_5D)-temp5D_SEM,mean(temp_6D)-temp6D_SEM,...
            mean(temp_4Z)-temp4Z_SEM,mean(temp_5Z)-temp5Z_SEM,mean(temp_6Z)-temp6Z_SEM]);
        temp_y_max456 = max([mean(temp_4D)+temp4D_SEM,mean(temp_5D)+temp5D_SEM,mean(temp_6D)+temp6D_SEM,...
            mean(temp_4Z)+temp4Z_SEM,mean(temp_5Z)+temp5Z_SEM,mean(temp_6Z)+temp6Z_SEM]);
        
        temp_y_min = min([temp_y_min123,temp_y_min456]);
        temp_y_max = max([temp_y_max123,temp_y_max456]);
       
        temp_lineWidth_errorbar = 1;%2
        temp_capSize = 7;%10
        temp_color_D = [1 1 1]*0.5;
        temp_color_Z = [1 1 1]*0;
        
        h_line = [];
        plot([1 2 3], [mean(temp_1D) mean(temp_2D) mean(temp_3D)],'color',temp_color_D,'linewidth',0.5);
        hold on
        h_line = [h_line errorbar([1 2 3], [mean(temp_1D) mean(temp_2D) mean(temp_3D)],[temp1D_SEM temp2D_SEM temp3D_SEM],...
            '.', 'Color', temp_color_D, 'MarkerFaceColor', temp_color_D, 'LineWidth', temp_lineWidth_errorbar, 'CapSize', temp_capSize)];
        hold on
        
        plot([4 5 6], [mean(temp_4D) mean(temp_5D) mean(temp_6D)],'color',temp_color_D,'linewidth',0.5);
        hold on
        errorbar([4 5 6], [mean(temp_4D) mean(temp_5D) mean(temp_6D)],[temp4D_SEM temp5D_SEM temp6D_SEM],...
            '.', 'Color', temp_color_D, 'MarkerFaceColor', temp_color_D, 'LineWidth', temp_lineWidth_errorbar, 'CapSize', temp_capSize);
        hold on
        
        
        plot([1 2 3], [mean(temp_1Z) mean(temp_2Z) mean(temp_3Z)],'color',temp_color_Z,'linewidth',0.5);
        hold on
        h_line = [h_line errorbar([1 2 3], [mean(temp_1Z) mean(temp_2Z) mean(temp_3Z)],[temp1Z_SEM temp2Z_SEM temp3Z_SEM],...
            '.', 'Color', temp_color_Z, 'MarkerFaceColor', temp_color_Z, 'LineWidth', temp_lineWidth_errorbar, 'CapSize', temp_capSize)];
        hold on
        
        plot([4 5 6], [mean(temp_4Z) mean(temp_5Z) mean(temp_6Z)],'color',temp_color_Z,'linewidth',0.5);
        hold on
        errorbar([4 5 6], [mean(temp_4Z) mean(temp_5Z) mean(temp_6Z)],[temp4Z_SEM temp5Z_SEM temp6Z_SEM],...
            '.', 'Color', temp_color_Z, 'MarkerFaceColor', temp_color_Z, 'LineWidth', temp_lineWidth_errorbar, 'CapSize', temp_capSize);
        hold on
        
        %le = legend(h_line,'Monkey D','Monkey Z','Location','northeast','fontsize',6.5);%10
        %le.ItemTokenSize = ones(1,2)*8;
        %legend('boxoff');
        
        set(gca,'linewidth',1.5)
        
        %xlim([0.3 6.7]);
        xlim([0.5 6.5]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);%0.3
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        xticks([1 2 3 4 5 6]);
        
        xtl = ["BSL-WM", "BSL-Meta","WM-Meta","BSL-BSL", "WM-WM","Meta-Meta"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=temp_y_min*ones(1,length(xt))+temp_y_min-(temp_y_max-temp_y_min)*0.14;
        xtext_yp=temp_y_min*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.23;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',6.5);%9-->12
        
        set(gca,'xticklabel','');
        
        
        ylabel('VAF ratio', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('VAF between subspaces'),'fontsize',9);
        temp_title = title(sprintf('Between subspaces'),'fontsize',9);
        
    end
end


%% Explained variance of neuron activity in three axes
if true
    EV_neuronPrior_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_EV_neuronPrior_summary;
    EV_neuronMemory_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_EV_neuronMemory_summary;
    EV_neuronMeta_summary = allPart_summary_mulFOV.part6_summary_mulFOV.temp_EV_neuronMeta_summary;
    
    EV_neuronThree_summary = [EV_neuronPrior_summary,EV_neuronMemory_summary,EV_neuronMeta_summary];
    
    
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name','EV','NumberTitle','off');
        set(gcf,'Position',[1100 45 240*0.80*0.93 336*1.11*0.67*0.9*1.2*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_1 = 100*EV_neuronPrior_summary';
        temp_2 = 100*EV_neuronMemory_summary';
        temp_3 = 100*EV_neuronMeta_summary';
        
        
        temp_y_min = min([temp_1 temp_2 temp_3]);
        temp_y_max = max([temp_1 temp_2 temp_3]);
        
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
        
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.3 3.7]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = ["BSL", "WM strength", "Meta"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.15;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
        
        set(gca,'xticklabel','');
        
        %yticks([0 40 80]);
        
        ylabel('EV (%)', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('Explained variance\nof neuron activity'),'fontsize',8);
        
    end
    
end


%% Explained variance of decoded variables (scores) in three axes
if true
    r2_prior_linearAxis_summary = allPart_summary_mulFOV.part6_summary_mulFOV.r2_prior_linearAxis_summary;
    r2_memory_linearAxis_summary = allPart_summary_mulFOV.part6_summary_mulFOV.r2_memory_linearAxis_summary;
    r2_meta_linearAxis_summary = allPart_summary_mulFOV.part6_summary_mulFOV.r2_meta_linearAxis_summary;
    
    r2_three_linearAxis_summary = [r2_prior_linearAxis_summary,r2_memory_linearAxis_summary,r2_meta_linearAxis_summary];
    
    
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name','EV_linearAxes','NumberTitle','off');
        %set(gcf,'Position',[1100 45 240*0.80*0.93 336*1.11*0.67*0.9*1.2*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[1100 45 240*0.80*0.93 336*1.11*0.67*0.9*1.2*1.2*0.7]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
        nexttile
        
        temp_1 = 100*r2_prior_linearAxis_summary';
        temp_2 = 100*r2_memory_linearAxis_summary';
        temp_3 = 100*r2_meta_linearAxis_summary';
        
        
        temp_y_min = min([temp_1 temp_2 temp_3]);
        temp_y_max = max([temp_1 temp_2 temp_3]);
        
        temp_y_min = 0;
        temp_y_max = 100;
        
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
        
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.3 3.7]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = ["Baseline", "WM strength", "Meta"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.15;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
        
        set(gca,'xticklabel','');
        
        %yticks([0 40 80]);
        
        ylabel('EV (%)', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('Explained variance'),'fontsize',9);
        
    end
    
end



%% PCA results of neuron activity
if true
    pca_explained_baseline_cell_summary = allPart_summary_mulFOV.part6_summary_mulFOV.pca_explained_baseline_cell_summary;
    pca_explained_delay_cell_summary = allPart_summary_mulFOV.part6_summary_mulFOV.pca_explained_delay_cell_summary;
    pca_r_eachPC_VS_eventVector_baseline_cell_summary = allPart_summary_mulFOV.part6_summary_mulFOV.pca_r_eachPC_VS_eventVector_baseline_cell_summary;
    pca_r_eachPC_VS_eventVector_delay_cell_summary = allPart_summary_mulFOV.part6_summary_mulFOV.pca_r_eachPC_VS_eventVector_delay_cell_summary;
    
    temp_FOVNum = length(pca_explained_baseline_cell_summary);
    
    totalEVofRelatedPC_eachFeature_baseline = nan(temp_FOVNum,4);
    totalEVofRelatedPC_eachFeature_delay = totalEVofRelatedPC_eachFeature_baseline;
    for tempi=1:temp_FOVNum
        temp_pca_explained_baseline = pca_explained_baseline_cell_summary{tempi};
        temp_pca_r_eachPC_VS_eventVector_baseline = pca_r_eachPC_VS_eventVector_baseline_cell_summary{tempi};
        
        temp_pca_explained_delay = pca_explained_delay_cell_summary{tempi};
        temp_pca_r_eachPC_VS_eventVector_delay = pca_r_eachPC_VS_eventVector_delay_cell_summary{tempi};
        
        
        tempBoolIndex_A = abs(temp_pca_r_eachPC_VS_eventVector_baseline) > 0;
        tempBoolIndex_B = abs(temp_pca_r_eachPC_VS_eventVector_delay) > 0;
        
        for tempj=1:4
            totalEVofRelatedPC_eachFeature_baseline(tempi,tempj) = sum(tempBoolIndex_A(:,tempj).*temp_pca_explained_baseline);
            totalEVofRelatedPC_eachFeature_delay(tempi,tempj) = sum(tempBoolIndex_B(:,tempj).*temp_pca_explained_delay);
        end
        
    end
    totalEVofRelatedPC_eachFeature_baseline;
    totalEVofRelatedPC_eachFeature_delay;
    
    
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name','totalEVofRelatedPC','NumberTitle','off');
        set(gcf,'Position',[1100 45 240*0.80*0.93*1.8 336*1.11*0.67*0.9*1.2*1.2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        t.Title.String = sprintf('Total EV of related PC (p < 0.001)');
        t.Title.FontSize = 9;
        
        
        %totalEVofRelatedPC_eachFeature_baseline
        nexttile
        
        temp_1 = totalEVofRelatedPC_eachFeature_baseline(:,1)';
        temp_2 = totalEVofRelatedPC_eachFeature_baseline(:,2)';
        temp_3 = totalEVofRelatedPC_eachFeature_baseline(:,3)';
        temp_4 = totalEVofRelatedPC_eachFeature_baseline(:,4)';
        
        
        temp_y_min = min([temp_1 temp_2 temp_3 temp_4]);
        temp_y_max = max([temp_1 temp_2 temp_3 temp_4]);
        
        temp_data = [temp_1';temp_2';temp_3';temp_4'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        g4 = repmat({'D'},length(temp_4),1);
        
        temp_label = [g1;g2;g3;g4];
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 4, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        h(3).ViolinPlot.FaceAlpha = 0.1;
        h(4).ViolinPlot.FaceAlpha = 0.1;
        
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.3 4.7]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = ["BSL", "WM strength", "Meta", "Time"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.15;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
        
        set(gca,'xticklabel','');
        
        
        ylabel('Total EV (%)', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('Baseline'),'fontsize',7.5);
        
        
        %totalEVofRelatedPC_eachFeature_delay
        nexttile
        
        temp_1 = totalEVofRelatedPC_eachFeature_delay(:,1)';
        temp_2 = totalEVofRelatedPC_eachFeature_delay(:,2)';
        temp_3 = totalEVofRelatedPC_eachFeature_delay(:,3)';
        temp_4 = totalEVofRelatedPC_eachFeature_delay(:,4)';
        
        
        temp_y_min = min([temp_1 temp_2 temp_3 temp_4]);
        temp_y_max = max([temp_1 temp_2 temp_3 temp_4]);
        
        temp_data = [temp_1';temp_2';temp_3';temp_4'];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        g3 = repmat({'C'},length(temp_3),1);
        g4 = repmat({'D'},length(temp_4),1);
        
        temp_label = [g1;g2;g3;g4];
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 4, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'};{'C'};{'D'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        h(3).ViolinPlot.FaceAlpha = 0.1;
        h(4).ViolinPlot.FaceAlpha = 0.1;
        
        
        
        set(gca,'linewidth',1.5)
        
        xlim([0.3 4.7]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.1 temp_y_max+(temp_y_max-temp_y_min)*0.1]);
        set(gca, 'FontSize', 8) %14
        set(gca,'box','off');% 取消右、上边框
        
        
        xtl = ["BSL", "WM strength", "Meta", "Time"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(temp_y_min)*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.15;
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%9-->12
        
        set(gca,'xticklabel','');
        
        
        %ylabel('Total EV (%)', 'FontSize', 9, 'FontWeight', 'normal');
        temp_title = title(sprintf('Delay'),'fontsize',7.5);
        
    end
    
    
    
    totalEVofRelatedPC_eachFeature_baseline;
    totalEVofRelatedPC_eachFeature_delay;
    
    temp_exampleFOVID = 4;
    temp_PCLimit = 25;
    temp_pca_explained_baseline_example = pca_explained_baseline_cell_summary{temp_exampleFOVID};
    temp_pca_r_eachPC_VS_eventVector_baseline_example = pca_r_eachPC_VS_eventVector_baseline_cell_summary{temp_exampleFOVID};
    temp_pca_explained_delay_example = pca_explained_delay_cell_summary{temp_exampleFOVID};
    temp_pca_r_eachPC_VS_eventVector_delay_example = pca_r_eachPC_VS_eventVector_delay_cell_summary{temp_exampleFOVID};
    
    %% Correlations between features and PCs of example FOV (baseline)
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name','Correlations between features and PCs of example FOV (baseline)','NumberTitle','off');
        %set(gcf,'Position',[10 50 390 103*1.9]);
        set(gcf,'Position',[10 50 390 103*1.9*0.85]);
        
        t = tiledlayout(4,31,'TileSpacing','tight','Padding','compact');
        
        nexttile([2 3])
        set(gca,'Visible','off');
        
        nexttile([2 27])
        
        temp_EV = temp_pca_explained_baseline_example(1:temp_PCLimit);
        
        C = temp_pca_r_eachPC_VS_eventVector_baseline_example(1:temp_PCLimit,:)';
        %C_max = 1;
        C = abs(C);
        %C_max = max(C,[],'all');
        %C_min = min(C,[],'all');
        C_max = 1;
        C_min = 0;
        
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
        hold on
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 8) %12
        
        
        %set(gca,'XTick',1:size(C,2));
        set(gca,'xticklabel','');
        
        set(gca,'YTick',1:size(C,1));
        ytl = ["BSL","WM strength","Meta","Time"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
        ytext_xp=(xt(1))*ones(1,length(yt))-4.8;%-1.6
        ytext_yp=yt;
        text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',7.5);%8,7
        set(gca,'yticklabel','');
        
        set(gca,'TickLength',[0 0]);
        
        x_lim = xlim;
        y_lim = ylim;
        xlabel(sprintf('PC (1~25, EV = %.0f%%)',sum(temp_EV)),'Position',[mean(x_lim) max(y_lim)+0.2],'FontSize', 9);
        temp_title = title(sprintf('Correlations with PCs of one FOV (baseline)'),'FontSize',9);%11
        
        c = colorbar('westoutside','FontSize',6.5);
        c.Position = c.Position+[0.832 0.056 -0.02 -0.065];
        c.Ticks = [0 1];
        
        colormap(coolwarm());
        
        
        
        nexttile([2 3])
        set(gca,'Visible','off');
        
        nexttile([2 27])
        
        temp_1 = temp_EV;
        
        temp_x = 1:length(temp_1);
        
        [temp_x_min,temp_x_max] = bounds(temp_x);
        [temp_y_min,temp_y_max] = bounds(temp_1);
        
        temp_width = 0.8;
        temp_x_min = temp_x_min - temp_width/2;
        temp_x_max = temp_x_max + temp_width/2;
        
        bar(temp_x, temp_1,temp_width, ...
            'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        
        
        xlim([temp_x_min temp_x_max]);
        ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.25]);
        
        set(gca,'xticklabel','');
        h=gca;
        h.XAxis.TickLength = [0 0];
        h.YAxis.TickLength = [0 0];
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 8) %14
        ylabel('EV (%)', 'FontSize', 9, 'FontWeight', 'normal');
        
    end
    
    
    %% Correlations between features and PCs of example FOV (delay)
    if if_plot == 1 && if_plot_part6 == 1
        fig = figure('Name','Correlations between features and PCs of example FOV (delay)','NumberTitle','off');
        set(gcf,'Position',[450 50 390 103*1.9*0.85]);
        
        t = tiledlayout(4,31,'TileSpacing','tight','Padding','compact');
        
        nexttile([2 3])
        set(gca,'Visible','off');
        
        nexttile([2 27])
        
        temp_EV = temp_pca_explained_delay_example(1:temp_PCLimit);
        
        C = temp_pca_r_eachPC_VS_eventVector_delay_example(1:temp_PCLimit,:)';
        %C_max = 1;
        C = abs(C);
        %C_max = max(C,[],'all');
        %C_min = min(C,[],'all');
        
        C_max = 1;
        C_min = 0;
        
        imagesc([1 size(C,2)],[1 size(C,1)],C,[C_min C_max]);
        hold on
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 8) %12
        
        
        %set(gca,'XTick',1:size(C,2));
        set(gca,'xticklabel','');
        
        set(gca,'YTick',1:size(C,1));
        ytl = ["BSL","WM strength","Meta","Time"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
        ytext_xp=(xt(1))*ones(1,length(yt))-4.8;%-1.6
        ytext_yp=yt;
        text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',7.5);%8,7
        set(gca,'yticklabel','');
        
        set(gca,'TickLength',[0 0]);
        
        x_lim = xlim;
        y_lim = ylim;
        xlabel(sprintf('PC (1~25, EV = %.0f%%)',sum(temp_EV)),'Position',[mean(x_lim) max(y_lim)+0.2],'FontSize', 9);
        temp_title = title(sprintf('Correlations with PCs of one FOV (delay)'),'FontSize',9);%11
        
        c = colorbar('westoutside','FontSize',6.5);
        c.Position = c.Position+[0.832 0.056 -0.02 -0.065];
        c.Ticks = [0 1];
        
        colormap(coolwarm());
        
        
        
        nexttile([2 3])
        set(gca,'Visible','off');
        
        nexttile([2 27])
        
        temp_1 = temp_EV;
        
        temp_x = 1:length(temp_1);
        
        [temp_x_min,temp_x_max] = bounds(temp_x);
        [temp_y_min,temp_y_max] = bounds(temp_1);
        
        temp_width = 0.8;
        temp_x_min = temp_x_min - temp_width/2;
        temp_x_max = temp_x_max + temp_width/2;
        
        bar(temp_x, temp_1,temp_width, ...
            'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        
        
        xlim([temp_x_min temp_x_max]);
        ylim([max(0,temp_y_min-(temp_y_max-temp_y_min)*0.25) temp_y_max+(temp_y_max-temp_y_min)*0.25]);
        
        set(gca,'xticklabel','');
        h=gca;
        h.XAxis.TickLength = [0 0];
        h.YAxis.TickLength = [0 0];
        
        set(gca,'linewidth',1.5)
        set(gca,'box','off');% 取消右、上边框
        set(gca, 'FontSize', 8) %14
        ylabel('EV (%)', 'FontSize', 9, 'FontWeight', 'normal');
        
    end
    
    
end


%% End