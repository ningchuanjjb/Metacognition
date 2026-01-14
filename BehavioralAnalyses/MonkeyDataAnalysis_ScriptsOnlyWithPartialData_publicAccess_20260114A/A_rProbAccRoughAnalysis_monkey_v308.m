%This script: To conduct behavioral analysis
%% Initialization
% clear all %#ok<CLALL>
close all

%% Some flags to control the following script

if_plot_choiceMinusNoChoice_inOne = 0;
% if_monkey_D0_Z1 = 1;% To decide whether dealing with Ding's data or Zelku's data
if_plot_partialAccuracy = 1;
if_lengthError_down0_up1_both2 = 2;
if_responseStructure_noChoice0_choice1 = 0;
if_ruleOut_deleteAddOne0_others1 = 1;
% if_numFrames8_ruleOut67 = 1;
if_newBoost = 0;
if_fig8_filter = 1;
if_fig8_boost_negative0_positive1_full2 = 2;
if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 = 2;% only work in negative boost.
if_additionFilter_lowGAcc = 1;
if_gAcc_errorBar = 1;

if_computePupil = 0;

if if_monkey_D0_Z1 == 0
    %     offloading_threshold = 0.977;% Ding, numFrames6
    if if_additionFilter_lowGAcc == 0
        offloading_threshold = 0.963;% Ding, numFrames8
    else
        %offloading_threshold = 0.914;% Ding, numFrames8
        %offloading_threshold = 0.95;% Ding, numFrames7
        offloading_threshold = 0.92;% Ding, numFrames6 2pB
        %offloading_threshold = 0.98;
    end
    if_numFrames8_ruleOut67 = 1;
elseif if_monkey_D0_Z1 == 1
    if if_additionFilter_lowGAcc == 0
        offloading_threshold = 0.86;% for Zelku, 0.85-->0.86
    elseif if_additionFilter_lowGAcc == 1
        %offloading_threshold = 0.987;% for Zelku, numFrames6A
        %offloading_threshold = 0.92;% for Zelku, numFrames6C
        %offloading_threshold = 0.98;
        offloading_threshold = 0.999;
    end
    if_numFrames8_ruleOut67 = 0;
end


ifFigure1 = 0;
ifFigure2 = 1;
ifFigure3 = 1;
ifFigure4 = 1;
ifFigure5 = 1;
ifFigure6 = 1;
ifFigure7 = 1;
ifFigure8 = 0;
ifFigure9 = 1;
ifFigure10 = 0;
ifFigure11 = 1;
ifFigure12 = 1;
ifFigure13 = 1;
ifFigure14 = 0;
ifFigure15 = 0;
ifFigure16 = 0;
ifFigure17 = 0;
ifFigure18 = 1;
ifFigure19 = 1;
ifFigure20 = 1;
ifFigure21 = 1;
ifFigure22 = 0;
ifFigure23 = 0;
ifFigure24 = 0;

ifSave_paperFig2Data = 0;

ifSave_fig = zeros(1, 30);
% ifSave_fig = ones(1, 3);

% ifSave_fig(1) = 1;
% ifSave_fig(2) = 1;
% ifSave_fig(3) = 1;
% ifSave_fig(8) = 1;

path_preprocessed = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\preprocessed';
path_results = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data';

% path_preprocessed = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\MonkeyDataAnalysis_ScriptsOnly_20260114A\Data\preprocessed';
% path_results = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\MonkeyDataAnalysis_ScriptsOnly_20260114A\Data';


if if_monkey_D0_Z1 == 0
    monkey_name = 'D';
    
    %     searchName = 'preprocessed_DnumFrames6_A_of_';
    %     searchName_gAcc = 'from21-12-06to21-12-13_D_gAcc_1';
    %     searchName_rProb = 'from21-12-06to21-12-13_D_offloadingProb_1';
    %     searchName_mappingParam = 'from21-12-06to21-12-13_D_mappingParam_1.mat';
    %     searchName_responseMatrix = 'D_responseMatrix_1.mat';
    
    
    
    %         searchName = 'preprocessed_DnumFrames7_A_of_';
    %         searchName_gAcc = 'from22-11-03to22-11-10_D_gAcc_1';
    %         searchName_rProb = 'from22-11-03to22-11-10_D_offloadingProb_1';
    %         searchName_mappingParam = 'from22-11-03to22-11-10_D_mappingParam_1';
    %         searchName_responseMatrix = 'D_responseMatrix_numFrames7_1';
    %         searchName_endingHold = 'endingHold_D_numFrames7';
    
    % ruleOutLength3 = [2,6,7,8,9,13,20,23,29,32,33];
    % isRuleOutLength3 = true(1,35);
    % isRuleOutLength3(ruleOutLength3) = false;
    
%     searchName = 'preprocessed_DnumFrames8_A_of_';
%     searchName_gAcc = 'from22-06-21to22-08-12_D_gAcc_1';
%     searchName_rProb = 'from22-06-21to22-08-12_D_offloadingProb_1';
%     searchName_mappingParam = 'from22-06-21to22-08-12_D_mappingParam_1';
%     searchName_responseMatrix = 'D_responseMatrix_2.mat';
%     searchName_endingHold = 'endingHold_D_A';
    
    searchName_pupil = 'D_pupil_all_filter_minusBaseline_mm';
    searchName_positive_negative_pupil = 'D_positive_negative_pupilSIze';
    searchName_positive_negative_seq_pupil = 'D_positive_negative_seq_pupilSIze';
    searchName_positive_negative_choiceMemory_seq_pupil = 'D_positive_negative_choiceMemory_seq_pupilSIze';
    searchName_positive_negative_choiceMemory_seq_pupil_timeCourse = ...
        'D_positive_negative_choiceMemory_seq_timecourse';
    
    
    %searchName = 'preprocessed_DnumFrames6_Train1';
    
    %     searchName = 'preprocessed_DnumFrames6_Train2';
    %     searchName_gAcc = 'from23-03-13to23-04-05_D_gAcc_1';
    %     searchName_rProb = 'from23-03-13to23-04-05_D_offloadingProb_1';
    %     searchName_mappingParam = 'from23-03-13to23-04-05_D_mappingParam_1';
    %     searchName_responseMatrix = 'D_responseMatrix_Train2';
    %     searchName_endingHold = 'endingHold_D_Train2';
    
    %     searchName = 'preprocessed_DnumFrames6_2p_of';
    %     searchName_gAcc = 'from23-04-26to23-06-12_D_gAcc_1';
    %     searchName_rProb = 'from23-04-26to23-06-12_D_offloadingProb_1';
    %     searchName_mappingParam = 'from23-04-26to23-06-12_D_mappingParam_1';
    %     searchName_responseMatrix = 'D_responseMatrix_2p_A.mat';
    %     searchName_endingHold = 'endingHold_D_2p_A';
    
    searchName = 'preprocessed_DnumFrames6_2p_of';
    searchName_gAcc = 'from23-04-26to23-06-20_D_gAcc_1';
    searchName_rProb = 'from23-04-26to23-06-20_D_offloadingProb_1';
    searchName_mappingParam = 'from23-04-26to23-06-20_D_mappingParam_1';
    searchName_responseMatrix = 'D_responseMatrix_2p_B.mat';
    searchName_endingHold = 'endingHold_D_2p_B';
    
    
elseif if_monkey_D0_Z1 == 1
    monkey_name = 'Z';
    
    %     searchName = 'preprocessed_ZnumFrames6_B_of_';
    %     searchName_gAcc = 'from22-06-20to22-07-01_Z_gAcc_1';
    %     searchName_rProb = 'from22-06-20to22-07-01_Z_offloadingProb_1';
    %     searchName_mappingParam = 'from22-06-20to22-07-01_Z_mappingParam_1';
    %     searchName_responseMatrix = 'Z_responseMatrix_B2';
    %     searchName_endingHold = 'endingHold_Z_B2';
    %     searchName_pupil = 'Z_pupil_all_filter_minusBaseline_mm';
    %     searchName_positive_negative_pupil = 'Z_positive_negative_pupilSIze';
    %     searchName_positive_negative_seq_pupil = 'Z_positive_negative_seq_pupilSIze';
    %     searchName_positive_negative_choiceMemory_seq_pupil = 'Z_positive_negative_choiceMemory_seq_pupilSIze';
    %     searchName_positive_negative_choiceMemory_seq_pupil_timeCourse = ...
    %         'Z_positive_negative_choiceMemory_seq_timecourse';
    
%     searchName = 'preprocessed_ZnumFrames6_A_of_';
%     searchName_gAcc = 'from21-06-16to21-08-04_Z_gAcc_1';
%     searchName_rProb = 'from21-06-16to21-08-04_Z_offloadingProb_1';
%     searchName_mappingParam = 'from21-06-16to21-08-04_Z_mappingParam_1';
%     searchName_responseMatrix = 'Z_responseMatrix_2';
%     searchName_endingHold = 'endingHold_Z_A';
    searchName_pupil = 'Z_pupil_all_filter_minusBaseline_mm';
    searchName_positive_negative_pupil = 'Z_positive_negative_pupilSIze';
    searchName_positive_negative_seq_pupil = 'Z_positive_negative_seq_pupilSIze';
    searchName_positive_negative_choiceMemory_seq_pupil = 'Z_positive_negative_choiceMemory_seq_pupilSIze';
    searchName_positive_negative_choiceMemory_seq_pupil_timeCourse = ...
        'Z_positive_negative_choiceMemory_seq_timecourse';
    
    %     searchName = 'preprocessed_ZnumFrames6_Train1';
    %     searchName_gAcc = 'from23-03-26to23-04-03_Z_gAcc_1';
    %     searchName_rProb = 'from23-03-26to23-04-03_Z_offloadingProb_1';
    %     searchName_mappingParam = 'from23-03-26to23-04-03_Z_mappingParam_1';
    %     searchName_responseMatrix = 'Z_responseMatrix_Train1';
    %     searchName_endingHold = 'endingHold_Z_Train1';
    
    %     searchName = 'preprocessed_ZnumFrames6_Train2';
    %     searchName_gAcc = 'from23-04-05to23-04-10_Z_gAcc_1';
    %     searchName_rProb = 'from23-04-05to23-04-10_Z_offloadingProb_1';
    %     searchName_mappingParam = 'from23-04-05to23-04-10_Z_mappingParam_1';
    %     searchName_responseMatrix = 'Z_responseMatrix_Train2';
    %     searchName_endingHold = 'endingHold_Z_Train2';
    
    %     searchName = 'preprocessed_ZnumFrames6_C';
    %     searchName_gAcc = 'from23-03-29to23-05-19_Z_gAcc_1';
    %     searchName_rProb = 'from23-03-29to23-05-19_Z_offloadingProb_1';
    %     searchName_mappingParam = 'from23-03-29to23-05-19_Z_mappingParam_1';
    %     searchName_responseMatrix = 'Z_responseMatrix_6C';
    %     searchName_endingHold = 'endingHold_Z_6C';
    
    %     searchName = 'preprocessed_ZnumFrames6_D';
    %     searchName_gAcc = 'from23-08-02to23-08-18_Z_gAcc_1';
    %     searchName_rProb = 'from23-08-02to23-08-18_Z_offloadingProb_1';
    %     searchName_mappingParam = 'from23-08-02to23-08-18_Z_mappingParam_1';
    %     searchName_responseMatrix = 'Z_responseMatrix_6D';
    %     searchName_endingHold = 'endingHold_Z_6D';
    
    
%     searchName = 'preprocessed_ZnumFrames6_early2p';
%     %searchName_gAcc = 'from01-08to01-23_Z_gAcc_1';
%     %searchName_rProb = 'from01-08to01-23_Z_offloadingProb_1';
%     searchName_gAcc = 'from01-08to01-31_Z_gAcc_1';
%     searchName_rProb = 'from01-08to01-31_Z_offloadingProb_1';
%     searchName_responseMatrix = 'Z_responseMatrix_early2p.mat';
%     searchName_mappingParam = 'from23-08-02to23-08-18_Z_mappingParam_1';
%     searchName_endingHold = 'endingHold_Z_6D';
    
    searchName = 'preprocessed_ZnumFrames6_all2p';
    searchName_gAcc = 'from01-08to04-11_Z_gAcc_1';
    searchName_rProb = 'from01-08to04-11_Z_offloadingProb_1';
    searchName_responseMatrix = 'Z_responseMatrix_all2p.mat';
    searchName_mappingParam = 'from23-08-02to23-08-18_Z_mappingParam_1';
    searchName_endingHold = 'endingHold_Z_6D';

%     searchName = 'preprocessed_ZnumFrames6_almost2p';
%     searchName_gAcc = 'from01-08to03-30_Z_gAcc_1';
%     searchName_rProb = 'from01-08to03-30_Z_offloadingProb_1';
%     searchName_responseMatrix = 'Z_responseMatrix_almost2p.mat';
%     searchName_mappingParam = 'from23-08-02to23-08-18_Z_mappingParam_1';
%     searchName_endingHold = 'endingHold_Z_6D';    
    
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


% Load other processed results
load_gAcc = loadMat_single(searchName_gAcc, path_results);

gAcc_noChoice = load_gAcc.gAcc_noChoice;
gAcc_choice = load_gAcc.gAcc_choice;
gAcc_noChoice_merged = load_gAcc.gAcc_noChoice_merged;
gAcc_choice_merged = load_gAcc.gAcc_choice_merged;
gAcc_noChoice_collapsed = load_gAcc.gAcc_noChoice_collapsed;
gAcc_choice_collapsed = load_gAcc.gAcc_choice_collapsed;
gAcc_noChoice_collapsed_inOne = load_gAcc.gAcc_noChoice_collapsed_inOne;
gCorrect_choice_trial_count_collapsed = load_gAcc.gCorrect_choice_trial_count_collapsed;
g_choice_trial_count_collapsed = load_gAcc.g_choice_trial_count_collapsed;
gCorrect_noChoice_trial_count_collapsed = load_gAcc.gCorrect_noChoice_trial_count_collapsed;
g_noChoice_trial_count_collapsed = load_gAcc.g_noChoice_trial_count_collapsed;


pointKindsNum = size(gAcc_noChoice_collapsed, 2);

gAcc_choice_collapsed_inOne = [];
for tempi=1:pointKindsNum
    gAcc_choice_collapsed_inOne = [gAcc_choice_collapsed_inOne gAcc_choice_collapsed{tempi}'];
end

numSeq = zeros(1, pointKindsNum);
for target_seqLength=1:pointKindsNum
    numSeq(target_seqLength) = size(gAcc_noChoice_collapsed{1, target_seqLength}, 1);
end
numFrames = numSeq(1);
if if_numFrames8_ruleOut67 == 0
    numSeq_extend = zeros(1, numFrames-1);
    %elseif numFrames == 8 && if_numFrames8_ruleOut67 == 1
elseif if_numFrames8_ruleOut67 == 1
    numSeq_extend = zeros(1, 5);
end
numSeq_extend(1:pointKindsNum) = numSeq;
for tempj=1:length(numSeq_extend)
    if numSeq_extend(tempj) == 0
        numSeq_extend(tempj) = numSeq_extend(numFrames - tempj);
    end
end

% Load other processed results
load_rProb = loadMat_single(searchName_rProb, path_results);
offloadingProb_collapsed = load_rProb.offloadingProb_all;

offloadingProb_inOne = [];
for tempi=1:pointKindsNum
    offloadingProb_inOne = [offloadingProb_inOne offloadingProb_collapsed{tempi}'];
end

% offloadingProb_inOne = offloadingProb_all_inOne;
offloadingProb_all = load_rProb.offloadingProb_all;

% Load other processed results
load_mapping = loadMat_single(searchName_mappingParam, path_results);
target_seqSet = load_mapping.target_seqSet;
seqSet_inOne = load_mapping.seqSet_inOne;

target_seqSet_inOne = [];
seqSet_inOne_inOne = [];
for tempi=1:pointKindsNum
    target_seqSet_inOne = [target_seqSet_inOne; target_seqSet{tempi}];
    seqSet_inOne_inOne = [seqSet_inOne_inOne; seqSet_inOne{tempi}'];
end

if if_numFrames8_ruleOut67 == 0
    target_seqSet_extend = cell(1, numFrames-1);
    %elseif numFrames == 8 && if_numFrames8_ruleOut67 == 1
elseif if_numFrames8_ruleOut67 == 1
    target_seqSet_extend = cell(1, 5);
end
target_seqSet_extend(1:pointKindsNum) = target_seqSet;
for tempj=1:length(numSeq_extend)
    if isempty(target_seqSet_extend{tempj})
        for tempk=1:numSeq_extend(tempj)
            %temp1 = target_seqSet_extend{numFrames - tempj}{numFrames-tempk+1};
            temp1 = target_seqSet_extend{numFrames - tempj}{numSeq(numFrames - tempj)-tempk+1};
            
            temp2 = 1:numFrames;
            target_seqSet_extend{tempj}{tempk, 1} = temp2(~ismember(temp2, temp1));
        end
    end
end

target_seqSet_extend_inOne = [];
for tempi=1:length(numSeq_extend)
    target_seqSet_extend_inOne = [target_seqSet_extend_inOne; target_seqSet_extend{tempi}];
end

if if_numFrames8_ruleOut67 == 0
    seqSet_inOne_extend = cell(1, numFrames-1);
elseif numFrames == 8 && if_numFrames8_ruleOut67 == 1
    seqSet_inOne_extend = cell(1, 5);
end

for order_index=1:length(numSeq_extend)
    current_numSeq = numSeq_extend(order_index);
    seqSet_inOne_extend{order_index} = zeros(1, current_numSeq);
    for tempi=1:current_numSeq
        for tempj=1:order_index
            seqSet_inOne_extend{order_index}(tempi) = 10^(order_index-tempj) * ...
                target_seqSet_extend{order_index}{tempi}(tempj) + seqSet_inOne_extend{order_index}(tempi);
        end
    end
end

seqSet_inOne_extend_inOne = [];
for tempi=1:length(numSeq_extend)
    seqSet_inOne_extend_inOne = [seqSet_inOne_extend_inOne; seqSet_inOne_extend{tempi}'];
end

% Load other processed results
load_responseMatrix = loadMat_single(searchName_responseMatrix, path_results);

if if_numFrames8_ruleOut67 == 0
    responseMatrix_noChoice = load_responseMatrix.responseMatrix_noChoice;
    responseMatrix_choice = load_responseMatrix.responseMatrix_choice;
    %elseif numFrames == 8 && if_numFrames8_ruleOut67 == 1
elseif if_numFrames8_ruleOut67 == 1
    responseMatrix_noChoice = load_responseMatrix.responseMatrix_noChoice(:, 1:sum(numSeq_extend(1:5)));
    responseMatrix_choice = load_responseMatrix.responseMatrix_choice(:, 1:sum(numSeq_extend(1:5)));
end


% Load other processed results
load_endingHold = loadMat_single(searchName_endingHold, path_results);
endingHold = load_endingHold.endingHold;

% Load other processed results
load_pupil = loadMat_single(searchName_pupil, path_results);
load_pupil = load_pupil.(cell2mat(fieldnames(load_pupil)));
% pupilSize_seqSet = load_pupil.Delaystart_mean;
pupilSize_seqSet = cell(1, length(load_pupil.Delaystart_mean));
load_pupil.Sample_mean;
load_pupil.Delaystart_mean;

pupilSize_mean_seqSet = zeros(1, length(pupilSize_seqSet));
for tempi=1:length(pupilSize_seqSet)
    %     pupilSize_mean_seqSet(tempi) = mean([load_pupil.Sample_mean{tempi} load_pupil.Delaystart_mean{tempi}]);
    %     pupilSize_mean_seqSet(tempi) = mean(load_pupil.Delaystart_mean{tempi});
    %     pupilSize_mean_seqSet(tempi) = mean(load_pupil.Delayend_mean{tempi});
    %     pupilSize_mean_seqSet(tempi) = mean(load_pupil.Sample_mean{tempi});
    %     pupilSize_mean_seqSet(tempi) = mean(load_pupil.Sample_mean{tempi}(end+1-600:end));
    %     pupilSize_mean_seqSet(tempi) = mean(load_pupil.Sample_mean{tempi}(end+1-300:end));
    %     pupilSize_mean_seqSet(tempi) = mean([load_pupil.Sample_mean{tempi}(end+1-300:end) load_pupil.Delaystart_mean{tempi}(1:225)]);
    pupilSize_mean_seqSet(tempi) = mean([load_pupil.Sample_mean{tempi}(end+1-300:end) load_pupil.Delaystart_mean{tempi}(1:300)]);
    %     pupilSize_mean_seqSet(tempi) = mean([load_pupil.Sample_mean{tempi}(end+1-300:end) load_pupil.Delaystart_mean{tempi}(1:450)]);
    %     pupilSize_mean_seqSet(tempi) = mean([load_pupil.Delaystart_mean{tempi}(25:375)]);
end
% for tempi=1:length(pupilSize_seqSet)
%     pupilSize_mean_seqSet(tempi) = mean(pupilSize_seqSet{tempi});
%     %pupilSize_mean_seqSet(tempi) = mean(pupilSize_seqSet{tempi}(85+1:end-85));
%     %pupilSize_mean_seqSet(tempi) = mean(pupilSize_seqSet{tempi}(85+1:300));
% end

load_positive_negative_pupil = ...
    loadMat_single(searchName_positive_negative_pupil, path_results);
positive_negative_pupilSize = load_positive_negative_pupil.positive_negative_pupilSize;

load_positive_negative_seq_pupil = ...
    loadMat_single(searchName_positive_negative_seq_pupil, path_results);
positive_negative_seq_pupilSize = load_positive_negative_seq_pupil.positive_negative_seq;

load_positive_negative_choiceMemory_seq_pupil = ...
    loadMat_single(searchName_positive_negative_choiceMemory_seq_pupil, path_results);
positive_negative_choiceMemory_seq_pupilSize = load_positive_negative_choiceMemory_seq_pupil.positive_negative_seq;

load_positive_negative_pupilSize_seqTimeCourse = ...
    loadMat_single(searchName_positive_negative_choiceMemory_seq_pupil_timeCourse, path_results);
positive_negative_pupilSize_seqTimeCourse = load_positive_negative_pupilSize_seqTimeCourse.positive_negative_pupilSize_seqTimeCourse;

positive_negative_choiceMemory_seq_pupilSize_ROI = positive_negative_choiceMemory_seq_pupilSize;

fields = {'positive_correct_sem','positive_error_sem','negative_correct_sem','negative_error_sem'};
positive_negative_choiceMemory_seq_pupilSize_ROI = ...
    rmfield(positive_negative_choiceMemory_seq_pupilSize_ROI,fields);

temp_length_positive = length(positive_negative_choiceMemory_seq_pupilSize_ROI.positive_correct_mean);
temp_length_negative = length(positive_negative_choiceMemory_seq_pupilSize_ROI.negative_correct_mean);

% temp_ROI = -300-450;% Last item and delay1 450ms.
% temp_ROI = -450;% Delay1 450ms.

% temp_ROI_start = -300-450;
% temp_ROI_end = -450;

temp_ROI_start = -450;
temp_ROI_end = -150;

temp_ROI = temp_ROI_start - temp_ROI_end;

for tempi=1:temp_length_positive
    temp_correct = positive_negative_pupilSize_seqTimeCourse.positive_correct_mean{tempi};
    temp_error = positive_negative_pupilSize_seqTimeCourse.positive_error_mean{tempi};
    positive_negative_choiceMemory_seq_pupilSize_ROI.positive_correct_mean(tempi) = nan;
    positive_negative_choiceMemory_seq_pupilSize_ROI.positive_error_mean(tempi) = nan;
    %if ~isnan(temp_correct)
    if length(temp_correct) >= abs(temp_ROI)
        positive_negative_choiceMemory_seq_pupilSize_ROI.positive_correct_mean(tempi) = ...
            mean(temp_correct(end+1+temp_ROI_start:end+temp_ROI_end));
    end
    %if ~isnan(temp_error)
    if length(temp_error) >= abs(temp_ROI)
        positive_negative_choiceMemory_seq_pupilSize_ROI.positive_error_mean(tempi) = ...
            mean(temp_error(end+1+temp_ROI_start:end+temp_ROI_end));
    end
end
for tempi=1:temp_length_negative
    temp_correct = positive_negative_pupilSize_seqTimeCourse.negative_correct_mean{tempi};
    temp_error = positive_negative_pupilSize_seqTimeCourse.negative_error_mean{tempi};
    positive_negative_choiceMemory_seq_pupilSize_ROI.negative_correct_mean(tempi) = nan;
    positive_negative_choiceMemory_seq_pupilSize_ROI.negative_error_mean(tempi) = nan;
    %if ~isnan(temp_correct)
    if length(temp_correct) >= abs(temp_ROI)
        positive_negative_choiceMemory_seq_pupilSize_ROI.negative_correct_mean(tempi) = ...
            mean(temp_correct(end+1+temp_ROI:end));
    end
    %if ~isnan(temp_error)
    if length(temp_error) >= abs(temp_ROI)
        positive_negative_choiceMemory_seq_pupilSize_ROI.negative_error_mean(tempi) = ...
            mean(temp_error(end+1+temp_ROI:end));
    end
end

if if_gAcc_errorBar == 1
    gAcc_noChoice_merged;
    fileMerged_size = size(gAcc_noChoice_merged, 1);
    gAcc_noChoice_merged_SEM = cell(1, pointKindsNum);
    gAcc_noChoice_merged_inSeq = cell(1, pointKindsNum);
    for tempi=1:pointKindsNum
        gAcc_noChoice_merged_inSeq{tempi} = cell(1, numSeq(tempi));
        for tempj=1:numSeq(tempi)
            temp_merged = [];
            for temp_mergedIndex=1:fileMerged_size
                temp_merged = [temp_merged gAcc_noChoice_merged{temp_mergedIndex, tempi}(tempj)];
            end
            temp_valid = ~isnan(temp_merged);
            gAcc_noChoice_merged_SEM{tempi}(tempj) = std(temp_merged(temp_valid))...
                /sqrt(sum(temp_valid));
            
            gAcc_noChoice_merged_inSeq{tempi}{tempj} = temp_merged;
        end
    end
    
    % some gAcc_choice_merged_SEM could be NaN, beacause offloading rate is
    % 100%
    gAcc_choice_merged;
    gAcc_choice_merged_SEM = cell(1, pointKindsNum);
    gAcc_choice_merged_inSeq = cell(1, pointKindsNum);
    for tempi=1:pointKindsNum
        gAcc_choice_merged_inSeq{tempi} = cell(1, numSeq(tempi));
        for tempj=1:numSeq(tempi)
            temp_merged = [];
            for temp_mergedIndex=1:fileMerged_size
                temp_merged = [temp_merged gAcc_choice_merged{temp_mergedIndex, tempi}(tempj)];
            end
            temp_valid = ~isnan(temp_merged);
            gAcc_choice_merged_SEM{tempi}(tempj) = std(temp_merged(temp_valid))...
                /sqrt(sum(temp_valid));
            
            gAcc_choice_merged_inSeq{tempi}{tempj} = temp_merged;
        end
    end
    
    p_choice_noChoice = cell(1, pointKindsNum);
    if_signifi_p_choice_noChoice = cell(1, pointKindsNum);
    gAcc_choice_merged_inSeq;
    gAcc_noChoice_merged_inSeq;
    for tempi=1:pointKindsNum
        p_choice_noChoice{tempi} = zeros(1, numSeq(tempi));
        if_signifi_p_choice_noChoice{tempi} = false(1, numSeq(tempi));
        for tempj=1:numSeq(tempi)
            %temp_gAcc_choice = gAcc_choice_merged_inSeq{tempi}{tempj};
            temp_valid = ~isnan(gAcc_choice_merged_inSeq{tempi}{tempj});
            temp_gAcc_choice = gAcc_choice_merged_inSeq{tempi}{tempj}(temp_valid);
            temp_gAcc_noChoice = gAcc_noChoice_merged_inSeq{tempi}{tempj};
            [h_choice_noChoice,p_choice_noChoice{tempi}(tempj)]=...
                ttest2(temp_gAcc_choice, temp_gAcc_noChoice);
            if p_choice_noChoice{tempi}(tempj) < 0.05
                if_signifi_p_choice_noChoice{tempi}(tempj) = true;
            end
        end
    end
    
end

% valid_index_noChoice = offloadingProb_inOne <= offloading_threshold;
% valid_index_choice = offloadingProb_inOne <= offloading_threshold;

if if_additionFilter_lowGAcc == 0
    valid_index = offloadingProb_inOne <= offloading_threshold;
elseif if_additionFilter_lowGAcc == 1
    valid_index = (offloadingProb_inOne <= offloading_threshold ...
        & gAcc_noChoice_collapsed_inOne > 0);
end
% valid_index = true(1, sum(numSeq));

valid_index = valid_index & (sum(responseMatrix_noChoice, 2)~=0)' & (sum(responseMatrix_choice, 2)~=0)';

valid_index_testingSet = offloadingProb_inOne <= 0.9999;

%% Get data from each file
initialization = 0;
subject_name = cell(1, fileSize);

if if_monkey_D0_Z1 == 1
    trial_num_valid = 63*3*2;
end

for i=1:fileSize
    temp_load = eval(['MAT_file_load.file',sprintf('%d',i)]);
    
    seq_length_rangeHead = temp_load.basic_para.seq_length_rangeHead;
    seq_length_rangeTail = temp_load.basic_para.seq_length_rangeTail;
    pointKindsNum = seq_length_rangeTail-seq_length_rangeHead+1;
    
    ifOffloadingInSeqLength = cell(1, pointKindsNum);
    %NumOffloadingInSeqLength = zeros(1, pointKindsNum);
    
    if initialization == 0
        NumChoiceInSeqLength = zeros(fileSize, pointKindsNum);
        NumOffloadingInSeqLength = zeros(fileSize, pointKindsNum);
        NumChoiceInSeqLength_afterError = zeros(fileSize, pointKindsNum);
        NumOffloadingInSeqLength_afterError = zeros(fileSize, pointKindsNum);
        NumInternalInSeqLength_afterError = zeros(fileSize, pointKindsNum);
        NumChoiceInSeqLength_afterCorrect = zeros(fileSize, pointKindsNum);
        NumOffloadingInSeqLength_afterCorrect = zeros(fileSize, pointKindsNum);
        NumInternalInSeqLength_afterCorrect = zeros(fileSize, pointKindsNum);
        
        ProbOffloadingInSeqLength = zeros(fileSize, pointKindsNum);
        ProbOffloadingInSeqLength_afterError = zeros(fileSize, pointKindsNum);
        ProbOffloadingInSeqLength_afterCorrect = zeros(fileSize, pointKindsNum);
        
        endingHold_correct = cell(fileSize, pointKindsNum);
        endingHold_error = cell(fileSize, pointKindsNum);
        
        internalNoChoiceAccuracy = zeros(fileSize, pointKindsNum);
        offloadNoChoiceAccuracy = zeros(fileSize, pointKindsNum);
        allNoChoiceAccuracy = zeros(fileSize, pointKindsNum);
        internalChoiceAccuracy = zeros(fileSize, pointKindsNum);
        offloadChoiceAccuracy = zeros(fileSize, pointKindsNum);
        allChoiceAccuracy = zeros(fileSize, pointKindsNum);
        
        internalChoiceAccuracy_afterError = zeros(fileSize, pointKindsNum);
        internalChoiceAccuracy_afterCorrect = zeros(fileSize, pointKindsNum);
        internalNoChoiceAccuracy_afterError = zeros(fileSize, pointKindsNum);
        internalNoChoiceAccuracy_afterCorrect = zeros(fileSize, pointKindsNum);
        
        
        selecting_RT = cell(fileSize, pointKindsNum);
        meanSelecting_RT = zeros(fileSize, pointKindsNum);
        
        initialization = 1;
    end
    
    
    trial_num = temp_load.trial_para.trial_count;
    
    if if_monkey_D0_Z1 == 1
        %trial_num = min([trial_num,trial_num_valid]);
    end
    
    cumulativeError = temp_load.trial_para.cumulativeError;
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
    
    trialNum_internalNoChoice = 0;
    trialNum_offloadNoChoice = 0;
    trialNum_allNoChoice = 0;
    trialNum_internalChoice = 0;
    trialNum_offloadChoice = 0;
    trialNum_allChoice = 0;
    trialNum_internalChoice_newSeq = 0;
    trialNum_offloadChoice_newSeq = 0;
    trialNum_allChoice_newSeq = 0;
    
    ifSelectOffloading = temp_load.trial_para.ifSelectOffloading;
    for j=1:trial_num
        if choiceCondition_flag(j) == 0
            trialNum_internalNoChoice = trialNum_internalNoChoice + 1;
            trialNum_allNoChoice = trialNum_allNoChoice + 1;
        elseif choiceCondition_flag(j) == 1
            trialNum_offloadNoChoice = trialNum_offloadNoChoice + 1;
            trialNum_allNoChoice = trialNum_allNoChoice + 1;
        elseif choiceCondition_flag(j) == 2
            trialNum_allChoice = trialNum_allChoice + 1;
            if ifSelectOffloading(j) == -1
                trialNum_internalChoice = trialNum_internalChoice + 1;
            elseif ifSelectOffloading(j) == 1
                trialNum_offloadChoice = trialNum_offloadChoice + 1;
            end
            if selectFlag_newSeq(j) == 1
                trialNum_allChoice_newSeq = trialNum_allChoice_newSeq + 1;
                if ifSelectOffloading(j) == -1
                    trialNum_internalChoice_newSeq = trialNum_internalChoice_newSeq + 1;
                elseif ifSelectOffloading(j) == 1
                    trialNum_offloadChoice_newSeq = trialNum_offloadChoice_newSeq + 1;
                end
            end
        end
    end
    
    
    for j=1:pointKindsNum
        tempIndex = find(temp_load.trial_para.seq_length(1: trial_num) == (j+seq_length_rangeHead-1));
        
        tempIndex_internalNoChoice = tempIndex(choiceCondition_flag(tempIndex) == 0);
        tempIndex_offloadNoChoice = tempIndex(choiceCondition_flag(tempIndex) == 1);
        tempIndex_allNoChoice = sort([tempIndex_internalNoChoice tempIndex_offloadNoChoice]);
        
        tempIndex_allChoice = tempIndex(choiceCondition_flag(tempIndex) == 2);
        tempIndex_internalChoice = tempIndex_allChoice(ifSelectOffloading(tempIndex_allChoice) == -1);
        tempIndex_offloadChoice = tempIndex_allChoice(ifSelectOffloading(tempIndex_allChoice) == 1);
        
        tempIndex_allChoice_newSeq = tempIndex_allChoice(selectFlag_newSeq(tempIndex_allChoice) == 1);
        tempIndex_internalChoice_newSeq = tempIndex_internalChoice(selectFlag_newSeq(tempIndex_internalChoice) == 1);
        tempIndex_offloadChoice_newSeq =  tempIndex_offloadChoice(selectFlag_newSeq(tempIndex_offloadChoice) == 1);
        
        if isempty(tempIndex_allChoice_newSeq) == 0
            if tempIndex_allChoice_newSeq(1) == 1
                tempIndex_allChoice_newSeq = tempIndex_allChoice_newSeq(2:end);
            end
        end
        
        
        if tempIndex_internalNoChoice(1) == 1
            tempIndex_internalNoChoice_2 = tempIndex_internalNoChoice(2:end);
        else
            tempIndex_internalNoChoice_2 = tempIndex_internalNoChoice;
        end
        if tempIndex_allChoice_newSeq(1) == 1
            tempIndex_allChoice_newSeq_2 = tempIndex_allChoice_newSeq(2:end);
        else
            tempIndex_allChoice_newSeq_2 = tempIndex_allChoice_newSeq;
        end
        
        tempIndex_internalNoChoice_afterError = tempIndex_internalNoChoice_2(isCorrect(tempIndex_internalNoChoice_2-1) == -1);
        tempIndex_internalNoChoice_afterCorrect = tempIndex_internalNoChoice_2(isCorrect(tempIndex_internalNoChoice_2-1) == 1);
        
        tempIndex_allChoice_afterError_newSeq = tempIndex_allChoice_newSeq_2(isCorrect(tempIndex_allChoice_newSeq_2-1) == -1);
        tempIndex_offloadChoice_afterError_newSeq = ...
            tempIndex_offloadChoice_newSeq(ismember(tempIndex_offloadChoice_newSeq, tempIndex_allChoice_afterError_newSeq));
        tempIndex_internalChoice_afterError_newSeq = ...
            tempIndex_internalChoice_newSeq(ismember(tempIndex_internalChoice_newSeq, tempIndex_allChoice_afterError_newSeq));
        
        tempIndex_allChoice_afterCorrect_newSeq = tempIndex_allChoice_newSeq(isCorrect(tempIndex_allChoice_newSeq-1) == 1);
        tempIndex_offloadChoice_afterCorrect_newSeq = ...
            tempIndex_offloadChoice_newSeq(ismember(tempIndex_offloadChoice_newSeq, tempIndex_allChoice_afterCorrect_newSeq));
        tempIndex_internalChoice_afterCorrect_newSeq = ...
            tempIndex_internalChoice_newSeq(ismember(tempIndex_internalChoice_newSeq, tempIndex_allChoice_afterCorrect_newSeq));
        
        
        NumChoiceInSeqLength(i, j) = length(tempIndex_allChoice_newSeq);
        NumOffloadingInSeqLength(i, j) = length(tempIndex_offloadChoice_newSeq);
        NumChoiceInSeqLength_afterError(i, j) = length(tempIndex_allChoice_afterError_newSeq);
        NumOffloadingInSeqLength_afterError(i, j) = length(tempIndex_offloadChoice_afterError_newSeq);
        NumInternalInSeqLength_afterError(i, j) = length(tempIndex_internalChoice_afterError_newSeq);
        NumChoiceInSeqLength_afterCorrect(i, j) = length(tempIndex_allChoice_afterCorrect_newSeq);
        NumOffloadingInSeqLength_afterCorrect(i, j) = length(tempIndex_offloadChoice_afterCorrect_newSeq);
        NumInternalInSeqLength_afterCorrect(i, j) = length(tempIndex_internalChoice_afterCorrect_newSeq);
        
        ProbOffloadingInSeqLength(i, j) = length(tempIndex_offloadChoice_newSeq)/length(tempIndex_allChoice_newSeq);
        
        ProbOffloadingInSeqLength_afterError(i, j) = length(tempIndex_offloadChoice_afterError_newSeq)/ ...
            length(tempIndex_allChoice_afterError_newSeq);
        
        ProbOffloadingInSeqLength_afterCorrect(i, j) = length(tempIndex_offloadChoice_afterCorrect_newSeq)/ ...
            length(tempIndex_allChoice_afterCorrect_newSeq);
        
        %
        internalNoChoiceAccuracy_afterError(i, j) = ...
            sum(isCorrect(tempIndex_internalNoChoice_afterError)==1)/length(tempIndex_internalNoChoice_afterError);
        
        internalNoChoiceAccuracy_afterCorrect(i, j) = ...
            sum(isCorrect(tempIndex_internalNoChoice_afterCorrect)==1)/length(tempIndex_internalNoChoice_afterCorrect);
        
        internalChoiceAccuracy_afterError(i, j) = ...
            sum(isCorrect(tempIndex_internalChoice_afterError_newSeq)==1)/length(tempIndex_internalChoice_afterError_newSeq);
        
        internalChoiceAccuracy_afterCorrect(i, j) = ...
            sum(isCorrect(tempIndex_internalChoice_afterCorrect_newSeq)==1)/length(tempIndex_internalChoice_afterCorrect_newSeq);
        %
        
        tempIndex_internal = sort([tempIndex_internalNoChoice tempIndex_internalChoice_newSeq]);
        tempIndex_internal_correct = tempIndex_internal(isCorrect(tempIndex_internal) == 1);
        tempIndex_internal_error = tempIndex_internal(isCorrect(tempIndex_internal) == -1);
        
        if length(temp_load.trial_para.endingHold_RT) < trial_num
            tempEnding = trial_num - length(temp_load.trial_para.endingHold_RT);
            
            temp_load.trial_para.endingHold_RT = [temp_load.trial_para.endingHold_RT zeros(1, tempEnding)];
        end
        
        endingHold_correct{i, j} = temp_load.trial_para.endingHold_RT(tempIndex_internal_correct);
        endingHold_error{i, j} = temp_load.trial_para.endingHold_RT(tempIndex_internal_error);
        
        internalNoChoiceAccuracy(i, j) = sum(isCorrect(tempIndex_internalNoChoice)==1)/length(tempIndex_internalNoChoice);
        offloadNoChoiceAccuracy(i, j) = sum(isCorrect(tempIndex_offloadNoChoice)==1)/length(tempIndex_offloadNoChoice);
        allNoChoiceAccuracy(i, j) = sum(isCorrect(tempIndex_allNoChoice)==1)/length(tempIndex_allNoChoice);
        
        internalChoiceAccuracy(i, j) = sum(isCorrect(tempIndex_internalChoice)==1)/length(tempIndex_internalChoice);
        offloadChoiceAccuracy(i, j) = sum(isCorrect(tempIndex_offloadChoice)==1)/length(tempIndex_offloadChoice);
        allChoiceAccuracy(i, j) = sum(isCorrect(tempIndex_allChoice)==1)/length(tempIndex_allChoice);
        
        selecting_RT{i, j} = temp_load.trial_para.selecting_RT(tempIndex_allChoice_newSeq);
        meanSelecting_RT(i, j) = sum(selecting_RT{i, j})/length(selecting_RT{i, j});
        
    end
    
    %     ProbOffloadingInSeqLength
    
    %     internalNoChoiceAccuracy
    %     offloadNoChoiceAccuracy
    %     allNoChoiceAccuracy
    %     internalChoiceAccuracy
    %     offloadChoiceAccuracy
    %     allChoiceAccuracy
    
    %     selecting_RT
    %     meanSelecting_RT
    
    file_name = cell2mat(MAT_file_name(i));
    temp_tail = strfind(file_name,'train');
    %subject_name(i) = {file_name(length(searchName)+1+5: temp_tail-1)};
    subject_name(i) = {file_name(length(searchName)+1+2: temp_tail-1)};
    
    %temp_tail2 = strfind(file_name,'train');
    monkey_name = file_name(temp_tail+4+1);
    
    if ifFigure1 == 1
        %% figure 1, every single subject
        fig1 = figure(1);
        
        set(gcf,'Position',[0 0 1920 1080]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
        
        subplot(3,fileSize+1,i);
        plot(seq_length_rangeHead: seq_length_rangeTail, internalNoChoiceAccuracy(i, :), '-o', 'Color', [0.25 0.25 0.25 0.75], 'MarkerSize', 5, 'LineWidth', 1);
        hold on
        plot(seq_length_rangeHead: seq_length_rangeTail, internalChoiceAccuracy(i, :), '-o', 'Color', [0.25 0.75 0.25 0.75], 'MarkerSize', 5, 'LineWidth', 1);
        hold on
        plot(seq_length_rangeHead: seq_length_rangeTail, offloadChoiceAccuracy(i, :), '-o', 'Color', [0.75 0.25 0.25 0.75], 'MarkerSize', 5, 'LineWidth', 1);
        hold on
        plot(seq_length_rangeHead: seq_length_rangeTail, allChoiceAccuracy(i, :), '-o', 'Color', [0.25 0.25 0.75 0.75], 'MarkerSize', 5, 'LineWidth', 1);
        
        ylim([0 1]);
        set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
        set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
        set(gca,'XTickLabelRotation',0);
        set(gca,'box','off');% 取消右、上边框
        %title(sprintf('Accuracy-%s',cell2mat(subject_name(i))));
        %title(sprintf('Acc-%s',cell2mat(subject_name(i))));
        title(sprintf('Accuracy\n%s',cell2mat(subject_name(i))))
        
        subplot(3,fileSize+1,i+fileSize+1);
        %plot(seq_length_rangeHead: seq_length_rangeTail, ProbOffloadingInSeqLength(i, :),'--o', 'LineWidth', 2);
        plot(seq_length_rangeHead: seq_length_rangeTail, ProbOffloadingInSeqLength(i, :),'-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'MarkerSize', 5);
        ylim([0 1]);
        set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
        set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
        set(gca,'XTickLabelRotation',0);
        set(gca,'box','off');% 取消右、上边框
        %title(sprintf('ProbOffload-%s',cell2mat(subject_name(i))));
        %     tempname = cell2mat(subject_name(i));
        %     tempa = tempname(1:2);
        %     tempb = tempname(end-1:end);
        %     title(sprintf('rProb\n%s',[tempa tempb]));
        %title(sprintf('rProb\n%s',cell2mat(subject_name(i))))
        title(sprintf('offload\n%s',cell2mat(subject_name(i))))
        
        
        subplot(3,fileSize+1,i+(fileSize+1)*2);
        bar(seq_length_rangeHead: seq_length_rangeTail, meanSelecting_RT(i, :), ...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        for k = 1: seq_length_rangeTail-seq_length_rangeHead+1
            scatter((k-1+seq_length_rangeHead)*ones(1, length(selecting_RT{i, k})), selecting_RT{i, k}, ...
                6, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
        end
        ylim([0 4]);
        set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
        set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
        set(gca,'XTickLabelRotation',0);
        set(gca,'box','off');% 取消右、上边框
        %title(sprintf('selectingRT-%s',cell2mat(subject_name(i))));
        title(sprintf('RT-%s',cell2mat(subject_name(i))));
        
    end
end

%% Figure 1, average analysis
if ifFigure1 == 1
    subplot(3,fileSize+1,fileSize+1);
    
    plot(seq_length_rangeHead: seq_length_rangeTail, mean(internalNoChoiceAccuracy), '-o', 'Color', [0.25 0.25 0.25 0.75], 'MarkerSize', 5, 'LineWidth', 1);
    hold on
    plot(seq_length_rangeHead: seq_length_rangeTail, mean(internalChoiceAccuracy), '-o', 'Color', [0.25 0.75 0.25 0.75], 'MarkerSize', 5, 'LineWidth', 1);
    hold on
    plot(seq_length_rangeHead: seq_length_rangeTail, mean(offloadChoiceAccuracy), '-o', 'Color', [0.75 0.25 0.25 0.75], 'MarkerSize', 5, 'LineWidth', 1);
    hold on
    plot(seq_length_rangeHead: seq_length_rangeTail, mean(allChoiceAccuracy), '-o', 'Color', [0.25 0.25 0.75 0.75], 'MarkerSize', 5, 'LineWidth', 1);
    
    legend('internal|noChoice',...%'offload|noChoice',...
        'internal|Choice','offload|Choice','choice',...
        'Location','southwest')
    
    
    ylim([0 1]);
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'XTickLabelRotation',0);
    set(gca,'box','off');% 取消右、上边框
    %title(sprintf('AverageAccuracy'));
    %title(sprintf('Acc-mean'));
    title(sprintf('Accuracy\nmean'));
    
    subplot(3,fileSize+1,2*(fileSize+1));
    % plot(seq_length_rangeHead: seq_length_rangeTail, mean(ProbOffloadingInSeqLength),'--o', 'LineWidth', 2);
    plot(seq_length_rangeHead: seq_length_rangeTail, mean(ProbOffloadingInSeqLength),'-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'MarkerSize', 5);
    
    
    
    ylim([0 1]);
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'XTickLabelRotation',0);
    set(gca,'box','off');% 取消右、上边框
    %title(sprintf('AverageProbOffload'));
    %title(sprintf('rProb\nmean'));
    title(sprintf('offload\nmean'));
    
    subplot(3,fileSize+1,3*(fileSize+1));
    bar(seq_length_rangeHead: seq_length_rangeTail, mean(meanSelecting_RT), ...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
    ylim([0 4]);
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'XTickLabelRotation',0);
    set(gca,'box','off');% 取消右、上边框
    %title(sprintf('AverageSelectingRT'));
    title(sprintf('RT-mean'));
    
    currentFigName = ['rProbAccRough', '_'];
    % to generate a unique file name for saving figure
    fileName_fig1 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(1) == 1
        exportgraphics(fig1, fileName_fig1, 'BackgroundColor', 'none');
    end
end



%% Figure2
if ifFigure2 == 1
    fig2 = figure(2);
    set(gcf,'Position',[0 50 1100 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    rewardRatio = temp_load.basic_para.InternalMemoryWaterRatio;
    totalReward = rewardRatio*1;
    % totalReward = rewardRatio*1 + 1*1;
    % totalReward = rewardRatio*(1-mean(ProbOffloadingInSeqLength)) + 1*mean(ProbOffloadingInSeqLength);
    % totalReward = mean(rewardRatio*internalNoChoiceAccuracy.*(1-ProbOffloadingInSeqLength) + 1*offloadChoiceAccuracy.*ProbOffloadingInSeqLength);
    % totalReward = mean(rewardRatio*internalNoChoiceAccuracy);
    % optimal_rProb = mean(offloadChoiceAccuracy)-rewardRatio*mean(internalNoChoiceAccuracy) > 0;
    % totalReward = rewardRatio*mean(internalNoChoiceAccuracy).*(1-optimal_rProb) +...
    %     1*mean(offloadChoiceAccuracy).*optimal_rProb;
    
    %     %% Plot average accuracy
    %     % subplot(2,2,1);
    %     nexttile
    %     errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)), '-o', 'Color', [0.25 0.25 0.25 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
    %     hold on
    %     errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(internalChoiceAccuracy),std(internalChoiceAccuracy)./sqrt(size(internalChoiceAccuracy, 1)), '-o', 'Color', [0.25 0.75 0.25 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
    %     hold on
    %     errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(offloadChoiceAccuracy),std(offloadChoiceAccuracy)./sqrt(size(offloadChoiceAccuracy, 1)), '-o', 'Color', [0.75 0.25 0.25 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
    %     hold on
    %     errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(allChoiceAccuracy),std(allChoiceAccuracy)./sqrt(size(allChoiceAccuracy, 1)), '-o', 'Color', [0.25 0.25 0.75 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
    %     hold on
    %     errorbar(seq_length_rangeHead: seq_length_rangeTail,...
    %         mean((rewardRatio*internalChoiceAccuracy.*(1-ProbOffloadingInSeqLength) + 1*offloadChoiceAccuracy.*ProbOffloadingInSeqLength)./totalReward), ...
    %         std((rewardRatio*internalChoiceAccuracy.*(1-ProbOffloadingInSeqLength) + 1*offloadChoiceAccuracy.*ProbOffloadingInSeqLength)./totalReward)./sqrt(size(allChoiceAccuracy, 1)), '-o', 'Color', [0.75 0.25 0.75 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
    %
    %     % legend('internal|noChoice',...%'offload|noChoice',...
    %     %     'internal|Choice','offload|Choice','choice',...
    %     %     'Location','southwest','fontsize',13)
    %     legend('internal|noChoice',...
    %         'internal|Choice','offload|Choice','choice','rewardAccuracy',...
    %         'Location','southwest','fontsize',13)
    %
    %
    %     ylim([0 1]);
    %     set(gca, 'FontSize', 20)
    %     set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    %     set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    %     set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    %     set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    %     set(gca,'box','off');% 取消右、上边框
    %     % set(gca, 'position', [0.1 0.1 0.4 0.4]);
    %     xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
    %     ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
    %     title(sprintf('Average accuracy'));
    %     % title('\fontsize{20}AverageAccuracy');
    
    
    
    %% Plot g|choice - g|noChoice
    % subplot(2,2,2);
    nexttile
    temp_delta = internalChoiceAccuracy-internalNoChoiceAccuracy;
    temp_mean = mean(temp_delta);
    temp_SEM = std(temp_delta)./sqrt(fileSize);
    errorbar(seq_length_rangeHead: seq_length_rangeTail, temp_mean,temp_SEM, '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'CapSize', 12, 'MarkerSize', 5);
    
    
    %plot(seq_length_rangeHead: seq_length_rangeTail, mean(internalChoiceAccuracy)-mean(internalNoChoiceAccuracy), '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'MarkerSize', 5);
    %hold on
    %plot([seq_length_rangeHead-1 seq_length_rangeHead: seq_length_rangeTail seq_length_rangeTail+1]...
    %    , zeros(1, pointKindsNum+2), '--');
    
    p_t = zeros(1, pointKindsNum);
    p = p_t;
    for j=1:pointKindsNum
        %paired-sample t-test
        [h_t,p_t(j)]=ttest(internalChoiceAccuracy(:,j), internalNoChoiceAccuracy(:,j));
        
        p(j) = p_t(j);
        %Bonferroni Correction
        if p(j) < 0.05/pointKindsNum(1)
            tempTxt = sprintf('*');
            if p(j) < 0.01/pointKindsNum(1)
                tempTxt = sprintf('**');
            elseif p(j) < 0.001/pointKindsNum(1)
                tempTxt = sprintf('***');
            end
            text(j+seq_length_rangeHead-1,(mean(internalChoiceAccuracy(:,j))-mean(internalNoChoiceAccuracy(:,j)))+0.05,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
                'HorizontalAlignment','center');
        end
        
        
    end
    
    ylim([0 0.5]);
    set(gca, 'FontSize', 20)
    set(gca,'YTick',[0:0.1:0.5]);%设置要显示坐标刻度的范围
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'box','off');% 取消右、上边框
    % title('\fontsize{20}{\color{green}internal|Choice} - internal|noChoice');
    xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
    title('\fontsize{20}{\color[rgb]{0.25,0.75,0.25}internal|Choice} - internal|noChoice');
    
    % title(['\fontsize{16}black {\color{magenta}magenta '...
    % '\color[rgb]{0 .5 .5}teal \color{red}red} black again'])
    
    %% Plot average offloading rate
    x = seq_length_rangeHead: seq_length_rangeTail;
    y = mean(ProbOffloadingInSeqLength);
    
    myfunc = inline('1./(beta(1)+beta(2).*exp(-x))','beta','x');%#ok<DINLN> %三个参数分别为：函数模型(注意霿要使用点除和点乘)，待定系数，自变釿
    
    beta0 = [0.2,0.2]';%待定系数的预估忿
    beta = nlinfit(x,y,myfunc,beta0);
    
    
    x = seq_length_rangeHead:0.1:seq_length_rangeTail;
    % y = 1./(beta(1)+beta(2).*exp(-x));
    % plot(x,y);
    ProbOffloadingInSeqLength_fit = 1./(beta(1)+beta(2).*exp(-x));
    
    % subplot(2,2,3);
    nexttile
    %errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),'-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'CapSize', 12, 'MarkerSize', 5);
    
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(ProbOffloadingInSeqLength_afterError),std(ProbOffloadingInSeqLength_afterError)./sqrt(length(ProbOffloadingInSeqLength_afterError)),'-o', 'LineWidth', 2, 'CapSize', 15);
    hold on
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),'-o', 'LineWidth', 2, 'CapSize', 15, 'Color', [0.7 0.7 0.7]);
    hold on
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(ProbOffloadingInSeqLength_afterCorrect),std(ProbOffloadingInSeqLength_afterCorrect)./sqrt(length(ProbOffloadingInSeqLength_afterCorrect)),'-o', 'LineWidth', 2, 'CapSize', 15);
    
    % errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength),'--o', 'LineWidth', 2, 'CapSize', 20);
    % hold on
    % errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(ProbOffloadingInSeqLength_afterError),std(ProbOffloadingInSeqLength_afterError),'--o', 'LineWidth', 2, 'CapSize', 20);
    
    for j=1:pointKindsNum
        %Two-sample t-test
        [h_t2(j),p_t2(j)]=ttest2(ProbOffloadingInSeqLength_afterCorrect(:,j), ProbOffloadingInSeqLength_afterError(:,j));
        %paired-sample t-test
        [h_t(j),p_t(j)]=ttest(ProbOffloadingInSeqLength_afterCorrect(:,j), ProbOffloadingInSeqLength_afterError(:,j));
        
        %fisher's exact test
        A = sum(NumChoiceInSeqLength_afterError(:,j)) - sum(NumOffloadingInSeqLength_afterError(:,j));
        B = sum(NumOffloadingInSeqLength_afterError(:,j));
        C = sum(NumChoiceInSeqLength_afterCorrect(:,j)) - sum(NumOffloadingInSeqLength_afterCorrect(:,j));
        D = sum(NumOffloadingInSeqLength_afterCorrect(:,j));
        [h_f(j),p_f(j),stats_f(j)] = fishertest([A B;C D]);
        
        %chi-squared test
        [p_chi(j), Q_chi(j)]= chi2test([A B;C D]);
        
        
        p(j) = p_t(j);
        %Bonferroni Correction
        if p(j) < 0.05/pointKindsNum(1)
            tempTxt = sprintf('*');
            if p(j) < 0.01/pointKindsNum(1)
                tempTxt = sprintf('**');
            elseif p(j) < 0.001/pointKindsNum(1)
                tempTxt = sprintf('***');
            end
            %             text(j+seq_length_rangeHead-1,mean(ProbOffloadingInSeqLength_afterError(:,j))+0.2,tempTxt,'Color','black','FontSize',50,'FontWeight','bold',...
            %                 'HorizontalAlignment','center');
            text(j+seq_length_rangeHead-1,mean(ProbOffloadingInSeqLength_afterError(:,j))+0.1,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
                'HorizontalAlignment','center');
        end
        
        
    end
    
    
    % hold on
    % plot(seq_length_rangeHead:0.1: seq_length_rangeTail, ProbOffloadingInSeqLength_fit, '-k', 'LineWidth', 3);
    % legend('data','fit data','Location','southwest')
    %legend('probOffload|afterError','probOffload','probOffload|afterCorrect','Location','southeast')
    %legend('probOffload','Location','southeast')
    legend('offloadingRate|afterError','offloadingRate','offloadingRate|afterCorrect','Location','southeast','fontsize',10)
    %legend('offloading rate','Location','northwest')
    
    ylim([0 1]);
    set(gca, 'FontSize', 20)
    set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'box','off');% 取消右、上边框
    xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');
    title(sprintf('Average offloading rate'));
    %title(sprintf('AverageProbOffload'));
    % title('\fontsize{20}AverageProbOffload');
    
    
    %% Plot g|choice with history
    nexttile
    
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(internalChoiceAccuracy_afterError),std(internalChoiceAccuracy_afterError)./sqrt(length(internalChoiceAccuracy_afterError)),'-o', 'LineWidth', 2, 'CapSize', 15);
    hold on
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(internalChoiceAccuracy),std(internalChoiceAccuracy)./sqrt(size(internalChoiceAccuracy, 1)),'-o', 'LineWidth', 2, 'CapSize', 15, 'Color', [0.7 0.7 0.7]);
    hold on
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(internalChoiceAccuracy_afterCorrect),std(internalChoiceAccuracy_afterCorrect)./sqrt(length(internalChoiceAccuracy_afterCorrect)),'-o', 'LineWidth', 2, 'CapSize', 15);
    
    for j=1:pointKindsNum
        %Two-sample t-test
        [h_t2(j),p_t2(j)]=ttest2(internalChoiceAccuracy_afterCorrect(:,j), internalChoiceAccuracy_afterError(:,j));
        %paired-sample t-test
        [h_t(j),p_t(j)]=ttest(internalChoiceAccuracy_afterCorrect(:,j), internalChoiceAccuracy_afterError(:,j));
        
        p(j) = p_t(j);
        %Bonferroni Correction
        if p(j) < 0.05/pointKindsNum(1)
            tempTxt = sprintf('*');
            if p(j) < 0.01/pointKindsNum(1)
                tempTxt = sprintf('**');
            elseif p(j) < 0.001/pointKindsNum(1)
                tempTxt = sprintf('***');
            end
            text(j+seq_length_rangeHead-1,mean(internalChoiceAccuracy_afterError(:,j))+0.1,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
                'HorizontalAlignment','center');
        end
        
    end
    
    legend('accuracy|afterError','accuracy','accuracy|afterCorrect','Location','southeast','fontsize',10)
    
    ylim([0 1]);
    set(gca, 'FontSize', 20)
    set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'box','off');% 取消右、上边框
    xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
    title(sprintf('Accuracy|choice'));
    
    
    %% Plot g|noChoice with history
    nexttile
    
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(internalNoChoiceAccuracy_afterError),std(internalNoChoiceAccuracy_afterError)./sqrt(length(internalNoChoiceAccuracy_afterError)),'-o', 'LineWidth', 2, 'CapSize', 15);
    hold on
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)),'-o', 'LineWidth', 2, 'CapSize', 15, 'Color', [0.7 0.7 0.7]);
    hold on
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(internalNoChoiceAccuracy_afterCorrect),std(internalNoChoiceAccuracy_afterCorrect)./sqrt(length(internalNoChoiceAccuracy_afterCorrect)),'-o', 'LineWidth', 2, 'CapSize', 15);
    
    for j=1:pointKindsNum
        %Two-sample t-test
        [h_t2(j),p_t2(j)]=ttest2(internalNoChoiceAccuracy_afterCorrect(:,j), internalNoChoiceAccuracy_afterError(:,j));
        %paired-sample t-test
        [h_t(j),p_t(j)]=ttest(internalNoChoiceAccuracy_afterCorrect(:,j), internalNoChoiceAccuracy_afterError(:,j));
        
        p(j) = p_t(j);
        %Bonferroni Correction
        if p(j) < 0.05/pointKindsNum(1)
            tempTxt = sprintf('*');
            if p(j) < 0.01/pointKindsNum(1)
                tempTxt = sprintf('**');
            elseif p(j) < 0.001/pointKindsNum(1)
                tempTxt = sprintf('***');
            end
            text(j+seq_length_rangeHead-1,mean(internalNoChoiceAccuracy_afterCorrect(:,j))+0.1,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
                'HorizontalAlignment','center');
        end
        
    end
    
    legend('accuracy|afterError','accuracy','accuracy|afterCorrect','Location','southeast','fontsize',10)
    
    ylim([0 1]);
    set(gca, 'FontSize', 20)
    set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'box','off');% 取消右、上边框
    xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
    title(sprintf('Accuracy|noChoice'));
    
    
    %     %% Plot g|choice - g|noChoice with history
    %     nexttile
    %
    %     deltaGACC = internalChoiceAccuracy - internalNoChoiceAccuracy;
    %     deltaGACC_afterError = internalChoiceAccuracy_afterError - internalNoChoiceAccuracy_afterError;
    %     deltaGACC_afterCorrect = internalChoiceAccuracy_afterCorrect - internalNoChoiceAccuracy_afterCorrect;
    %
    %     errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(deltaGACC_afterError),std(deltaGACC_afterError)./sqrt(length(deltaGACC_afterError)),'-o', 'LineWidth', 2, 'CapSize', 15);
    %     hold on
    %     errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(deltaGACC),std(deltaGACC)./sqrt(size(deltaGACC, 1)),'-o', 'LineWidth', 2, 'CapSize', 15, 'Color', [0.7 0.7 0.7]);
    %     hold on
    %     errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(deltaGACC_afterCorrect),std(deltaGACC_afterCorrect)./sqrt(length(deltaGACC_afterCorrect)),'-o', 'LineWidth', 2, 'CapSize', 15);
    %
    %     for j=1:pointKindsNum
    %         %Two-sample t-test
    %         [h_t2(j),p_t2(j)]=ttest2(deltaGACC_afterCorrect(:,j), deltaGACC_afterError(:,j));
    %         %paired-sample t-test
    %         [h_t(j),p_t(j)]=ttest(deltaGACC_afterCorrect(:,j), deltaGACC_afterError(:,j));
    %
    %         p(j) = p_t(j);
    %         %Bonferroni Correction
    %         if p(j) < 0.05/pointKindsNum(1)
    %             tempTxt = sprintf('*');
    %             if p(j) < 0.01/pointKindsNum(1)
    %                 tempTxt = sprintf('**');
    %             elseif p(j) < 0.001/pointKindsNum(1)
    %                 tempTxt = sprintf('***');
    %             end
    %             text(j+seq_length_rangeHead-1,mean(deltaGACC_afterCorrect(:,j))+0.1,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
    %                 'HorizontalAlignment','center');
    %         end
    %
    %     end
    %
    %     legend('accuracy|afterError','accuracy','accuracy|afterCorrect','Location','southeast','fontsize',10)
    %
    %     ylim([0 1]);
    %     set(gca, 'FontSize', 20)
    %     set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    %     set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    %     set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    %     set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    %     set(gca,'box','off');% 取消右、上边框
    %     xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
    %     ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
    %     title(sprintf('Accuracy|choice - no choice'));
    
    
    currentFigName = ['rProbAccRough_average', '_'];
    % to generate a unique file name for saving figure
    fileName_fig2 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(2) == 1
        exportgraphics(fig2, fileName_fig2, 'BackgroundColor', 'none');
    end
    
end

%% Figure3, plot rProb, accuracy, gAcc|choice-gAcc|noChoice in a line
if ifFigure3 == 1
    fig3 = figure(3);
    % set(gcf,'color',[231 231 230]/255);
    set(gcf,'Position',[0 50 1650 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    
    %% Plot rProb
    nexttile
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(ProbOffloadingInSeqLength),std(ProbOffloadingInSeqLength)./sqrt(size(ProbOffloadingInSeqLength, 1)),'-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'CapSize', 12, 'MarkerSize', 5);
    legend('offloading rate','Location','northwest')
    % legend('offloading rate','Location','northwest','Color',[231 231 230]/255)
    
    ylim([0 1]);
    % set(gca,'color',[231 231 230]/255);
    set(gca, 'FontSize', 20)
    set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'box','off');% 取消右、上边框
    xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');
    title(sprintf('Offloading rate'));
    
    % rewardRatio = temp_load.basic_para.InternalMemoryWaterRatio;
    % totalReward = rewardRatio*1 + 1*1;
    
    %% Plot accuracy
    nexttile
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(internalNoChoiceAccuracy),std(internalNoChoiceAccuracy)./sqrt(size(internalNoChoiceAccuracy, 1)), '-o', 'Color', [0.25 0.25 0.25 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
    hold on
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(internalChoiceAccuracy),std(internalChoiceAccuracy)./sqrt(size(internalChoiceAccuracy, 1)), '-o', 'Color', [0.25 0.75 0.25 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
    hold on
    errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(offloadChoiceAccuracy),std(offloadChoiceAccuracy)./sqrt(size(offloadChoiceAccuracy, 1)), '-o', 'Color', [0.75 0.25 0.25 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
    %     hold on
    %     errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(allChoiceAccuracy),std(allChoiceAccuracy)./sqrt(size(allChoiceAccuracy, 1)), '-o', 'Color', [0.25 0.25 0.75 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
    % hold on
    % errorbar(seq_length_rangeHead: seq_length_rangeTail, mean((rewardRatio*internalChoiceAccuracy + 1*offloadChoiceAccuracy)/totalReward),...
    %     std((rewardRatio*internalChoiceAccuracy + 1*offloadChoiceAccuracy)/totalReward)./sqrt(size(allChoiceAccuracy, 1)), '-o', 'Color', [0.75 0.25 0.75 0.75], 'LineWidth', 1, 'CapSize', 12, 'MarkerSize', 5);
    
    %     hl = legend('internal|noChoice',...
    %         'internal|Choice','offload|Choice','choice',...
    %         'Location','southwest','fontsize',13);
    hl = legend('internal|noChoice',...
        'internal|Choice','offload|Choice',...
        'Location','southwest','fontsize',13);
    
    set(hl,'box','off');
    
    
    
    % legend('internal|noChoice',...
    %     'internal|Choice','offload|Choice','choice','rewardExpect',...
    %     'Location','southwest','fontsize',13)
    ylim([0 1]);
    % set(gca,'color',[231 231 230]/255);
    set(gca, 'FontSize', 20)
    set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'box','off');% 取消右、上边框
    xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
    %ylabel('Correct rate', 'FontSize', 18, 'FontWeight', 'bold');
    % title(sprintf('Average correct rate'));
    ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
    title(sprintf('Accuracy'));
    
end

%% Calculate gAcc|choice-gAcc|noChoice
% Now choiceMinusNoChoice is calculated in a precise way, which means
% calculating gAcc|choice-gAcc|noChoice in a same seq

% choiceMinusNoChoice = internalChoiceAccuracy-internalNoChoiceAccuracy;
% errorbar(seq_length_rangeHead: seq_length_rangeTail, mean(choiceMinusNoChoice),std(choiceMinusNoChoice)./sqrt(size(choiceMinusNoChoice, 1)), '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'CapSize', 12, 'MarkerSize', 5);

choiceDivedeNoChoice = cell(1, pointKindsNum);
choiceNoChoiceContrast = cell(1, pointKindsNum);
choiceMinusNoChoiceNormalized = cell(1, pointKindsNum);

choiceMinusNoChoice = cell(1, pointKindsNum);
choiceMinusNoChoice_weighted = cell(1, pointKindsNum);
choiceMinusNoChoice_rawWeight = cell(1, pointKindsNum);
choiceMinusNoChoice_weight = cell(1, pointKindsNum);

choiceMinusNoChoice_inOne = zeros(sum(numSeq), 1);
% choiceMinusNoChoice_weighted_inOne = zeros(sum(numSeq), 1);
choiceMinusNoChoice_rawWeight_inOne = zeros(sum(numSeq), 1);
% choiceMinusNoChoice_weight_inOne = zeros(sum(numSeq), 1);

gCorrect_choice_trial_count_collapsed;
g_choice_trial_count_collapsed;
gCorrect_noChoice_trial_count_collapsed;
g_noChoice_trial_count_collapsed;
tempj = 0;
for target_seqLength=1:pointKindsNum
    choiceMinusNoChoiceNormalized{target_seqLength} = zeros(numSeq(target_seqLength), 1);
    choiceNoChoiceContrast{target_seqLength} = zeros(numSeq(target_seqLength), 1);
    choiceDivedeNoChoice{target_seqLength} = zeros(numSeq(target_seqLength), 1);
    choiceMinusNoChoice{target_seqLength} = zeros(numSeq(target_seqLength), 1);
    choiceMinusNoChoice_rawWeight{target_seqLength} = zeros(numSeq(target_seqLength), 1);
    for tempi=1:numSeq(target_seqLength)
        tempj = tempj + 1;
        
        temp_gAcc_choice = gCorrect_choice_trial_count_collapsed{target_seqLength}(tempi) / g_choice_trial_count_collapsed{target_seqLength}(tempi);
        temp_gAcc_noChoice = gCorrect_noChoice_trial_count_collapsed{target_seqLength}(tempi) / g_noChoice_trial_count_collapsed{target_seqLength}(tempi);
        
        choiceMinusNoChoiceNormalized{target_seqLength}(tempi) = ...
            (temp_gAcc_choice - temp_gAcc_noChoice) / temp_gAcc_noChoice;
        
        choiceNoChoiceContrast{target_seqLength}(tempi) = ...
            (temp_gAcc_choice - temp_gAcc_noChoice) / (temp_gAcc_choice + temp_gAcc_noChoice);
        
        choiceDivedeNoChoice{target_seqLength}(tempi) = ...
            (gCorrect_choice_trial_count_collapsed{target_seqLength}(tempi) / g_choice_trial_count_collapsed{target_seqLength}(tempi))...
            /...
            (gCorrect_noChoice_trial_count_collapsed{target_seqLength}(tempi) / g_noChoice_trial_count_collapsed{target_seqLength}(tempi));
        
        [choiceMinusNoChoice{target_seqLength}(tempi), ...
            choiceMinusNoChoice_rawWeight{target_seqLength}(tempi)] = ...
            getBalcancedChoiceMinusNoChoice_oneSeq(...
            gCorrect_choice_trial_count_collapsed{target_seqLength}(tempi), ...
            g_choice_trial_count_collapsed{target_seqLength}(tempi), ...
            gCorrect_noChoice_trial_count_collapsed{target_seqLength}(tempi), ...
            g_noChoice_trial_count_collapsed{target_seqLength}(tempi)...
            );
        choiceMinusNoChoice_inOne(tempj) = choiceMinusNoChoice{target_seqLength}(tempi);
        choiceMinusNoChoice_rawWeight_inOne(tempj) = choiceMinusNoChoice_rawWeight{target_seqLength}(tempi);
    end
    choiceMinusNoChoice_weight{target_seqLength} = choiceMinusNoChoice_rawWeight{target_seqLength}...
        / (sum(choiceMinusNoChoice_rawWeight{target_seqLength}) / numSeq(target_seqLength));
    choiceMinusNoChoice_weighted{target_seqLength} = choiceMinusNoChoice{target_seqLength} ...
        .* choiceMinusNoChoice_weight{target_seqLength};
end
choiceMinusNoChoice_weight_inOne = choiceMinusNoChoice_rawWeight_inOne ...
    / (sum(choiceMinusNoChoice_rawWeight_inOne) / sum(numSeq));
choiceMinusNoChoice_weighted_inOne = choiceMinusNoChoice_inOne ...
    .* choiceMinusNoChoice_weight_inOne;

choiceMinusNoChoice;
choiceMinusNoChoice_weight;
choiceMinusNoChoice_weighted;
choiceMinusNoChoice_inOne;
choiceMinusNoChoice_weight_inOne;
choiceMinusNoChoice_weighted_inOne;

choiceDivedeNoChoice;
choiceNoChoiceContrast;
choiceMinusNoChoiceNormalized;

% validLowBound_gAcc = 0.05;
% validLowBound_rProb = 0.2;
% offloading_threshold = 0.75;
offloadingProb_collapsed;

% if_plot_choiceMinusNoChoice_inOne = 0;

if ifFigure3 == 1 && if_plot_partialAccuracy == 0
    %% Plot gAcc|choice-gAcc|noChoice
    nexttile
    if if_plot_choiceMinusNoChoice_inOne == 1
        offloadingProb_collapsed_inOne = [];
        for target_seqLength=1:pointKindsNum
            offloadingProb_collapsed_inOne = ...
                [offloadingProb_collapsed_inOne offloadingProb_collapsed{target_seqLength}']; %#ok<*AGROW>
        end
        
        nonZero_index_inOne = offloadingProb_collapsed_inOne' <= offloading_threshold;
        %nonZero_index_inOne = gAcc_noChoice_collapsed_inOne' > validLowBound_gAcc;
        
        %choiceMinusNoChoice_valid_inOne = choiceMinusNoChoice_weighted_inOne(nonZero_index_inOne);
        choiceMinusNoChoice_valid_inOne = choiceMinusNoChoice_inOne(nonZero_index_inOne);
        tempZero_inOne = zeros(length(choiceMinusNoChoice_valid_inOne), 1);
        [h_t_inOne,p_t_inOne]=ttest(choiceMinusNoChoice_valid_inOne, tempZero_inOne);
        
        choiceMinusNoChoice_valid_inOne;
        p_t_inOne;
        %boxplot(choiceMinusNoChoice_valid_inOne);
        violinplot(choiceMinusNoChoice_valid_inOne);
        hold on
        temp_xlim = xlim;
        plot([temp_xlim(1)-0.1 temp_xlim(2)], [0 0], '--', 'Color', 'k');
        
        if p_t_inOne < 0.001
            tempTxt = sprintf('***');
        elseif p_t_inOne < 0.01
            tempTxt = sprintf('**');
        elseif p_t_inOne < 0.05
            tempTxt = sprintf('*');
        end
        
        %text(1.1,0.05,[tempTxt sprintf(' p = %.4f',p_t_inOne)], 'fontsize',16,'FontWeight','bold');
        %text(1.1,0.4,[tempTxt sprintf(' p = %.4f',p_t_inOne)], 'fontsize',16,'FontWeight','bold');
        text(1.1,0.3,[tempTxt sprintf(' p = %.4f',p_t_inOne)], 'fontsize',16,'FontWeight','bold');
        
        
        %     temp_mean = mean(choiceMinusNoChoice_valid_inOne);
        %     temp_SEM = std(choiceMinusNoChoice_valid_inOne)./sqrt(sum(nonZero_index_inOne));
        %     bar(0,mean(choiceMinusNoChoice_valid_inOne),0.6, ...
        %         'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
        %     hold on
        %     errorbar(0, temp_mean,temp_SEM, '-o', 'Color', [0 0 0], 'LineWidth', 2, 'CapSize', 12, 'MarkerSize', 5);
        
        set(gca,'XTickLabel', " ");%给坐标加标签
        set(gca, 'FontSize', 20)
        set(gca,'box','off');% 取消右、上边框
        ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
        %title('\fontsize{16}{Weighted \color[rgb]{0.25,0.75,0.25}internal|Choice} - internal|noChoice');
        title('\fontsize{16}{\color[rgb]{0.25,0.75,0.25}internal|Choice} - internal|noChoice');
        
    elseif if_plot_choiceMinusNoChoice_inOne == 0
        
        
        % temp_count = 0;
        % for tempi=1:pointKindsNum
        %     choiceMinusNoChoice_partial{tempi} = zeros(numSeq(tempi), 1);
        %     for tempj=1:numSeq(tempi)
        %         temp_count = temp_count + 1;
        %         choiceMinusNoChoice_partial{tempi}(tempj) = choiceMinusNoChoice_partial_inOne(temp_count);
        %
        %         if valid_index(temp_count) == 1
        %             choiceMinusNoChoice_partial_valid{tempi} = [choiceMinusNoChoice_partial_valid{tempi};...
        %                 choiceMinusNoChoice_partial_inOne(temp_count)];
        %         end
        %     end
        % end
        
        choiceMinusNoChoice_valid = cell(1, pointKindsNum);
        temp_mean = zeros(1, pointKindsNum);
        temp_SEM = zeros(1, pointKindsNum);
        for target_seqLength=1:pointKindsNum
            %         nonZero_index = gAcc_noChoice_collapsed{target_seqLength} > validLowBound_gAcc;
            %         nonZero_index = true(numSeq(target_seqLength), 1);
            
            %nonZero_index = offloadingProb_collapsed{target_seqLength} <= offloading_threshold;
            
            if target_seqLength == 1
                nonZero_index = valid_index(1:numSeq(1));
            else
                nonZero_index = valid_index( sum(numSeq(1:target_seqLength-1))+...
                    (1:numSeq(target_seqLength)) );
            end
            
            
            
            %         choiceMinusNoChoice_valid{target_seqLength} = ...
            %             choiceMinusNoChoice_weighted{target_seqLength}(nonZero_index);
            choiceMinusNoChoice_valid{target_seqLength} = ...
                choiceMinusNoChoice{target_seqLength}(nonZero_index);
            %         choiceMinusNoChoice_valid{target_seqLength} = ...
            %             choiceDivedeNoChoice{target_seqLength}(nonZero_index);
            %         choiceMinusNoChoice_valid{target_seqLength} = ...
            %             choiceNoChoiceContrast{target_seqLength}(nonZero_index);
            %         choiceMinusNoChoice_valid{target_seqLength} = ...
            %             choiceMinusNoChoiceNormalized{target_seqLength}(nonZero_index);
            
            %     choiceMinusNoChoice_valid{target_seqLength} = ...
            %         gAcc_choice_collapsed{target_seqLength}(nonZero_index) - gAcc_noChoice_collapsed{target_seqLength}(nonZero_index);
            temp_mean(target_seqLength) = mean(choiceMinusNoChoice_valid{target_seqLength});
            temp_SEM(target_seqLength) = ...
                std(choiceMinusNoChoice_valid{target_seqLength})./sqrt(sum(nonZero_index));
        end
        errorbar(seq_length_rangeHead: seq_length_rangeTail, temp_mean,temp_SEM, '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'CapSize', 12, 'MarkerSize', 5);
        
        hold on
        plot([seq_length_rangeHead-1 seq_length_rangeHead: seq_length_rangeTail seq_length_rangeTail+1]...
            , zeros(1, pointKindsNum+2), '--');
        %     plot([seq_length_rangeHead-1 seq_length_rangeHead: seq_length_rangeTail seq_length_rangeTail+1]...
        %         , ones(1, pointKindsNum+2), '--');
        
        choiceMinusNoChoice;
        choiceMinusNoChoice_valid;
        for j=1:pointKindsNum
            tempZero = zeros(length(choiceMinusNoChoice_valid{j}), 1);
            tempOne = ones(length(choiceMinusNoChoice_valid{j}), 1);
            %paired-sample t-test
            [h_t(j),p_t(j)]=ttest(choiceMinusNoChoice_valid{j}, tempZero);
            %[h_t(j),p_t(j)]=ttest(choiceMinusNoChoice_valid{j}, tempOne);
            % Wilcoxon rank sum test
            p_rs(j) = ...
                ranksum(choiceMinusNoChoice_valid{j}, tempZero); %#ok<*SAGROW>
            
            p(j) = p_t(j);
            %p(j) = p_rs(j);
            
            %Bonferroni Correction
            
            if p(j) < 0.001/pointKindsNum
                tempTxt = sprintf('***');
            elseif p(j) < 0.01/pointKindsNum
                tempTxt = sprintf('**');
            elseif p(j) < 0.05/pointKindsNum
                tempTxt = sprintf('*');
            end
            if p(j) < 0.05/pointKindsNum
                text(j+seq_length_rangeHead-1,mean(choiceMinusNoChoice_valid{j})-0.08,tempTxt,'Color','black','FontSize',30,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                %             text(j+seq_length_rangeHead-1,mean(choiceMinusNoChoice_valid{j})-0.02,tempTxt,'Color','black','FontSize',30,'FontWeight','bold',...
                %                 'HorizontalAlignment','center');
            end
        end
        
        
        ylim([-0.3 0.3]);
        % set(gca,'color',[231 231 230]/255);
        set(gca, 'FontSize', 20)
        set(gca,'YTick',[-0.3:0.1:0.3]);%设置要显示坐标刻度的范围
        set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
        set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        % title('\fontsize{20}{\color{green}internal|Choice} - internal|noChoice');
        xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
        %ylabel('Contrast', 'FontSize', 18, 'FontWeight', 'bold');
        %ylabel('Normalized accuracy', 'FontSize', 18, 'FontWeight', 'bold');
        %title('\fontsize{20}{Weighted \color[rgb]{0.25,0.75,0.25}internal|Choice} - internal|noChoice');
        title('\fontsize{20}{\color[rgb]{0.25,0.75,0.25}internal|Choice} - internal|noChoice');
        %title('Contrast, \fontsize{20}{\color[rgb]{0.25,0.75,0.25}internal|Choice} and internal|noChoice');
        %title('Normalized, \fontsize{20}{\color[rgb]{0.25,0.75,0.25}internal|Choice} and internal|noChoice');
        
    end
    
    currentFigName = ['rProbAccRough_average_inLine', '_'];
    % to generate a unique file name for saving figure
    fileName_fig3 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(3) == 1
        exportgraphics(fig3, fileName_fig3, 'BackgroundColor', 'none');
    end
    
end

%% Process response maxtrix
responseMatrix_noChoice;
responseMatrix_choice;
target_seqSet_inOne;
target_seqSet_extend_inOne;

responseMatrix_noChoice_n11n = responseMatrix_noChoice ./ ....
    sum(responseMatrix_noChoice, 2);
responseMatrix_choice_n11n = responseMatrix_choice ./ ....
    sum(responseMatrix_choice, 2);

responseMatrix = struct;
responseMatrix.noChoice = responseMatrix_noChoice;
responseMatrix.choice = responseMatrix_choice;
responseMatrix.noChoice_n11n = responseMatrix_noChoice_n11n;
responseMatrix.choice_n11n = responseMatrix_choice_n11n;

valid_index;


% Get relation matrix of neighbor length error
target_seqSet_inOne;
neighbor_lengthError_seqSet = cell(sum(numSeq), 1);
for tempi=1:sum(numSeq)
    temp_seq = target_seqSet_inOne{tempi};
    temp_length = length(temp_seq);
    
    temp_neighbor = cell(numFrames+1,1);
    temp_neighbor{1} = temp_seq;
    % down length
    for tempj=1:temp_length
        temp_downSeq = temp_seq(~ismember(temp_seq, temp_seq(temp_length-tempj+1)));
        temp_neighbor{1+tempj} = temp_downSeq;
    end
    
    % up length
    temp_isFill = false(1, numFrames);
    temp_isFill(temp_seq) = true;
    tempk = 0;
    for tempj=1:numFrames
        if temp_isFill(tempj) == 0
            tempk = tempk + 1;
            temp_false = false(1, numFrames);
            temp_false(tempj) = true;
            
            temp_isFill_up = temp_isFill | temp_false;
            temp_upSeq = find(temp_isFill_up==true);
            temp_neighbor{1+temp_length+tempk} = temp_upSeq;
        end
    end
    
    neighbor_lengthError_seqSet{tempi} = temp_neighbor;
end
neighbor_lengthError_seqSet_valid = neighbor_lengthError_seqSet(valid_index);
seqSet_valid = target_seqSet_inOne(valid_index);
size_stimuli = sum(valid_index);
size_response = size(target_seqSet_extend_inOne, 1);

% Get neighbor length error matrix
temp_valid_count = 0;
neighbor_lengthError_indexMatrix = false(sum(numSeq), sum(numSeq_extend));
for tempj=1:sum(numSeq)
    %neighbor_lengthError_seqSet;
    temp_seq = target_seqSet_inOne{tempj};
    temp_length = length(temp_seq);
    temp_stimuli_index = tempj;
    
    for tempk=1:(numFrames+1)
        % for length1 seq, down length would be empty, so skip the
        % empty one.
        if isempty(neighbor_lengthError_seqSet{tempj}{tempk}) == 0
            temp_seq_res = neighbor_lengthError_seqSet{tempj}{tempk};
            temp_length_res = length(temp_seq_res);
            
            temp_response_index = 0;
            for tempkk=1:sum(numSeq_extend)
                %target_seqSet_inOne{tempkk};
                if length(target_seqSet_extend_inOne{tempkk}) == temp_length_res
                    if  sum(temp_seq_res == target_seqSet_extend_inOne{tempkk}) == temp_length_res
                        temp_response_index = tempkk;
                    end
                end
            end
            temp_response_index;
            temp_stimuli_index;
            
            if temp_response_index ~= 0
                if if_lengthError_down0_up1_both2 == 0
                    if temp_length_res <= temp_length
                        neighbor_lengthError_indexMatrix(temp_stimuli_index, temp_response_index) = true;
                    end
                elseif if_lengthError_down0_up1_both2 == 1
                    if temp_length_res >= temp_length
                        neighbor_lengthError_indexMatrix(temp_stimuli_index, temp_response_index) = true;
                    end
                elseif if_lengthError_down0_up1_both2 == 2
                    neighbor_lengthError_indexMatrix(temp_stimuli_index, temp_response_index) = true;
                end
                
            end
        end
    end
    
end
neighbor_lengthError_indexMatrix;

%% Get partial accuracy
% With regard both neighbor length error and exact correct as partial correct
neighbor_lengthError_indexMatrix;
responseMatrix_noChoice;
responseMatrix_choice;

gAcc_partial_noChoice_inOne = zeros(1, sum(numSeq));
gAcc_partial_choice_inOne = zeros(1, sum(numSeq));
for tempi=1:sum(numSeq)
    temp_neighborIndex = neighbor_lengthError_indexMatrix(tempi, :); %#ok<*NOEFF>
    temp_response_noChoice = responseMatrix_noChoice(tempi, :);
    temp_response_choice = responseMatrix_choice(tempi, :);
    
    gAcc_partial_noChoice_inOne(tempi) = ...
        sum(temp_response_noChoice(temp_neighborIndex))/sum(temp_response_noChoice);
    gAcc_partial_choice_inOne(tempi) = ...
        sum(temp_response_choice(temp_neighborIndex))/sum(temp_response_choice);
end
gAcc_partial_choice_inOne;
gAcc_partial_noChoice_inOne;

choiceMinusNoChoice_partial_inOne = gAcc_partial_choice_inOne - gAcc_partial_noChoice_inOne;
% choiceMinusNoChoice_partial_inOne_valid = choiceMinusNoChoice_partial_inOne(valid_index);
choiceMinusNoChoice_partial = cell(1, pointKindsNum);
choiceMinusNoChoice_partial_valid = cell(1, pointKindsNum);
choiceMinusNoChoice;
% choiceMinusNoChoice_valid
temp_count = 0;
for tempi=1:pointKindsNum
    choiceMinusNoChoice_partial{tempi} = zeros(numSeq(tempi), 1);
    for tempj=1:numSeq(tempi)
        temp_count = temp_count + 1;
        choiceMinusNoChoice_partial{tempi}(tempj) = choiceMinusNoChoice_partial_inOne(temp_count);
        
        if valid_index(temp_count) == 1
            %if isnan(choiceMinusNoChoice_partial_inOne(temp_count))
            %    a = 1;
            %end
            choiceMinusNoChoice_partial_valid{tempi} = [choiceMinusNoChoice_partial_valid{tempi};...
                choiceMinusNoChoice_partial_inOne(temp_count)];
        end
    end
end
choiceMinusNoChoice_partial;
choiceMinusNoChoice_partial_valid;


if ifFigure3 == 1 && if_plot_partialAccuracy == 1
    %% Plot gAcc|choice-gAcc|noChoice, with partial accuracy
    nexttile
    
    temp_mean = zeros(1, pointKindsNum);
    temp_SEM = zeros(1, pointKindsNum);
    for target_seqLength=1:pointKindsNum
        temp_mean(target_seqLength) = mean(choiceMinusNoChoice_partial_valid{target_seqLength});
        temp_SEM(target_seqLength) = ...
            std(choiceMinusNoChoice_partial_valid{target_seqLength})./...
            sqrt(length(choiceMinusNoChoice_partial_valid{target_seqLength}));
    end
    errorbar(seq_length_rangeHead: seq_length_rangeTail, temp_mean,temp_SEM, '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'CapSize', 12, 'MarkerSize', 5);
    
    hold on
    plot([seq_length_rangeHead-1 seq_length_rangeHead: seq_length_rangeTail seq_length_rangeTail+1]...
        , zeros(1, pointKindsNum+2), '--');
    
    choiceMinusNoChoice_partial_valid;
    for j=1:pointKindsNum
        tempZero = zeros(length(choiceMinusNoChoice_partial_valid{j}), 1);
        tempOne = ones(length(choiceMinusNoChoice_partial_valid{j}), 1);
        %paired-sample t-test
        [h_t(j),p_t(j)]=ttest(choiceMinusNoChoice_partial_valid{j}, tempZero);
        % Wilcoxon rank sum test
        p_rs(j) = ...
            ranksum(choiceMinusNoChoice_partial_valid{j}, tempZero); %#ok<*SAGROW>
        
        p(j) = p_t(j);
        %p(j) = p_rs(j);
        
        %Bonferroni Correction
        
        if p(j) < 0.001/pointKindsNum
            tempTxt = sprintf('***');
        elseif p(j) < 0.01/pointKindsNum
            tempTxt = sprintf('**');
        elseif p(j) < 0.05/pointKindsNum
            tempTxt = sprintf('*');
        end
        if p(j) < 0.05/pointKindsNum
            text(j+seq_length_rangeHead-1,mean(choiceMinusNoChoice_partial_valid{j})-0.08,tempTxt,...
                'Color','black','FontSize',30,'FontWeight','bold',...
                'HorizontalAlignment','center');
        end
    end
    
    
    ylim([-0.3 0.3]);
    set(gca, 'FontSize', 20)
    set(gca,'YTick',[-0.3:0.1:0.3]);%设置要显示坐标刻度的范围
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'box','off');% 取消右、上边框
    xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Partial accuracy', 'FontSize', 18, 'FontWeight', 'bold');
    title('\fontsize{20}{\color[rgb]{0.25,0.75,0.25}internal|Choice} - internal|noChoice');
    
    
    
    currentFigName = ['rProbAccRough_average_inLine', '_'];
    % to generate a unique file name for saving figure
    fileName_fig3 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(3) == 1
        exportgraphics(fig3, fileName_fig3, 'BackgroundColor', 'none');
    end
end

%% Figure4, plot choice response matrix
if ifFigure4 == 1
    fig4 = figure(4);
    set(gcf,'Position',[0 50 1500 700]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    nexttile
    
    current_responseMatrix = responseMatrix.choice_n11n;
    
    % Plot heat map
    imagesc(current_responseMatrix(valid_index, :), [0 1]);
    my_gray = gray;
    my_gray = my_gray(end:-1:1,:);
    colormap(my_gray);
    hold on
    
    % Plot neighbor length errer
    if_plot_neighbor_lengthError = ...
        current_responseMatrix ~= 0 & neighbor_lengthError_indexMatrix;
    [temp_row, temp_col] = find(if_plot_neighbor_lengthError(valid_index, :)==1);
    squareSize = 130;
    if numFrames == 8
        squareSize = 20;
    end
    scatter(temp_col,temp_row,squareSize,[0 1 1]*0.7,'s');
    hold on
    
    size_stimuli;
    size_response;
    % Plot length bound in response
    for tempj=1:size_response
        if tempj > 1
            if length(target_seqSet_extend_inOne{tempj}) > ...
                    length(target_seqSet_extend_inOne{tempj-1})
                plot([tempj tempj]-0.5, [0.5 size_stimuli+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in stimuli
    for tempj=1:size_stimuli
        if tempj > 1
            if length(seqSet_valid{tempj}) > ...
                    length(seqSet_valid{tempj-1})
                plot([0.5 size_response+0.5], [tempj tempj]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    axis equal
    xlim([0.5 size_response+0.5]);
    ylim([0.5 size_stimuli+0.5]);
    x_lim = xlim;
    y_lim = ylim;
    
    set(gca, 'FontSize', 20)
    set(gca,'YDir','reverse');% Reverse the direction of y axis
    set(gca,'box','off');% 取消右、上边框
    
    set(gca,'XTick',1:size_response);
    %set(gca,'XTickLabel',seqSet_inOne_extend_inOne,'FontSize',9);%给坐标加标签
    set(gca,'YTick',1:size_stimuli);
    %set(gca,'YTickLabel',seqSet_inOne_inOne(valid_index),'FontSize',9);%给坐标加标签
    
    %xtl=get(gca,'XTickLabel');
    %ytl=get(gca,'YTickLabel');
    xtl=string(seqSet_inOne_extend_inOne);
    ytl=string(seqSet_inOne_inOne(valid_index));
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    xtext_xp=xt;
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(yt))-0.75;
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    tempFontSize = 9;
    if numFrames == 8
        tempFontSize = 3;
    end
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',90,'fontsize',tempFontSize);%9
    text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',tempFontSize);
    % 取消原始ticklabel
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    xlabel('Response seq','Position',[mean(x_lim) max(y_lim)+2],'FontSize',18,'FontWeight','bold');
    ylabel('Stimuli seq','Position',[min(x_lim)-2.5 mean(y_lim)],'FontSize',18,'FontWeight','bold');
    title('\fontsize{20}{choice}');
    
    colorbar('Ticks', [0:0.2:1]);
end

%% Figure5, plot no choice response matrix
if ifFigure5 == 1
    fig5 = figure(5);
    set(gcf,'Position',[0 50 1500 700]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    nexttile
    
    current_responseMatrix = responseMatrix.noChoice_n11n;
    
    % Plot heat map
    imagesc(current_responseMatrix(valid_index, :), [0 1]);
    my_gray = gray;
    my_gray = my_gray(end:-1:1,:);
    colormap(my_gray);
    hold on
    
    % Plot neighbor length errer
    if_plot_neighbor_lengthError = ...
        current_responseMatrix ~= 0 & neighbor_lengthError_indexMatrix;
    [temp_row, temp_col] = find(if_plot_neighbor_lengthError(valid_index, :)==1);
    squareSize = 130;
    if numFrames == 8
        squareSize = 20;
    end
    scatter(temp_col,temp_row,squareSize,[0 1 1]*0.7,'s');
    hold on
    
    size_stimuli;
    size_response;
    % Plot length bound in response
    for tempj=1:size_response
        if tempj > 1
            if length(target_seqSet_extend_inOne{tempj}) > ...
                    length(target_seqSet_extend_inOne{tempj-1})
                plot([tempj tempj]-0.5, [0.5 size_stimuli+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in stimuli
    for tempj=1:size_stimuli
        if tempj > 1
            if length(seqSet_valid{tempj}) > ...
                    length(seqSet_valid{tempj-1})
                plot([0.5 size_response+0.5], [tempj tempj]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    axis equal
    xlim([0.5 size_response+0.5]);
    ylim([0.5 size_stimuli+0.5]);
    x_lim = xlim;
    y_lim = ylim;
    
    set(gca, 'FontSize', 20)
    set(gca,'YDir','reverse');% Reverse the direction of y axis
    set(gca,'box','off');% 取消右、上边框
    
    set(gca,'XTick',1:size_response);
    %set(gca,'XTickLabel',seqSet_inOne_extend_inOne,'FontSize',9);%给坐标加标签
    set(gca,'YTick',1:size_stimuli);
    %set(gca,'YTickLabel',seqSet_inOne_inOne(valid_index),'FontSize',9);%给坐标加标签
    
    %xtl=get(gca,'XTickLabel');
    %ytl=get(gca,'YTickLabel');
    xtl=string(seqSet_inOne_extend_inOne);
    ytl=string(seqSet_inOne_inOne(valid_index));
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    xtext_xp=xt;
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(yt))-0.75;
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    tempFontSize = 9;
    if numFrames == 8
        tempFontSize = 3;
    end
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',90,'fontsize',tempFontSize);%9
    text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',tempFontSize);
    % 取消原始ticklabel
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    xlabel('Response seq','Position',[mean(x_lim) max(y_lim)+2],'FontSize',18,'FontWeight','bold');
    ylabel('Stimuli seq','Position',[min(x_lim)-2.5 mean(y_lim)],'FontSize',18,'FontWeight','bold');
    title('\fontsize{20}{no choice}');
    
    
    colorbar('Ticks', [0:0.2:1]);
end

%% Figure6, plot response matrix difference between choice and no choice
if ifFigure6 == 1
    fig6 = figure(6);
    set(gcf,'Position',[0 50 1500 700]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    nexttile
    
    responseMatrix.choice_n11n;
    responseMatrix.noChoice_n11n;
    error_responseMatrix = responseMatrix.choice_n11n-responseMatrix.noChoice_n11n;
    
    %size_stimuli = size(target_seqSet_inOne, 1);
    %size_stimuli = sum(valid_index);
    %size_response = size(target_seqSet_extend_inOne, 1);
    
    % Plot heat map
    imagesc(error_responseMatrix(valid_index, :));
    %     my_gray = gray;
    %     my_gray = my_gray(end:-1:1,:);
    %     colormap(my_gray);
    
    temp_max = max(max(error_responseMatrix(valid_index, :)));
    temp_min = min(min(error_responseMatrix(valid_index, :)));
    
    my_color = cool;
    %my_color(ceil((0-temp_min)/(temp_max-temp_min)*256), :) = [1 1 1];
    
    colormap(my_color);
    %colormap(jet);
    hold on
    
    % Plot neighbor length errer
    if_plot_neighbor_lengthError = ...
        error_responseMatrix ~= 0 & neighbor_lengthError_indexMatrix;
    [temp_row, temp_col] = find(if_plot_neighbor_lengthError(valid_index, :)==1);
    %scatter(temp_col,temp_row,70,[1 1 1]*0.2,'s');
    %scatter(temp_col,temp_row,70,[1 1 1]*0.5,'s');
    
    squareSize = 130;
    if numFrames == 8
        squareSize = 20;
    end
    scatter(temp_col,temp_row,squareSize,[1 1 1]*0.2,'s');
    
    hold on
    
    size_stimuli;
    size_response;
    % Plot length bound in response
    for tempj=1:size_response
        if tempj > 1
            if length(target_seqSet_extend_inOne{tempj}) > ...
                    length(target_seqSet_extend_inOne{tempj-1})
                plot([tempj tempj]-0.5, [0.5 size_stimuli+0.5], 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    % Plot length bound in stimuli
    %seqSet_valid = target_seqSet_inOne(valid_index);
    for tempj=1:size_stimuli
        if tempj > 1
            if length(seqSet_valid{tempj}) > ...
                    length(seqSet_valid{tempj-1})
                plot([0.5 size_response+0.5], [tempj tempj]-0.5, 'color', [0.25 0.25 0.25]);
                hold on
            end
        end
    end
    
    axis equal
    xlim([0.5 size_response+0.5]);
    ylim([0.5 size_stimuli+0.5]);
    x_lim = xlim;
    y_lim = ylim;
    
    set(gca, 'FontSize', 20)
    set(gca,'YDir','reverse');% Reverse the direction of y axis
    set(gca,'box','off');% 取消右、上边框
    
    set(gca,'XTick',1:size_response);
    set(gca,'YTick',1:size_stimuli);
    
    xtl=string(seqSet_inOne_extend_inOne);
    ytl=string(seqSet_inOne_inOne(valid_index));
    xt=get(gca,'XTick');
    yt=get(gca,'YTick');
    % 设置xtext的x坐标位置
    xtext_xp=xt;
    % 设置xtext的y坐标位置
    xtext_yp=(yt(end))*ones(1,length(xt))+0.75;
    % 设置ttext的x坐标位置
    ytext_xp=(xt(1))*ones(1,length(yt))-0.75;
    % 设置ttext的y坐标位置
    ytext_yp=yt;
    % rotation，正的旋转角度代表?时针旋转，旋转轴可以由HorizontalAlignment属濧来设定＿
    % 朿3个属性忼：left，right，center
    tempFontSize = 9;
    if numFrames == 8
        tempFontSize = 3;
    end
    text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',90,'fontsize',tempFontSize);%9
    text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','right','rotation',0,'fontsize',tempFontSize);
    % 取消原始ticklabel
    set(gca,'xticklabel','');
    set(gca,'yticklabel','');
    
    xlabel('Response seq','Position',[mean(x_lim) max(y_lim)+2],'FontSize',18,'FontWeight','bold');
    ylabel('Stimuli seq','Position',[min(x_lim)-2.5 mean(y_lim)],'FontSize',18,'FontWeight','bold');
    title('\fontsize{20}{choice - no choice}');
    
    colorbar('Ticks', [-1:0.1:1]);
    %colorbar;
end



%% Save data of paper fig2
if ifSave_paperFig2Data == 1
    currentName = ['fig2_data', '_'];
    fileName_figData = getFileNameMonkey_MAT(subject_name, currentName, monkey_name);
    
    fig2_panelA_data = struct;
    fig2_panelA_data.internalAccuracy_noChoice = internalNoChoiceAccuracy;
    fig2_panelA_data.internalAccuracy_choice = internalChoiceAccuracy;
    fig2_panelA_data.offloadingAccuracy_choice = offloadChoiceAccuracy;
    fig2_panelA_data.allAccuracy_choice = allChoiceAccuracy;
    
    fig2_panelB_data = struct;
    fig2_panelB_data.offloadingRate = ProbOffloadingInSeqLength;
    
    save(fileName_figData, 'fig2_panelA_data');
    save(fileName_figData, 'fig2_panelB_data','-append');
end


%% Plot and compare between gAcc_choice and gAcc_noChoice, resorted
% To resort sequences as boost effect
nonZero_index_collapsed = gAcc_choice_collapsed;
temp_subtract_collapsed = gAcc_choice_collapsed;
nonZero_index_sorted_collapsed= gAcc_choice_collapsed;

for target_seqLength=1:pointKindsNum
    if target_seqLength == 1
        temp_valid = valid_index(1:sum(numSeq(1)));
    else
        temp_valid = valid_index(sum(numSeq(1:target_seqLength-1))+1:sum(numSeq(1:target_seqLength)));
    end
    nonZero_index_collapsed{target_seqLength} = find(temp_valid==1);
    
    
    temp_subtract_collapsed{target_seqLength} = ...
        ( ( gAcc_choice_collapsed{target_seqLength}(nonZero_index_collapsed{target_seqLength}) ) ...
        - gAcc_noChoice_collapsed{target_seqLength}(nonZero_index_collapsed{target_seqLength}));
    
    [B,I] = sort(temp_subtract_collapsed{target_seqLength}, 'descend');
    nonZero_index_sorted_collapsed{target_seqLength} = nonZero_index_collapsed{target_seqLength}(I);
    for tempi=1:size(target_seqSet{target_seqLength}, 1)
        if ismember(tempi, nonZero_index_sorted_collapsed{target_seqLength}) == 0
            nonZero_index_sorted_collapsed{target_seqLength} = ...
                [nonZero_index_sorted_collapsed{target_seqLength}, tempi];
        end
    end
end
nonZero_index_sorted_collapsed;



%% Figure7
if ifFigure7 == 1
    fig7 = figure(7);
    set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    
    nonZero_index_sorted_collapsed;
    seqSet_inOne;
    gAcc_choice_collapsed;
    gAcc_noChoice_collapsed;
    
    for target_seqLength=1:pointKindsNum
        temp_sorted_index = nonZero_index_sorted_collapsed{target_seqLength};
        
        current_numSeq = length(target_seqSet{target_seqLength});
        
        subplot(pointKindsNum, 1, target_seqLength);
        plot(1:current_numSeq, gAcc_choice_collapsed{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1.5,'MarkerSize', 5, 'MarkerFaceColor', 'k', 'Color', [0.25 0.75 0.25]);
        hold on
        plot(1:current_numSeq, gAcc_noChoice_collapsed{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'Color', [0.25 0.25 0.25]);
        hold on
        
        if length(nonZero_index_collapsed{target_seqLength}) < numSeq(target_seqLength)
            plot(0.5+([1 1]*length(nonZero_index_collapsed{target_seqLength})),...
                [0 1], '-', 'LineWidth', 3, 'Color', [0.25 0.25 0.75]);
            hold on
        end
        
        temp_boost = gAcc_choice_collapsed{target_seqLength}(temp_sorted_index) - gAcc_noChoice_collapsed{target_seqLength}(temp_sorted_index);
        temp_zeroBound = find(temp_boost<0, 1);
        if isempty(temp_zeroBound) == 0
            plot(([1 1]*temp_zeroBound)-0.5,...
                [0 1], '-', 'LineWidth', 3, 'Color', [0.75 0.25 0.25]);
            hold on
        end
        
        if if_gAcc_errorBar == 1
            temp_valid_length = length(nonZero_index_collapsed{target_seqLength});
            % plot shaded erro bar for choice
            temp_SEM = gAcc_choice_merged_SEM{target_seqLength}(temp_sorted_index(1:temp_valid_length));
            y_CI = temp_SEM;
            x = 1:temp_valid_length;
            y = gAcc_choice_collapsed{target_seqLength}(temp_sorted_index(1:temp_valid_length))';
            x_conf = [x x(end:-1:1)] ;
            y_conf = [y+y_CI y(end:-1:1)-y_CI(end:-1:1)];
            h = fill(x_conf,y_conf,[0.25 0.75 0.25]);
            set(h,'edgealpha',0,'facealpha',0.2)
            %其中，'edgealpha'值为1时表示显示边界，为0时表示没有边界；
            %'facealpha'的值表示填充颜色的透明度，值为1时表示颜色最深。
            hold on
            
            % plot shaded erro bar for noChoice
            temp_SEM = gAcc_noChoice_merged_SEM{target_seqLength}(temp_sorted_index(1:temp_valid_length));
            y_CI = temp_SEM;
            x = 1:temp_valid_length;
            y = gAcc_noChoice_collapsed{target_seqLength}(temp_sorted_index(1:temp_valid_length))';
            x_conf = [x x(end:-1:1)] ;
            y_conf = [y+y_CI y(end:-1:1)-y_CI(end:-1:1)];
            h = fill(x_conf,y_conf,[0.25 0.25 0.25]);
            set(h,'edgealpha',0,'facealpha',0.2)
            %其中，'edgealpha'值为1时表示显示边界，为0时表示没有边界；
            %'facealpha'的值表示填充颜色的透明度，值为1时表示颜色最深。
            hold on
            
            temp_significance = if_signifi_p_choice_noChoice{target_seqLength}(temp_sorted_index(1:temp_valid_length));
            temp_significance_sorted_index = find(temp_significance==true);
            
            text(temp_significance_sorted_index,ones(1, length(temp_significance_sorted_index))-0.05,'*','Color','black','FontSize',18,'FontWeight','bold',...
                'HorizontalAlignment','center');
        end
        
        if target_seqLength == 1
            legend('choice accuracy',...
                'no choice accuracy','Location','southwest','fontsize',9)
        end
        
        ylim([0 1]);
        set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
        set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
        
        set(gca,'XTickLabel',seqSet_inOne{target_seqLength}(temp_sorted_index));%给坐标加标签
        ylabel(sprintf('Accuracy'), 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Choice and no choice accuracy, resorted, length=%d', target_seqLength), 'FontSize',15);
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
    
    currentFigName = ['gAcc_resorted', '_'];
    % to generate a unique file name for saving figure
    fileName_fig7 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(7) == 1
        exportgraphics(fig7, fileName_fig7, 'BackgroundColor', 'none');
    end
end


%% Plot contribution of delete one & add one, and interference effect
% Only deal with no choice condition!
if if_responseStructure_noChoice0_choice1 == 0
    current_responseMatrix = responseMatrix.noChoice_n11n;
elseif if_responseStructure_noChoice0_choice1 == 1
    current_responseMatrix = responseMatrix.choice_n11n;
end
neighbor_lengthError_indexMatrix;
neighbor_lengthError_withoutDiagonal_indexMatrix = neighbor_lengthError_indexMatrix;
diagonal_indexMatrix = zeros(size(current_responseMatrix,1), size(current_responseMatrix,2));
for tempi=1:sum(numSeq)
    neighbor_lengthError_withoutDiagonal_indexMatrix(tempi, tempi) = 0;
    diagonal_indexMatrix(tempi, tempi) = 1;
end

lengthCorrect_indexMatrix = ...
    false(size(current_responseMatrix,1), size(current_responseMatrix,2));
for tempi=1:pointKindsNum
    if tempi == 1
        temp_range = 1:numSeq(1);
    else
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    end
    temp_range;
    lengthCorrect_indexMatrix(temp_range, temp_range) = true;
end
lengthCorrect_indexMatrix;
lengthCorrect_withoutDiagonal_indexMatrix = lengthCorrect_indexMatrix;
for tempi=1:sum(numSeq)
    lengthCorrect_withoutDiagonal_indexMatrix(tempi, tempi) = 0;
end

interference_indexMatrix = zeros(size(current_responseMatrix,1), size(current_responseMatrix,2));
% Get interference_indexMatrix
for tempi=1:sum(numSeq)
    interference_indexMatrix(tempi, :);
    temp_seq = target_seqSet_inOne{tempi};
    temp_length = length(temp_seq);
    
    for tempj=1:numSeq(temp_length)
        %temp_seq2 = target_seqSet{temp_length}{tempj};
        if temp_length == 1
            tempi2 = tempj;
        else
            tempi2 = sum(numSeq(1:temp_length-1)) + tempj;
        end
        
        if tempi2 == tempi
            continue
        end
        
        temp_seq2 = target_seqSet_inOne{tempi2};
        temp_if_interference = 1;
        
        %         for tempk=1:temp_length
        %             if abs(temp_seq(tempk)-temp_seq2(tempk)) > 1
        %                 temp_if_interference = 0;
        %             end
        %         end
        if sum(abs(temp_seq-temp_seq2)) ~= 1
            temp_if_interference = 0;
        end
        
        if temp_if_interference == 1
            interference_indexMatrix(tempi,tempi2) = 1;
        end
    end
end

current_responseMatrix;
diagonal_indexMatrix;
neighbor_lengthError_withoutDiagonal_indexMatrix;
interference_indexMatrix;
valid_index;

condition_1_indexMatrix = diagonal_indexMatrix;
condition_2_indexMatrix = diagonal_indexMatrix | neighbor_lengthError_withoutDiagonal_indexMatrix;
condition_3_indexMatrix = diagonal_indexMatrix | interference_indexMatrix;
condition_full_indexMatrix = diagonal_indexMatrix | neighbor_lengthError_withoutDiagonal_indexMatrix ...
    | interference_indexMatrix;


condition_1_responseMatrix = current_responseMatrix .* condition_1_indexMatrix;
condition_2_responseMatrix = current_responseMatrix .* condition_2_indexMatrix;
condition_3_responseMatrix = current_responseMatrix .* condition_3_indexMatrix;
condition_full_responseMatrix = current_responseMatrix .* condition_full_indexMatrix;

condition_1_responseProportion = sum(condition_1_responseMatrix,2) ./ sum(current_responseMatrix,2);
condition_2_responseProportion = sum(condition_2_responseMatrix,2) ./ sum(current_responseMatrix,2);
condition_3_responseProportion = sum(condition_3_responseMatrix,2) ./ sum(current_responseMatrix,2);
condition_full_responseProportion = sum(condition_full_responseMatrix,2) ./ sum(current_responseMatrix,2);

condition_responseProportion = ...
    {condition_1_responseProportion';    % diagonal
    condition_2_responseProportion';     % diagonal + delete/add one
    condition_3_responseProportion';     % diagonal + interference
    condition_full_responseProportion'}';% full

condition_responseProportion;
valid_index;

lengthCorrect_withoutDiagonal_responseMatrix = current_responseMatrix .* lengthCorrect_withoutDiagonal_indexMatrix;
interference_responseMatrix = current_responseMatrix .* interference_indexMatrix;
interference_responseProportion_underSameLengthButError = ...
    sum(interference_responseMatrix,2) ./ sum(lengthCorrect_withoutDiagonal_responseMatrix,2);



%% Figure9
if ifFigure9 == 1
    fig9 = figure(9);
    set(gcf,'Position',[0 50 700 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    length_rangeHead = seq_length_rangeHead;
    length_rangeTail = seq_length_rangeTail;
    
    temp_color = [0.00 0.00 0.00;
        0.28 0.28 0.28;
        0.56 0.56 0.56;
        0.84 0.84 0.84];
    temp_color = temp_color(end:-1:1, :);
    
    
    temp_mean = zeros(pointKindsNum, 4);
    temp_SEM = zeros(pointKindsNum, 4);
    for tempi=1:4
        condition_responseProportion{tempi};
        for target_seqLength=1:pointKindsNum
            
            if target_seqLength == 1
                temp_range = 1:numSeq(1);
            else
                temp_range = sum(numSeq(1:target_seqLength-1))+1:sum(numSeq(1:target_seqLength));
            end
            
            temp_valid_index = valid_index(temp_range);
            temp_responseProportion = condition_responseProportion{tempi}(temp_range);
            temp_responseProportion_valid = temp_responseProportion(temp_valid_index);
            
            temp_mean(target_seqLength, tempi) = mean(temp_responseProportion_valid);
            temp_SEM(target_seqLength, tempi) = std(temp_responseProportion_valid) / sqrt(sum(temp_valid_index));
            
        end
        errorbar(length_rangeHead:length_rangeTail,temp_mean(:, tempi),temp_SEM(:, tempi),...
            '-o','Color', temp_color(tempi, :),'LineWidth',1,'CapSize',12,'MarkerSize',5);
        hold on
    end
    
    legend('diagonal','diagonal + delete/add one','diagonal + interference','full',...
        'Location','southeast','fontsize',11);
    
    
    ylim([0 1]);
    set(gca, 'FontSize', 20)
    % set(gca,'YTick',[-0.3:0.1:0.3]);%设置要显示坐标刻度的范围
    set(gca,'XLim',[length_rangeHead-1 length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[length_rangeHead:1:length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[length_rangeHead:1:length_rangeTail]);%给坐标加标签
    set(gca,'box','off');% 取消右、上边框
    xlabel('Length', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Probability', 'FontSize', 18, 'FontWeight', 'bold');
    if if_responseStructure_noChoice0_choice1 == 0
        title({'\fontsize{20}{\bf No choice response structure}'});
    elseif if_responseStructure_noChoice0_choice1 == 1
        title({'\fontsize{20}{\bf Choice response structure}'});
    end
    
end

%% Rule out (delete/add one) or (others) trials and normalize responseMatrix, then calculate boost again
if_ruleOut_deleteAddOne0_others1;

responseMatrix.noChoice_n11n;
responseMatrix.choice_n11n;
new_responseMatrix = struct;

current_responseMatrix = responseMatrix.choice_n11n;
if if_ruleOut_deleteAddOne0_others1 == 0
    tempNew_responseMatrix = current_responseMatrix .* (neighbor_lengthError_withoutDiagonal_indexMatrix==0);
elseif if_ruleOut_deleteAddOne0_others1 == 1
    otherError_indexMatrix = (current_responseMatrix>0) & ...
        (neighbor_lengthError_withoutDiagonal_indexMatrix==0) & ...
        (interference_indexMatrix==0) & ...
        (diagonal_indexMatrix==0);
    tempNew_responseMatrix = current_responseMatrix .* (otherError_indexMatrix==0);
end
new_responseMatrix_n11n = tempNew_responseMatrix ./ sum(tempNew_responseMatrix,2);
new_responseMatrix.choice_n11n = new_responseMatrix_n11n;

current_responseMatrix = responseMatrix.noChoice_n11n;
if if_ruleOut_deleteAddOne0_others1 == 0
    tempNew_responseMatrix = current_responseMatrix .* (neighbor_lengthError_withoutDiagonal_indexMatrix==0);
elseif if_ruleOut_deleteAddOne0_others1 == 1
    otherError_indexMatrix = (current_responseMatrix>0) & ...
        (neighbor_lengthError_withoutDiagonal_indexMatrix==0) & ...
        (interference_indexMatrix==0) & ...
        (diagonal_indexMatrix==0);
    tempNew_responseMatrix = current_responseMatrix .* (otherError_indexMatrix==0);
end
new_responseMatrix_n11n = tempNew_responseMatrix ./ sum(tempNew_responseMatrix,2);
new_responseMatrix.noChoice_n11n = new_responseMatrix_n11n;

new_gAcc_choice_collapsed = cell(1, pointKindsNum);
new_gAcc_noChoice_collapsed = cell(1, pointKindsNum);
for tempi=1:pointKindsNum
    new_gAcc_choice_collapsed{tempi} = zeros(numSeq(tempi), 1);
    new_gAcc_noChoice_collapsed{tempi} = zeros(numSeq(tempi), 1);
    for tempj=1:numSeq(tempi)
        if tempi==1
            tempTotal = tempj;
        else
            tempTotal = sum(numSeq(1:tempi-1)) + tempj;
        end
        
        new_gAcc_choice_collapsed{tempi}(tempj) = new_responseMatrix.choice_n11n(tempTotal, tempTotal);
        new_gAcc_noChoice_collapsed{tempi}(tempj) = new_responseMatrix.noChoice_n11n(tempTotal, tempTotal);
        if tempi==4
            a=1;
        end
    end
end
new_gAcc_choice_collapsed;
new_gAcc_noChoice_collapsed;

nonZero_index_collapsed;
new_temp_subtract_collapsed = gAcc_choice_collapsed;
new_nonZero_index_sorted_collapsed= gAcc_choice_collapsed;

for target_seqLength=1:pointKindsNum
    
    new_temp_subtract_collapsed{target_seqLength} = ...
        ( ( new_gAcc_choice_collapsed{target_seqLength}(nonZero_index_collapsed{target_seqLength}) ) ...
        - new_gAcc_noChoice_collapsed{target_seqLength}(nonZero_index_collapsed{target_seqLength}));
    
    
    [B,I] = sort(new_temp_subtract_collapsed{target_seqLength},'descend','MissingPlacement','last');
    new_nonZero_index_sorted_collapsed{target_seqLength} = nonZero_index_collapsed{target_seqLength}(I);
    for tempi=1:size(target_seqSet{target_seqLength}, 1)
        if ismember(tempi, new_nonZero_index_sorted_collapsed{target_seqLength}) == 0
            new_nonZero_index_sorted_collapsed{target_seqLength} = ...
                [new_nonZero_index_sorted_collapsed{target_seqLength}, tempi];
        end
    end
end
new_nonZero_index_sorted_collapsed;

%% Figure10
if ifFigure10 == 1
    fig10 = figure(10);
    set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    
    new_nonZero_index_sorted_collapsed;
    seqSet_inOne;
    new_gAcc_choice_collapsed;
    new_gAcc_noChoice_collapsed;
    
    for target_seqLength=1:pointKindsNum
        temp_sorted_index = new_nonZero_index_sorted_collapsed{target_seqLength};
        
        current_numSeq = length(target_seqSet{target_seqLength});
        
        subplot(pointKindsNum, 1, target_seqLength);
        plot(1:current_numSeq, new_gAcc_choice_collapsed{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 2,'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.25 0.75 0.25]);
        hold on
        plot(1:current_numSeq, new_gAcc_noChoice_collapsed{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Color', [0.25 0.25 0.25]);
        hold on
        if length(nonZero_index_collapsed{target_seqLength}) < numSeq(target_seqLength)
            plot(0.5+([1 1]*length(nonZero_index_collapsed{target_seqLength})),...
                [0 1], '-', 'LineWidth', 3, 'Color', [0.25 0.25 0.75]);
            hold on
        end
        
        if target_seqLength == 1
            legend('choice accuracy',...
                'no choice accuracy','Location','southwest','fontsize',9)
        end
        
        ylim([0 1]);
        set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
        set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
        
        set(gca,'XTickLabel',seqSet_inOne{target_seqLength}(temp_sorted_index));%给坐标加标签
        ylabel(sprintf('Accuracy'), 'FontSize', 12, 'FontWeight', 'bold');
        if if_ruleOut_deleteAddOne0_others1 == 0
            title(sprintf('Choice and no choice accuracy, resorted, rule out delete/add one, length=%d',target_seqLength), 'FontSize',15);
        elseif if_ruleOut_deleteAddOne0_others1 == 1
            title(sprintf('Choice and no choice accuracy, resorted, rule out others, length=%d',target_seqLength), 'FontSize',15);
        end
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
    
    if if_ruleOut_deleteAddOne0_others1 == 0
        currentFigName = ['gAcc_resorted_ruleOut_deleteAddOne', '_'];
    elseif if_ruleOut_deleteAddOne0_others1 == 1
        currentFigName = ['gAcc_resorted_ruleOut_others', '_'];
    end
    % to generate a unique file name for saving figure
    fileName_fig10 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(10) == 1
        exportgraphics(fig10, fileName_fig10, 'BackgroundColor', 'none');
    end
end

%% Compute proportion error types of  delete/add one, interference, others in each seq.
responseMatrix.choice_n11n;
responseMatrix.noChoice_n11n;

errorType_indexMatrix = struct;
errorType_indexMatrix.deleteAddOne = neighbor_lengthError_withoutDiagonal_indexMatrix==1;
errorType_indexMatrix.interference = interference_indexMatrix==1;
errorType_indexMatrix.others = struct;

current_responseMatrix = responseMatrix.choice_n11n;
otherError_indexMatrix = (current_responseMatrix>0) & ...
    (neighbor_lengthError_withoutDiagonal_indexMatrix==0) & ...
    (interference_indexMatrix==0) & ...
    (diagonal_indexMatrix==0);
errorType_indexMatrix.others.choice = otherError_indexMatrix==1;

current_responseMatrix = responseMatrix.noChoice_n11n;
otherError_indexMatrix = (current_responseMatrix>0) & ...
    (neighbor_lengthError_withoutDiagonal_indexMatrix==0) & ...
    (interference_indexMatrix==0) & ...
    (diagonal_indexMatrix==0);
errorType_indexMatrix.others.noChoice = otherError_indexMatrix==1;

errorType_indexMatrix;
errorType_seqSummary = struct;

errorType_seqSummary.choice = struct;
current_responseMatrix = responseMatrix.choice_n11n;
errorType_seqSummary.choice.deleteAddOne = sum(current_responseMatrix .* errorType_indexMatrix.deleteAddOne, 2)';
errorType_seqSummary.choice.interference = sum(current_responseMatrix .* errorType_indexMatrix.interference, 2)';
errorType_seqSummary.choice.others = sum(current_responseMatrix .* errorType_indexMatrix.others.choice, 2)';

errorType_seqSummary.noChoice = struct;
current_responseMatrix = responseMatrix.noChoice_n11n;
errorType_seqSummary.noChoice.deleteAddOne = sum(current_responseMatrix .* errorType_indexMatrix.deleteAddOne, 2)';
errorType_seqSummary.noChoice.interference = sum(current_responseMatrix .* errorType_indexMatrix.interference, 2)';
errorType_seqSummary.noChoice.others = sum(current_responseMatrix .* errorType_indexMatrix.others.noChoice, 2)';

errorType_seqSummary.delta = struct;
errorType_seqSummary.delta.deleteAddOne = errorType_seqSummary.choice.deleteAddOne - errorType_seqSummary.noChoice.deleteAddOne;
errorType_seqSummary.delta.interference = errorType_seqSummary.choice.interference - errorType_seqSummary.noChoice.interference;
errorType_seqSummary.delta.others = errorType_seqSummary.choice.others - errorType_seqSummary.noChoice.others;

errorType_seqSummary_collapsed = struct;
errorType_seqSummary_collapsed.choice = struct;
errorType_seqSummary_collapsed.choice.deleteAddOne = cell(1, pointKindsNum);
errorType_seqSummary_collapsed.choice.interference = cell(1, pointKindsNum);
errorType_seqSummary_collapsed.choice.others = cell(1, pointKindsNum);
errorType_seqSummary_collapsed.noChoice = struct;
errorType_seqSummary_collapsed.noChoice.deleteAddOne = cell(1, pointKindsNum);
errorType_seqSummary_collapsed.noChoice.interference = cell(1, pointKindsNum);
errorType_seqSummary_collapsed.noChoice.others = cell(1, pointKindsNum);
errorType_seqSummary_collapsed.delta = struct;
errorType_seqSummary_collapsed.delta.deleteAddOne = cell(1, pointKindsNum);
errorType_seqSummary_collapsed.delta.interference = cell(1, pointKindsNum);
errorType_seqSummary_collapsed.delta.others = cell(1, pointKindsNum);

for tempi=1:pointKindsNum
    if tempi == 1
        temp_range = 1:numSeq(1);
    else
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    end
    errorType_seqSummary_collapsed.choice.deleteAddOne{tempi} = errorType_seqSummary.choice.deleteAddOne(temp_range);
    errorType_seqSummary_collapsed.choice.interference{tempi} = errorType_seqSummary.choice.interference(temp_range);
    errorType_seqSummary_collapsed.choice.others{tempi} = errorType_seqSummary.choice.others(temp_range);
    
    errorType_seqSummary_collapsed.noChoice.deleteAddOne{tempi} = errorType_seqSummary.noChoice.deleteAddOne(temp_range);
    errorType_seqSummary_collapsed.noChoice.interference{tempi} = errorType_seqSummary.noChoice.interference(temp_range);
    errorType_seqSummary_collapsed.noChoice.others{tempi} = errorType_seqSummary.noChoice.others(temp_range);
    
    errorType_seqSummary_collapsed.delta.deleteAddOne{tempi} = errorType_seqSummary.delta.deleteAddOne(temp_range);
    errorType_seqSummary_collapsed.delta.interference{tempi} = errorType_seqSummary.delta.interference(temp_range);
    errorType_seqSummary_collapsed.delta.others{tempi} = errorType_seqSummary.delta.others(temp_range);
end

%% Figure 11, choice error type
fig11 = figure(11);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

for target_seqLength=1:pointKindsNum
    temp_sorted_index = nonZero_index_sorted_collapsed{target_seqLength};
    current_numSeq = length(target_seqSet{target_seqLength});
    
    subplot(pointKindsNum, 1, target_seqLength);
    plot(1:current_numSeq, errorType_seqSummary_collapsed.choice.deleteAddOne{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1.5,'MarkerSize', 4, 'MarkerFaceColor', 'k', 'Color', '#0072BD');
    hold on
    plot(1:current_numSeq, errorType_seqSummary_collapsed.choice.interference{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'Color', '#D95319');
    hold on
    plot(1:current_numSeq, errorType_seqSummary_collapsed.choice.others{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'Color', '#7E2F8E');
    hold on
    if length(nonZero_index_collapsed{target_seqLength}) < numSeq(target_seqLength)
        plot(0.5+([1 1]*length(nonZero_index_collapsed{target_seqLength})),...
            [0 1], '-', 'LineWidth', 3, 'Color', [0.25 0.25 0.75]);
        hold on
    end
    
    temp_boost = gAcc_choice_collapsed{target_seqLength}(temp_sorted_index) - gAcc_noChoice_collapsed{target_seqLength}(temp_sorted_index);
    temp_zeroBound = find(temp_boost<0, 1);
    if isempty(temp_zeroBound) == 0
        plot(([1 1]*temp_zeroBound)-0.5,...
            [0 1], '-', 'LineWidth', 3, 'Color', [0.75 0.25 0.25]);
        hold on
    end
    
    %if target_seqLength == 1
    if target_seqLength == 1 || target_seqLength == 2
        legend('delete/add One','interference','others','Location','southwest','fontsize',9)
    end
    
    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
    
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength}(temp_sorted_index));%给坐标加标签
    ylabel(sprintf('Proportion'), 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Choice, error types in each sequence, length=%d',target_seqLength), 'FontSize',15);
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

currentFigName = ['choice_errorType', '_'];
% to generate a unique file name for saving figure
fileName_fig11 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(11) == 1
    exportgraphics(fig11, fileName_fig11, 'BackgroundColor', 'none');
end

%% Figure 12, noChoice error type
fig12 = figure(12);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

for target_seqLength=1:pointKindsNum
    temp_sorted_index = nonZero_index_sorted_collapsed{target_seqLength};
    current_numSeq = length(target_seqSet{target_seqLength});
    
    subplot(pointKindsNum, 1, target_seqLength);
    plot(1:current_numSeq, errorType_seqSummary_collapsed.noChoice.deleteAddOne{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1.5,'MarkerSize', 4, 'MarkerFaceColor', 'k', 'Color', '#0072BD');
    hold on
    plot(1:current_numSeq, errorType_seqSummary_collapsed.noChoice.interference{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'Color', '#D95319');
    hold on
    plot(1:current_numSeq, errorType_seqSummary_collapsed.noChoice.others{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'Color', '#7E2F8E');
    hold on
    if length(nonZero_index_collapsed{target_seqLength}) < numSeq(target_seqLength)
        plot(0.5+([1 1]*length(nonZero_index_collapsed{target_seqLength})),...
            [0 1], '-', 'LineWidth', 3, 'Color', [0.25 0.25 0.75]);
        hold on
    end
    
    temp_boost = gAcc_choice_collapsed{target_seqLength}(temp_sorted_index) - gAcc_noChoice_collapsed{target_seqLength}(temp_sorted_index);
    temp_zeroBound = find(temp_boost<0, 1);
    if isempty(temp_zeroBound) == 0
        plot(([1 1]*temp_zeroBound)-0.5,...
            [0 1], '-', 'LineWidth', 3, 'Color', [0.75 0.25 0.25]);
        hold on
    end
    
    
    if target_seqLength == 1
        legend('delete/add One','interference','others','Location','southwest','fontsize',9)
    end
    
    ylim([0 1]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
    
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength}(temp_sorted_index));%给坐标加标签
    ylabel(sprintf('Proportion'), 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('No choice, error types in each sequence, length=%d',target_seqLength), 'FontSize',15);
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

currentFigName = ['noChoice_errorType', '_'];
% to generate a unique file name for saving figure
fileName_fig12 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(12) == 1
    exportgraphics(fig12, fileName_fig12, 'BackgroundColor', 'none');
end

plot_yMax = max([errorType_seqSummary.delta.deleteAddOne,...
    errorType_seqSummary.delta.interference,errorType_seqSummary.delta.others]);
plot_yMin = min([errorType_seqSummary.delta.deleteAddOne,...
    errorType_seqSummary.delta.interference,errorType_seqSummary.delta.others]);

%% Figure 13, choice - noChoice error type
fig13 = figure(13);
set(gcf,'Position',[0 50 1280 720]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

for target_seqLength=1:pointKindsNum
    temp_sorted_index = nonZero_index_sorted_collapsed{target_seqLength};
    current_numSeq = length(target_seqSet{target_seqLength});
    
    subplot(pointKindsNum, 1, target_seqLength);
    plot(1:current_numSeq, errorType_seqSummary_collapsed.delta.deleteAddOne{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1,'MarkerSize', 3, 'MarkerFaceColor', 'k', 'Color', '#0072BD');
    hold on
    plot(1:current_numSeq, errorType_seqSummary_collapsed.delta.interference{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor', 'k', 'Color', '#D95319');
    hold on
    plot(1:current_numSeq, errorType_seqSummary_collapsed.delta.others{target_seqLength}(temp_sorted_index),'-s', 'LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor', 'k', 'Color', '#7E2F8E');
    hold on
    if length(nonZero_index_collapsed{target_seqLength}) < numSeq(target_seqLength)
        plot(0.5+([1 1]*length(nonZero_index_collapsed{target_seqLength})),...
            [plot_yMin plot_yMax], '-', 'LineWidth', 3, 'Color', [0.25 0.25 0.75]);
        hold on
    end
    
    temp_boost = gAcc_choice_collapsed{target_seqLength}(temp_sorted_index) - gAcc_noChoice_collapsed{target_seqLength}(temp_sorted_index);
    temp_zeroBound = find(temp_boost<0, 1);
    if isempty(temp_zeroBound) == 0
        plot(([1 1]*temp_zeroBound)-0.5,...
            [plot_yMin plot_yMax], '-', 'LineWidth', 3, 'Color', [0.75 0.25 0.25]);
        hold on
    end
    
    
    plot(1:current_numSeq,zeros(1, current_numSeq),'--','Color','k');
    hold on
    
    %if target_seqLength == 1
    if target_seqLength == 1 || target_seqLength == 2
        legend('delete/add One','interference','others','Location','southwest','fontsize',8)
    end
    
    %ylim([0 1]);
    ylim([plot_yMin plot_yMax]);
    set(gca,'XLim',[0 current_numSeq+1]);%X轴的数据显示范围
    set(gca,'XTick',[1:1:current_numSeq]);%设置要显示坐标刻度的范围
    
    set(gca,'XTickLabel',seqSet_inOne{target_seqLength}(temp_sorted_index));%给坐标加标签
    ylabel(sprintf('Proportion'), 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Choice - no choice, error types in each sequence, length=%d',target_seqLength), 'FontSize',15);
    % 获取xticklabel的忿
    xtl=get(gca,'XTickLabel');
    % 获取xtick的忿
    xt=get(gca,'XTick');
    % 获取ytick的忿
    yt=get(gca,'YTick');
    % 设置text的x坐标位置仿
    xtextp=xt;
    % 设置text的y坐标位置仿
    %ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.08;
    if numFrames < 8
        ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.15;
    elseif numFrames == 8
        ytextp=(yt(1)-0.2*(yt(2)-yt(1)))*ones(1,length(xt))-0.5;
    end
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

currentFigName = ['boost_errorType', '_'];
% to generate a unique file name for saving figure
fileName_fig13 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(13) == 1
    exportgraphics(fig13, fileName_fig13, 'BackgroundColor', 'none');
end

%% Calculate and plot average boost in length that more than 0.
% if_newBoost = 1;

nonZero_index_collapsed;
valid_index;
choiceMinusNoChoice_inOne;

if if_newBoost == 0
    current_choiceMinusNoChoice_inOne = choiceMinusNoChoice_inOne';
elseif if_newBoost == 1
    new_gAcc_choice_inOne = [];
    new_gAcc_noChoice_inOne = [];
    for tempi=1:pointKindsNum
        new_gAcc_choice_inOne = [new_gAcc_choice_inOne new_gAcc_choice_collapsed{tempi}'];
        new_gAcc_noChoice_inOne = [new_gAcc_noChoice_inOne new_gAcc_noChoice_collapsed{tempi}'];
    end
    current_choiceMinusNoChoice_inOne = new_gAcc_choice_inOne - new_gAcc_noChoice_inOne;
end

% valid_boostMoreThan0_index = (current_choiceMinusNoChoice_inOne>0)' & valid_index;
valid_boostMoreThan0_index = (current_choiceMinusNoChoice_inOne>0) & valid_index;
choiceMinusNoChoice_moreThan0_valid = cell(1, pointKindsNum);
temp_mean = zeros(1, pointKindsNum);
temp_SEM = zeros(1, pointKindsNum);

%% Figure 14, choiceMinusNoChoice, more than 0
if ifFigure14 == 1
    fig14 = figure(14);
    % set(gcf,'Position',[0 50 800 600]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    for target_seqLength=1:pointKindsNum
        if target_seqLength == 1
            temp_range = 1:numSeq(1);
        else
            temp_range = sum(numSeq(1:target_seqLength-1))+1:sum(numSeq(1:target_seqLength));
        end
        temp_valid = valid_boostMoreThan0_index(temp_range);
        temp_choiceMinusNoChoice = current_choiceMinusNoChoice_inOne(temp_range)';
        
        temp_mean(target_seqLength) = mean(temp_choiceMinusNoChoice(temp_valid));
        temp_SEM(target_seqLength) = ...
            std(temp_choiceMinusNoChoice(temp_valid))./sqrt(sum(temp_valid));
        
        choiceMinusNoChoice_moreThan0_valid{target_seqLength} = temp_choiceMinusNoChoice(temp_valid);
    end
    errorbar(seq_length_rangeHead: seq_length_rangeTail, temp_mean,temp_SEM, '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'CapSize', 12, 'MarkerSize', 5);
    
    hold on
    plot([seq_length_rangeHead-1 seq_length_rangeHead: seq_length_rangeTail seq_length_rangeTail+1]...
        , zeros(1, pointKindsNum+2), '--');
    
    choiceMinusNoChoice_moreThan0_valid;
    for j=1:pointKindsNum
        tempZero = zeros(length(choiceMinusNoChoice_moreThan0_valid{j}), 1);
        tempOne = ones(length(choiceMinusNoChoice_moreThan0_valid{j}), 1);
        %paired-sample t-test
        [h_t(j),p_t(j)]=ttest(choiceMinusNoChoice_moreThan0_valid{j}, tempZero);
        %[h_t(j),p_t(j)]=ttest(choiceMinusNoChoice_moreThan0_valid{j}, tempOne);
        % Wilcoxon rank sum test
        p_rs(j) = ...
            ranksum(choiceMinusNoChoice_moreThan0_valid{j}, tempZero); %#ok<*SAGROW>
        
        p(j) = p_t(j);
        %p(j) = p_rs(j);
        
        %Bonferroni Correction
        
        tempTxt = sprintf('');
        if p(j) < 0.001/pointKindsNum
            tempTxt = sprintf('***');
        elseif p(j) < 0.01/pointKindsNum
            tempTxt = sprintf('**');
        elseif p(j) < 0.05/pointKindsNum
            tempTxt = sprintf('*');
        end
        if p(j) < 0.05/pointKindsNum
            text(j+seq_length_rangeHead-1,mean(choiceMinusNoChoice_moreThan0_valid{j})-0.08,tempTxt,'Color','black','FontSize',30,'FontWeight','bold',...
                'HorizontalAlignment','center');
        end
    end
    
    % regression by 'fitlm'
    seq_length_inOne = [];
    for tempi=1:3
        seq_length_inOne = [seq_length_inOne tempi*ones(1, numSeq(tempi))];
    end
    x = seq_length_inOne(valid_boostMoreThan0_index( 1:sum(numSeq(1:3)) ));
    y = current_choiceMinusNoChoice_inOne(valid_boostMoreThan0_index( 1:sum(numSeq(1:3)) ));
    mdl_boost = fitlm(x, y);
    beta_boost = mdl_boost.Coefficients.Estimate;
    r2_boost = mdl_boost.Rsquared.Adjusted;
    p_boost = coefTest(mdl_boost);
    
    
    
    ylim([-0.1 0.4]);
    % set(gca,'color',[231 231 230]/255);
    set(gca, 'FontSize', 20)
    set(gca,'YTick',[-0.1:0.1:0.4]);%设置要显示坐标刻度的范围
    set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
    set(gca,'box','off');% 取消右、上边框
    xlabel('Sequence length', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Accuracy', 'FontSize', 18, 'FontWeight', 'bold');
    if if_newBoost == 0
        title({'\fontsize{20}{\color[rgb]{0.25,0.75,0.25}internal|Choice} - internal|noChoice';
            'only more than 0'});
    elseif if_newBoost == 1
        if if_ruleOut_deleteAddOne0_others1 == 0
            title({'\fontsize{20}{\color[rgb]{0.25,0.75,0.25}internal|Choice} - internal|noChoice';
                'only more than 0, rule out delete/add one'});
        elseif if_ruleOut_deleteAddOne0_others1 == 1
            title({'\fontsize{20}{\color[rgb]{0.25,0.75,0.25}internal|Choice} - internal|noChoice';
                'only more than 0, rule out others'});
        end
    end
    
end

%% Compute and plot boost error types in length level
% temp_deleteAddOne = errorType_seqSummary_collapsed.delta.deleteAddOne;
% temp_interference = errorType_seqSummary_collapsed.delta.interference;
% temp_others = errorType_seqSummary_collapsed.delta.others;
valid_index;
valid_boostLessThan0_index = (choiceMinusNoChoice_inOne'<0) & valid_index;
% valid_boostLessThan0_index = (choiceMinusNoChoice_inOne'<-0.03) & valid_index;
% valid_boostLessThan0_index = (choiceMinusNoChoice_inOne'<-0.01) & valid_index;

length_rangeHead = seq_length_rangeHead;
length_rangeTail = seq_length_rangeTail;


temp_mean = zeros(pointKindsNum, 3);
temp_SEM = zeros(pointKindsNum, 3);

for target_seqLength=1:pointKindsNum
    if target_seqLength == 1
        temp_range = 1:numSeq(1);
    else
        temp_range = sum(numSeq(1:target_seqLength-1))+1:sum(numSeq(1:target_seqLength));
    end
    temp_valid = valid_boostLessThan0_index(temp_range);
    
    temp_deleteAddOne = errorType_seqSummary_collapsed.delta.deleteAddOne{target_seqLength}(temp_valid);
    temp_interference = errorType_seqSummary_collapsed.delta.interference{target_seqLength}(temp_valid);
    temp_others = errorType_seqSummary_collapsed.delta.others{target_seqLength}(temp_valid);
    
    temp_mean(target_seqLength, 1) = mean(temp_deleteAddOne);
    temp_mean(target_seqLength, 2) = mean(temp_interference);
    temp_mean(target_seqLength, 3) = mean(temp_others);
    
    if sum(temp_valid) > 0
        temp_SEM(target_seqLength, 1) = std(temp_deleteAddOne) / sqrt(sum(temp_valid));
        temp_SEM(target_seqLength, 2) = std(temp_interference) / sqrt(sum(temp_valid));
        temp_SEM(target_seqLength, 3) = std(temp_others) / sqrt(sum(temp_valid));
    end
    
end
length_valid_range = ~isnan(temp_mean(:,1))';
temp_length = length_rangeHead:length_rangeTail;

% temp_color = [0.00 0.00 0.00;
%     0.28 0.28 0.28;
%     0.56 0.56 0.56;
%     0.84 0.84 0.84];
% temp_color = temp_color(end:-1:1, :);

temp_color = {'#0072BD';'#D95319';'#7E2F8E'};

%% Figure 15
if ifFigure15 == 1
    fig15 = figure(15);
    % set(gcf,'Position',[0 50 700 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    for tempi=1:3
        errorbar(temp_length(length_valid_range),temp_mean(length_valid_range, tempi),temp_SEM(length_valid_range, tempi),...
            '-o','Color', temp_color{tempi},'LineWidth',1,'CapSize',12,'MarkerSize',5);
        hold on
    end
    
    plot([seq_length_rangeHead-1 seq_length_rangeHead: seq_length_rangeTail seq_length_rangeTail+1]...
        , zeros(1, pointKindsNum+2), '--');
    
    legend('delete/add One','interference','others','Location','southwest','fontsize',11)
    
    ylim([-0.4 0.4]);
    set(gca, 'FontSize', 20)
    set(gca,'YTick',[-0.4:0.2:0.4]);%设置要显示坐标刻度的范围
    set(gca,'XLim',[length_rangeHead-1 length_rangeTail+1]);%X轴的数据显示范围
    set(gca,'XTick',[length_rangeHead:1:length_rangeTail]);%设置要显示坐标刻度的范围
    set(gca,'XTickLabel',[length_rangeHead:1:length_rangeTail]);%给坐标加标签
    set(gca,'box','off');% 取消右、上边框
    xlabel('Length', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Proportion', 'FontSize', 18, 'FontWeight', 'bold');
    title({'\fontsize{20}{\bf Boost error types in length}'});
    
end

a=seqSet_inOne_inOne(valid_boostLessThan0_index)';



%% Calculate and plot heatmap of boost (gAcc|choice - gAcc|noChoice) in each length and in every location
nonZero_index_collapsed;
target_seqSet;
gAcc_choice_collapsed;
gAcc_noChoice_collapsed;
choiceMinusNoChoice;

% if_fig8_filter = 1;
% if_fig8_boost_negative0_positive1_full2 = 1;

temp_index_collapsed = cell(1, pointKindsNum);

% if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 = 0;% only work in negative boost.

if if_fig8_boost_negative0_positive1_full2 == 2
    if if_fig8_filter == 1
        temp_index_collapsed = nonZero_index_collapsed;
        
    elseif if_fig8_filter == 0
        for tempi=1:pointKindsNum
            temp_index_collapsed{tempi} = 1:numSeq(tempi);
        end
    end
elseif if_fig8_boost_negative0_positive1_full2 == 0
    if if_fig8_filter == 1
        current_valid = (choiceMinusNoChoice_inOne'<0) & valid_index;
    elseif if_fig8_filter == 0
        current_valid = (choiceMinusNoChoice_inOne'<0);
    end
    for tempi=1:pointKindsNum
        if tempi == 1
            temp_range = 1:numSeq(1);
        else
            temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
        end
        temp_valid = current_valid(temp_range);
        temp_index_collapsed{tempi} = find(temp_valid==1);
    end
elseif if_fig8_boost_negative0_positive1_full2 == 1
    if if_fig8_filter == 1
        current_valid = (choiceMinusNoChoice_inOne'>0) & valid_index;
    elseif if_fig8_filter == 0
        current_valid = (choiceMinusNoChoice_inOne'>0);
    end
    for tempi=1:pointKindsNum
        if tempi == 1
            temp_range = 1:numSeq(1);
        else
            temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
        end
        temp_valid = current_valid(temp_range);
        temp_index_collapsed{tempi} = find(temp_valid==1);
    end
end

% for if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2
choiceMinusNoChoice_2 = choiceMinusNoChoice;
if if_fig8_boost_negative0_positive1_full2 == 0
    if  if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 0
        gAcc_choice_collapsed;
        errorType_seqSummary_collapsed.choice.deleteAddOne;
        for tempi=1:pointKindsNum
            choiceMinusNoChoice_2{tempi} = ...
                gAcc_choice_collapsed{tempi}./(gAcc_choice_collapsed{tempi}+errorType_seqSummary_collapsed.choice.deleteAddOne{tempi}')...
                - ...
                gAcc_noChoice_collapsed{tempi}./(gAcc_noChoice_collapsed{tempi}+errorType_seqSummary_collapsed.noChoice.deleteAddOne{tempi}');
            choiceMinusNoChoice_2{tempi}(isnan(choiceMinusNoChoice_2{tempi})) = 0;
        end
        
    elseif if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 1
        gAcc_choice_collapsed;
        errorType_seqSummary_collapsed.choice.interference;
        for tempi=1:pointKindsNum
            choiceMinusNoChoice_2{tempi} = ...
                gAcc_choice_collapsed{tempi}./(gAcc_choice_collapsed{tempi}+errorType_seqSummary_collapsed.choice.interference{tempi}')...
                - ...
                gAcc_noChoice_collapsed{tempi}./(gAcc_noChoice_collapsed{tempi}+errorType_seqSummary_collapsed.noChoice.interference{tempi}');
            choiceMinusNoChoice_2{tempi}(isnan(choiceMinusNoChoice_2{tempi})) = 0;
        end
    end
    
end
choiceMinusNoChoice_2;

boost_weightedRaw = cell(1, pointKindsNum);
temp_count = cell(1, pointKindsNum);

for target_seqLength=1:pointKindsNum
    boost_weightedRaw{target_seqLength} = zeros(1, numFrames);
    temp_count{target_seqLength} = zeros(1, numFrames);
    temp_index_collapsed{target_seqLength};
    for tempi=1:length(temp_index_collapsed{target_seqLength})
        temp_seqIndex = temp_index_collapsed{target_seqLength}(tempi);
        
        temp_seq = target_seqSet{target_seqLength}{temp_seqIndex};
        temp_length = length(temp_seq);
        %choiceMinusNoChoice_2{target_seqLength}(temp_seqIndex);
        for tempj=1:temp_length
            boost_weightedRaw{target_seqLength}(temp_seq(tempj)) = ...
                boost_weightedRaw{target_seqLength}(temp_seq(tempj)) + ...
                choiceMinusNoChoice_2{target_seqLength}(temp_seqIndex);
            temp_count{target_seqLength}(temp_seq(tempj)) = ...
                temp_count{target_seqLength}(temp_seq(tempj)) + 1;
        end
    end
end
boost_weighted = boost_weightedRaw;
for target_seqLength=1:pointKindsNum
    for tempi=1:numFrames
        boost_weighted{target_seqLength}(tempi) = ...
            boost_weightedRaw{target_seqLength}(tempi) / temp_count{target_seqLength}(tempi);
        
    end
end

boost_weighted;
boost_weighted_inOne = [];
for target_seqLength=1:pointKindsNum
    boost_weighted_inOne = [boost_weighted_inOne boost_weighted{target_seqLength}];
end


% Get some size for heatmap painting
showRadius = max(temp_load.basic_para.showBaseRect)/2;

radius = temp_load.basic_para.radius;
InitAngle_arc = temp_load.basic_para.InitAngle_arc;
screenXpixels = 1920;
screenYpixels = 1080;
screenXpixels_bias = temp_load.basic_para.screenXpixels_bias;
screenYpixels_bias = temp_load.basic_para.screenYpixels_bias;
% Screen positions of our  rectangles
squareXpos = ones(1, numFrames);
squareYpos = ones(1, numFrames);
%Arc-like points metrix
theta = zeros(1, numFrames);
for i = 1:numFrames
    theta(i) = (i-1)/(numFrames-1) * (180-2*InitAngle_arc)/180 * pi + InitAngle_arc/180 * pi;
    squareXpos(i) = -1 * radius * cos(theta(i));
    squareYpos(i) = tan(theta(i)) * squareXpos(i) + screenYpixels* 0.7;
    squareXpos(i) = squareXpos(i) + screenXpixels/2;
end
squareXpos = squareXpos + screenXpixels_bias;
squareYpos = squareYpos + screenYpixels_bias;
crossCenter = [screenXpixels/2+screenXpixels_bias screenYpixels/2+screenYpixels_bias+205-20-25+30];
layout_x = [squareXpos crossCenter(1)];
layout_y = [squareYpos crossCenter(2)];


boost_weighted;
%% Figure8
if ifFigure8 == 1
    fig8 = figure(8);
    set(gcf,'Position',[0 50 800 600]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
    
    if if_fig8_boost_negative0_positive1_full2 == 2
        my_color = cool;
        backgroundColor = [1 1 1];
    else
        %my_color = gray;
        %my_color = my_color(end:-1:1,:);
        %backgroundColor = [0.3010 0.7450 0.9330]*0.5;
        
        %my_color = parula;
        my_color = gray;
        my_color = my_color(end:-1:1,:);
        %my_color = my_color(1:200, :);
        backgroundColor = [1 1 1];
    end
    
    %my_color = jet;
    colormap(my_color);
    
    
    temp_max = max(boost_weighted_inOne);
    temp_min = min(boost_weighted_inOne);
    
    boost_weighted;
    for target_seqLength=1:pointKindsNum
        ax = nexttile;
        set(gca,'color',backgroundColor);
        
        % Plot frames as layout
        for tempj=1:numFrames
            showRadius;
            x = layout_x(tempj);
            y = layout_y(tempj);
            r = showRadius;
            pos = [x-r, y-r, 2*r, 2*r];
            %temp_color = [0 0 0];
            
            if isnan(boost_weighted{target_seqLength}(tempj)) == 0 %#ok<*COMPNOT>
                %temp_colorIndex = round((boost_weighted{target_seqLength}(tempj)-temp_min)/(temp_max-temp_min)*255)+1;
                temp_colorIndex = ceil((boost_weighted{target_seqLength}(tempj)-temp_min)/(temp_max-temp_min)*(size(my_color,1)-1))+1;
                temp_color = my_color(temp_colorIndex, :);
                
                if boost_weighted{target_seqLength}(tempj)<0
                    %temp_edgeColor =
                    LineWidth = 4;%5
                    %EdgeColor = '#F0C7A2';
                else
                    LineWidth = 1;
                    %EdgeColor = [0 0 0];
                end
                %EdgeColor = [0.6 0.6 0.6];
                %EdgeColor = [0 0 0];
                EdgeColor = '#F0C7A2';
                
                %                 if boost_weighted{target_seqLength}(tempj)<0
                %                     %temp_edgeColor =
                %                     LineWidth = 5;
                %                 else
                %                     LineWidth = 1;
                %                 end
                %                 if if_fig8_boost_negative0_positive1_full2 == 2
                %                     %EdgeColor = [0.6 0.6 0.6];
                %                     EdgeColor = [0 0 0];
                %                 elseif if_fig8_boost_negative0_positive1_full2 == 1
                %                     LineWidth = 1;
                %                     %EdgeColor = backgroundColor*1.5;
                %                     EdgeColor = [1 1 1];
                %                 elseif if_fig8_boost_negative0_positive1_full2 == 0
                %                     if if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 ~= 2
                %                         %EdgeColor = [0 0 0];
                %                         %EdgeColor = [1 1 1];
                %                         EdgeColor = [0.6 0.6 0.6];
                %                     else
                %                         LineWidth = 1;
                %                         EdgeColor = [1 1 1];
                %                     end
                %                 end
                
                %rectangle('Position',pos,'Curvature',[1 1],'FaceColor',temp_color,'LineWidth',LineWidth,'EdgeColor',[0.6 0.6 0.6]);
                rectangle('Position',pos,'Curvature',[1 1],'FaceColor',temp_color,'LineWidth',LineWidth,'EdgeColor',EdgeColor);
            elseif isnan(boost_weighted{target_seqLength}(tempj)) == 1 %#ok<*COMPNOP>
                LineWidth = 1;
                FaceColor = backgroundColor;
                EdgeColor = '#F9E9DB';
                %EdgeColor = [0 0 0];
                %                 if if_fig8_boost_negative0_positive1_full2 == 2
                %                     %EdgeColor = [0.6 0.6 0.6];
                %                     EdgeColor = [0 0 0];
                %                 else
                %                     %EdgeColor = [1 1 1];
                %                     EdgeColor = [0 0 0];
                %                 end
                %rectangle('Position',pos,'Curvature',[1 1],'FaceColor',FaceColor,'LineWidth',LineWidth,'EdgeColor',EdgeColor);
                %rectangle('Position',pos,'Curvature',[1 1],'FaceColor','#F9E9DB','LineWidth',LineWidth,'LineStyle', '--','EdgeColor',EdgeColor);
                %rectangle('Position',pos,'Curvature',[1 1],'FaceColor','#F9E9DB','LineWidth',LineWidth,'EdgeColor',EdgeColor);
                hold on
                plot(layout_x(tempj),layout_y(tempj),'x','LineWidth',2,'MarkerSize',20,'MarkerEdgeColor',[0 0 0]);
            end
        end
        % Plot central cross
        layout_x;
        layout_y;
        hold on
        plot(layout_x(end),layout_y(end),'+','LineWidth',1,'MarkerSize',20,'MarkerEdgeColor',[0 0 0]);
        
        ax.YAxis.Visible = 'off';
        ax.XAxis.Visible = 'off';
        axis equal
        set(gca,'YDir','reverse');% Reverse the direction of y axis
        xlim([min(layout_x)-150 max(layout_x)+150]);
        ylim([min(layout_y)-150 max(layout_y)+150]);
        
        title_str = ['Length=',num2str(target_seqLength)];
        if if_fig8_boost_negative0_positive1_full2 == 2
            title_str = [title_str,', full boost'];
        elseif if_fig8_boost_negative0_positive1_full2 == 0
            title_str = [title_str,', negative boost'];
            if if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 2
                title_str = [title_str,', full error types'];
            elseif if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 0
                title_str = [title_str,', only delete/add one'];
            elseif if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 1
                title_str = [title_str,', only interference'];
            end
        elseif if_fig8_boost_negative0_positive1_full2 == 1
            title_str = [title_str,', positive boost'];
        end
        title(title_str);
        
    end
    
    %colorbar('Ticks', [-1:0.2:1]);
    
    temp_stepNum = 6;
    temp_step = (temp_max-temp_min)/(temp_stepNum-1) ;
    temp_ticks = [temp_min:temp_step:temp_max];
    temp_ticksLabels = cell(1, temp_stepNum);
    for tempi=1:temp_stepNum
        %temp_ticksLabels{tempi} = num2str(temp_ticks(tempi));
        temp_ticksLabels{tempi} = sprintf('%.3f',temp_ticks(tempi));
    end
    colorbar('TickLabels',temp_ticksLabels)
    
    if if_fig8_boost_negative0_positive1_full2 == 2
        currentFigName = sprintf('BoostHeatmap_eachLength_eachLocation_');
    elseif if_fig8_boost_negative0_positive1_full2 == 0
        if if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 2
            currentFigName = sprintf('BoostHeatmap_negative_eachLength_eachLocation_');
        elseif if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 0
            currentFigName = sprintf('BoostHeatmap_negative_eachLength_eachLocation_onlyDeleteAddone_');
        elseif if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 1
            currentFigName = sprintf('BoostHeatmap_negative_eachLength_eachLocation_onlyInterference_');
        end
        
    elseif if_fig8_boost_negative0_positive1_full2 == 1
        currentFigName = sprintf('BoostHeatmap_positive_eachLength_eachLocation_');
    end
    
    % to generate a unique file name for saving figure
    fileName_fig8 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(8) == 1
        exportgraphics(fig8, fileName_fig8, 'BackgroundColor', 'none', 'ContentType', 'vector');
    end
end

%% Left offset seq vs right offset seq, comparing by boost effect
temp_index_collapsed;
choiceMinusNoChoice_2;

temp_centroid = cell(1, pointKindsNum);
temp_midline = (1+numFrames)/2;
temp_if_rightOffset = cell(1, pointKindsNum);

for tempi=1:pointKindsNum
    temp_centroid{tempi} = zeros(1, numSeq(tempi));
    temp_if_rightOffset{tempi} = nan(1, numSeq(tempi));
    for tempj=1:numSeq(tempi)
        temp_seq = target_seqSet{tempi}{tempj};
        temp_centroid{tempi}(tempj) = mean(temp_seq);
        if mean(temp_seq) > temp_midline
            temp_if_rightOffset{tempi}(tempj) = true;
        elseif mean(temp_seq) < temp_midline
            temp_if_rightOffset{tempi}(tempj) = false;
        end
    end
end
temp_if_rightOffset;% Seqs which are exactly of middle line is nan.
temp_index_collapsed;
choiceMinusNoChoice_2;
choiceMinusNoChoice_offset = struct;
choiceMinusNoChoice_offset.leftOffset = cell(1, pointKindsNum);
choiceMinusNoChoice_offset.rightOffset = cell(1, pointKindsNum);

for tempi=1:pointKindsNum
    choiceMinusNoChoice_offset.leftOffset{tempi};
    for tempj=1:length(temp_index_collapsed{tempi})
        tempk = temp_index_collapsed{tempi}(tempj);
        
        if temp_if_rightOffset{tempi}(tempk) == false
            choiceMinusNoChoice_offset.leftOffset{tempi} = ...
                [choiceMinusNoChoice_offset.leftOffset{tempi} choiceMinusNoChoice_2{tempi}(tempk)];
            
        elseif temp_if_rightOffset{tempi}(tempk) == true
            choiceMinusNoChoice_offset.rightOffset{tempi} = ...
                [choiceMinusNoChoice_offset.rightOffset{tempi} choiceMinusNoChoice_2{tempi}(tempk)];
            
            %         elseif isnan(temp_if_rightOffset{tempi}(tempk))
            %             choiceMinusNoChoice_offset.leftOffset{tempi} = ...
            %                 [choiceMinusNoChoice_offset.leftOffset{tempi} choiceMinusNoChoice_2{tempi}(tempk)];
            %             choiceMinusNoChoice_offset.rightOffset{tempi} = ...
            %                 [choiceMinusNoChoice_offset.rightOffset{tempi} choiceMinusNoChoice_2{tempi}(tempk)];
        end
    end
end
choiceMinusNoChoice_offset;

choiceMinusNoChoice_offset_inOne = struct;
choiceMinusNoChoice_offset_inOne.leftOffset = [];
choiceMinusNoChoice_offset_inOne.rightOffset = [];
for tempi=1:pointKindsNum
    choiceMinusNoChoice_offset_inOne.leftOffset = ...
        [choiceMinusNoChoice_offset_inOne.leftOffset choiceMinusNoChoice_offset.leftOffset{tempi}];
    choiceMinusNoChoice_offset_inOne.rightOffset = ...
        [choiceMinusNoChoice_offset_inOne.rightOffset choiceMinusNoChoice_offset.rightOffset{tempi}];
end
choiceMinusNoChoice_offset_inOne.max = max(...
    [choiceMinusNoChoice_offset_inOne.leftOffset, choiceMinusNoChoice_offset_inOne.rightOffset]);
choiceMinusNoChoice_offset_inOne.min = min(...
    [choiceMinusNoChoice_offset_inOne.leftOffset, choiceMinusNoChoice_offset_inOne.rightOffset]);

%% Figure16
if ifFigure16 == 1
    fig16 = figure(16);
    set(gcf,'Position',[0 50 800 600]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout('flow','TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    
    backgroundColor = [1 1 1];
    
    p_t_offset = zeros(1, pointKindsNum);
    for tempi=1:pointKindsNum
        if isempty(choiceMinusNoChoice_offset.leftOffset{tempi}) ||...
                isempty(choiceMinusNoChoice_offset.rightOffset{tempi})
            continue
        end
        ax = nexttile;
        set(gca,'color',backgroundColor);
        
        [h_t_offset,p_t_offset(tempi)]=ttest2(choiceMinusNoChoice_offset.leftOffset{tempi}, ...
            choiceMinusNoChoice_offset.rightOffset{tempi});
        
        temp_data = [choiceMinusNoChoice_offset.leftOffset{tempi}, choiceMinusNoChoice_offset.rightOffset{tempi}];
        temp_category1 = 'Left offset';
        temp_category2 = 'Right offset';
        temp_label = cell(1, length(temp_data));
        for tempj=1:length(temp_data)
            if tempj <= length(choiceMinusNoChoice_offset.leftOffset{tempi})
                temp_label{tempj} = temp_category1;
            else
                temp_label{tempj} = temp_category2;
            end
        end
        
        violinplot(temp_data, temp_label);
        hold on
        
        tempTxt = sprintf('');
        if p_t_offset(tempi) < 0.001/pointKindsNum
            tempTxt = sprintf('***');
        elseif p_t_offset(tempi) < 0.01/pointKindsNum
            tempTxt = sprintf('**');
        elseif p_t_offset(tempi) < 0.05/pointKindsNum
            tempTxt = sprintf('*');
        end
        
        text(1.5,choiceMinusNoChoice_offset_inOne.max-0.08,[tempTxt sprintf('p = %.4f',p_t_offset(tempi))], 'fontsize',16,'FontWeight','bold');
        
        xlim([0.5 2.5]);
        ylim([choiceMinusNoChoice_offset_inOne.min choiceMinusNoChoice_offset_inOne.max]);
        set(gca, 'FontSize', 16)
        set(gca,'box','off');% 取消右、上边框
        ylabel('Accuracy', 'FontSize', 16, 'FontWeight', 'bold');
        
        
        
        title_str = ['Length=',num2str(tempi)];
        if if_fig8_boost_negative0_positive1_full2 == 2
            title_str = [title_str,', full boost'];
        elseif if_fig8_boost_negative0_positive1_full2 == 0
            title_str = [title_str,', negative boost'];
            if if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 2
                title_str = [title_str,', full error types'];
            elseif if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 0
                title_str = [title_str,', only delete/add one'];
            elseif if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 1
                title_str = [title_str,', only interference'];
            end
        elseif if_fig8_boost_negative0_positive1_full2 == 1
            title_str = [title_str,', positive boost'];
        end
        title(title_str,'FontSize',12);
        
    end
    
    
    
    if if_fig8_boost_negative0_positive1_full2 == 2
        currentFigName = sprintf('BoostHeatmap_eachLength_eachLocation_');
    elseif if_fig8_boost_negative0_positive1_full2 == 0
        if if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 2
            currentFigName = sprintf('BoostHeatmap_negative_eachLength_eachLocation_');
        elseif if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 0
            currentFigName = sprintf('BoostHeatmap_negative_eachLength_eachLocation_onlyDeleteAddone_');
        elseif if_negativeBoostWithErrorType_deleteAddOne0_interference1_full2 == 1
            currentFigName = sprintf('BoostHeatmap_negative_eachLength_eachLocation_onlyInterference_');
        end
        
    elseif if_fig8_boost_negative0_positive1_full2 == 1
        currentFigName = sprintf('BoostHeatmap_positive_eachLength_eachLocation_');
    end
    
    % to generate a unique file name for saving figure
    fileName_fig16 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(16) == 1
        exportgraphics(fig16, fileName_fig16, 'BackgroundColor', 'none', 'ContentType', 'vector');
    end
end


%% Compute and plot no choice error types in length level
valid_index;
valid_boostMoreThan0_index = (choiceMinusNoChoice_inOne'>0) & valid_index;
valid_boostLessThan0_index = (choiceMinusNoChoice_inOne'<0) & valid_index;


length_rangeHead = seq_length_rangeHead;
length_rangeTail = seq_length_rangeTail;

% temp_mean = zeros(pointKindsNum, 3);
% temp_SEM = zeros(pointKindsNum, 3);

temp_mean = struct;
temp_mean.boostMoreThan0 = zeros(pointKindsNum, 3);
temp_mean.boostLessThan0 = zeros(pointKindsNum, 3);
temp_SEM = struct;
temp_SEM.boostMoreThan0 = zeros(pointKindsNum, 3);
temp_SEM.boostLessThan0 = zeros(pointKindsNum, 3);

for target_seqLength=1:pointKindsNum
    if target_seqLength == 1
        temp_range = 1:numSeq(1);
    else
        temp_range = sum(numSeq(1:target_seqLength-1))+1:sum(numSeq(1:target_seqLength));
    end
    for temp_conditionIndex=1:2
        if temp_conditionIndex == 1
            temp_valid = valid_boostMoreThan0_index(temp_range);
            temp_str = 'boostMoreThan0';
        elseif temp_conditionIndex == 2
            temp_valid = valid_boostLessThan0_index(temp_range);
            temp_str = 'boostLessThan0';
        end
        
        temp_deleteAddOne = errorType_seqSummary_collapsed.noChoice.deleteAddOne{target_seqLength}(temp_valid);
        temp_interference = errorType_seqSummary_collapsed.noChoice.interference{target_seqLength}(temp_valid);
        temp_others = errorType_seqSummary_collapsed.noChoice.others{target_seqLength}(temp_valid);
        
        temp_mean.(temp_str)(target_seqLength, 1) = mean(temp_deleteAddOne);
        temp_mean.(temp_str)(target_seqLength, 2) = mean(temp_interference);
        temp_mean.(temp_str)(target_seqLength, 3) = mean(temp_others);
        
        if sum(temp_valid) > 0
            temp_SEM.(temp_str)(target_seqLength, 1) = std(temp_deleteAddOne) / sqrt(sum(temp_valid));
            temp_SEM.(temp_str)(target_seqLength, 2) = std(temp_interference) / sqrt(sum(temp_valid));
            temp_SEM.(temp_str)(target_seqLength, 3) = std(temp_others) / sqrt(sum(temp_valid));
        end
        
    end
end
temp_length = length_rangeHead:length_rangeTail;


temp_color = {'#0072BD';'#D95319';'#7E2F8E'};

%% Figure 17
if ifFigure17 == 1
    fig17 = figure(17);
    set(gcf,'Position',[0 50 1200 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
    
    for temp_conditionIndex=1:2
        nexttile
        
        if temp_conditionIndex == 1
            temp_str = 'boostMoreThan0';
            temp_titleAppend = 'with seqs that boost>0';
        elseif temp_conditionIndex == 2
            temp_str = 'boostLessThan0';
            temp_titleAppend = 'with seqs that boost<0';
        end
        length_valid_range = ~isnan(temp_mean.(temp_str)(:,1))';
        
        for tempi=1:3
            errorbar(temp_length(length_valid_range),temp_mean.(temp_str)(length_valid_range, tempi),temp_SEM.(temp_str)(length_valid_range, tempi),...
                '-o','Color', temp_color{tempi},'LineWidth',1,'CapSize',12,'MarkerSize',5);
            hold on
        end
        
        plot([seq_length_rangeHead-1 seq_length_rangeHead: seq_length_rangeTail seq_length_rangeTail+1]...
            , zeros(1, pointKindsNum+2), '--');
        
        legend('delete/add One','interference','others','Location','northwest','fontsize',11)
        
        %ylim([-0.4 0.4]);
        ylim([0 0.5]);
        set(gca, 'FontSize', 20)
        %set(gca,'YTick',[-0.4:0.2:0.4]);%设置要显示坐标刻度的范围
        set(gca,'XLim',[length_rangeHead-1 length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XTick',[length_rangeHead:1:length_rangeTail]);%设置要显示坐标刻度的范围
        set(gca,'XTickLabel',[length_rangeHead:1:length_rangeTail]);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        xlabel('Length', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Proportion', 'FontSize', 18, 'FontWeight', 'bold');
        title({'\fontsize{20}{\bf No choice error types in length}';temp_titleAppend});
        
    end
end

a=seqSet_inOne_inOne(valid_boostMoreThan0_index)';
a2=seqSet_inOne_inOne(valid_boostLessThan0_index)';

%% Compute delete one, add one and interference proportion with each order in each length
current_responseMatrix = responseMatrix.noChoice_n11n;
valid_index;
errorType_eachOrder_eachLength = struct;
errorTransitionMatrix_eachLength = struct;
errorTransitionMatrix_eachLength_validCount = struct;

errorType_eachOrder_eachLength_negative = struct;
errorTransitionMatrix_eachLength_negative = struct;
errorTransitionMatrix_eachLength_validCount_negative = struct;

% for delete one
errorType_eachOrder_eachLength.deleteOne = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength.deleteOne = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength_validCount.deleteOne = cell(1, pointKindsNum);

errorType_eachOrder_eachLength_negative.deleteOne = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength_negative.deleteOne = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength_validCount_negative.deleteOne = cell(1, pointKindsNum);
for tempi=1:pointKindsNum
    if tempi == 1
        continue
    end
    errorType_eachOrder_eachLength.deleteOne{tempi} = zeros(1, tempi);
    errorTransitionMatrix_eachLength.deleteOne{tempi} = zeros(tempi, numFrames);
    errorTransitionMatrix_eachLength_validCount.deleteOne{tempi} = zeros(tempi, numFrames);
    
    errorType_eachOrder_eachLength_negative.deleteOne{tempi} = zeros(1, tempi);
    errorTransitionMatrix_eachLength_negative.deleteOne{tempi} = zeros(tempi, numFrames);
    errorTransitionMatrix_eachLength_validCount_negative.deleteOne{tempi} = zeros(tempi, numFrames);
    response_length = tempi - 1;
    for tempj=1:numSeq(tempi)
        stimuli_seq = target_seqSet{tempi}{tempj};
        stimuli_index = sum(numSeq(1:(tempi-1))) + tempj;
        if valid_index(stimuli_index) == 0
            continue
        end
        for delete_order=tempi:-1:1
            %errorType_eachOrder_eachLength.deleteOne{tempi}(delete_order);
            
            response_seq = stimuli_seq;
            response_seq(delete_order) = [];
            response_index = 0;
            for tempk=sum(numSeq(1:(response_length-1)))+1:sum(numSeq(1:response_length))
                
                target_seqSet_inOne{tempk};
                if ismember(response_seq, target_seqSet_inOne{tempk}) == 1
                    response_index = tempk;
                    break;
                end
            end
            stimuli_index;
            response_index;
            errorType_eachOrder_eachLength.deleteOne{tempi}(delete_order) = ...
                errorType_eachOrder_eachLength.deleteOne{tempi}(delete_order) + ...
                current_responseMatrix(stimuli_index, response_index);
            
            errorTransitionMatrix_eachLength_validCount.deleteOne{tempi}(delete_order,stimuli_seq(delete_order)) = ...
                errorTransitionMatrix_eachLength_validCount.deleteOne{tempi}(delete_order,stimuli_seq(delete_order)) + 1;
            
            errorTransitionMatrix_eachLength.deleteOne{tempi}(delete_order, stimuli_seq(delete_order)) = ...
                errorTransitionMatrix_eachLength.deleteOne{tempi}(delete_order, stimuli_seq(delete_order)) + ...
                current_responseMatrix(stimuli_index, response_index);
            
            if valid_boostLessThan0_index(stimuli_index) == 1
                errorType_eachOrder_eachLength_negative.deleteOne{tempi}(delete_order) = ...
                    errorType_eachOrder_eachLength_negative.deleteOne{tempi}(delete_order) + ...
                    current_responseMatrix(stimuli_index, response_index);
                
                errorTransitionMatrix_eachLength_validCount_negative.deleteOne{tempi}(delete_order,stimuli_seq(delete_order)) = ...
                    errorTransitionMatrix_eachLength_validCount_negative.deleteOne{tempi}(delete_order,stimuli_seq(delete_order)) + 1;
                
                errorTransitionMatrix_eachLength_negative.deleteOne{tempi}(delete_order, stimuli_seq(delete_order)) = ...
                    errorTransitionMatrix_eachLength_negative.deleteOne{tempi}(delete_order, stimuli_seq(delete_order)) + ...
                    current_responseMatrix(stimuli_index, response_index);
            end
            
        end
    end
    if tempi == 1
        temp_range = 1:numSeq(1);
    else
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    end
    errorType_eachOrder_eachLength.deleteOne{tempi} = ...
        errorType_eachOrder_eachLength.deleteOne{tempi} ./ sum(valid_index(temp_range));
    
    errorTransitionMatrix_eachLength.deleteOne{tempi} = ...
        errorTransitionMatrix_eachLength.deleteOne{tempi} ./ ...
        errorTransitionMatrix_eachLength_validCount.deleteOne{tempi};
    
    temp_range;
    valid_boostLessThan0_range = find(valid_boostLessThan0_index==1);
    temp_validRange_negative = valid_boostLessThan0_range(ismember(valid_boostLessThan0_range,temp_range));
    
    errorType_eachOrder_eachLength_negative.deleteOne{tempi} = ...
        errorType_eachOrder_eachLength_negative.deleteOne{tempi} ./ length(temp_validRange_negative);
    
    errorTransitionMatrix_eachLength_negative.deleteOne{tempi} = ...
        errorTransitionMatrix_eachLength_negative.deleteOne{tempi} ./ ...
        errorTransitionMatrix_eachLength_validCount_negative.deleteOne{tempi};
    
end

% for add one
errorType_eachOrder_eachLength.addOne = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength.addOne = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength_validCount.addOne = cell(1, pointKindsNum);

errorType_eachOrder_eachLength_negative.addOne = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength_negative.addOne = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength_validCount_negative.addOne = cell(1, pointKindsNum);
for tempi=1:pointKindsNum
    response_length = tempi + 1;
    errorType_eachOrder_eachLength.addOne{tempi} = zeros(1, response_length);
    errorTransitionMatrix_eachLength.addOne{tempi} = zeros(tempi+1, numFrames);
    errorTransitionMatrix_eachLength_validCount.addOne{tempi} = zeros(tempi+1, numFrames);
    
    errorType_eachOrder_eachLength_negative.addOne{tempi} = zeros(1, response_length);
    errorTransitionMatrix_eachLength_negative.addOne{tempi} = zeros(tempi+1, numFrames);
    errorTransitionMatrix_eachLength_validCount_negative.addOne{tempi} = zeros(tempi+1, numFrames);
    for tempj=1:numSeq(tempi)
        stimuli_seq = target_seqSet{tempi}{tempj};
        target_seqSet_inOne{stimuli_index};
        stimuli_index = sum(numSeq(1:(tempi-1))) + tempj;
        if valid_index(stimuli_index) == 0
            continue
        end
        
        add_candidate = 1:numFrames;
        add_candidate(stimuli_seq) = [];
        for add_index=1:length(add_candidate)
            add_candidate(add_index);
            add_order = sum(stimuli_seq<add_candidate(add_index))+1;
            response_seq = sort([stimuli_seq add_candidate(add_index)]);
            response_index = 0;
            for tempk=sum(numSeq_extend(1:(response_length-1)))+1:sum(numSeq_extend(1:response_length))
                
                target_seqSet_extend_inOne{tempk}; %#ok<*VUNUS>
                if ismember(response_seq, target_seqSet_extend_inOne{tempk}) == 1
                    response_index = tempk;
                    break;
                end
            end
            stimuli_index;
            response_index;
            errorType_eachOrder_eachLength.addOne{tempi}(add_order) = ...
                errorType_eachOrder_eachLength.addOne{tempi}(add_order) + ...
                current_responseMatrix(stimuli_index, response_index);
            
            errorTransitionMatrix_eachLength_validCount.addOne{tempi}(add_order, add_candidate(add_index)) = ...
                errorTransitionMatrix_eachLength_validCount.addOne{tempi}(add_order, add_candidate(add_index)) + 1;
            
            errorTransitionMatrix_eachLength.addOne{tempi}(add_order, add_candidate(add_index)) = ...
                errorTransitionMatrix_eachLength.addOne{tempi}(add_order, add_candidate(add_index)) + ...
                current_responseMatrix(stimuli_index, response_index);
            
            if valid_boostLessThan0_index(stimuli_index) == 1
                errorType_eachOrder_eachLength_negative.addOne{tempi}(add_order) = ...
                    errorType_eachOrder_eachLength_negative.addOne{tempi}(add_order) + ...
                    current_responseMatrix(stimuli_index, response_index);
                
                errorTransitionMatrix_eachLength_validCount_negative.addOne{tempi}(add_order, add_candidate(add_index)) = ...
                    errorTransitionMatrix_eachLength_validCount_negative.addOne{tempi}(add_order, add_candidate(add_index)) + 1;
                
                errorTransitionMatrix_eachLength_negative.addOne{tempi}(add_order, add_candidate(add_index)) = ...
                    errorTransitionMatrix_eachLength_negative.addOne{tempi}(add_order, add_candidate(add_index)) + ...
                    current_responseMatrix(stimuli_index, response_index);
            end
            
        end
    end
    if tempi == 1
        temp_range = 1:numSeq(1);
    else
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    end
    errorType_eachOrder_eachLength.addOne{tempi} = ...
        errorType_eachOrder_eachLength.addOne{tempi} / sum(valid_index(temp_range));
    
    errorTransitionMatrix_eachLength.addOne{tempi} = ...
        errorTransitionMatrix_eachLength.addOne{tempi} ./ ...
        errorTransitionMatrix_eachLength_validCount.addOne{tempi};
    
    temp_range;
    valid_boostLessThan0_range = find(valid_boostLessThan0_index==1);
    temp_validRange_negative = valid_boostLessThan0_range(ismember(valid_boostLessThan0_range,temp_range));
    
    errorType_eachOrder_eachLength_negative.addOne{tempi} = ...
        errorType_eachOrder_eachLength_negative.addOne{tempi} ./ length(temp_validRange_negative);
    
    errorTransitionMatrix_eachLength_negative.addOne{tempi} = ...
        errorTransitionMatrix_eachLength_negative.addOne{tempi} ./ ...
        errorTransitionMatrix_eachLength_validCount_negative.addOne{tempi};
    
end

% for interference
errorType_indexMatrix.interference;
errorType_eachOrder_eachLength.interference = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength.interference = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength_validCount.interference = cell(1, pointKindsNum);

errorType_eachOrder_eachLength_negative.interference = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength_negative.interference = cell(1, pointKindsNum);
errorTransitionMatrix_eachLength_validCount_negative.interference = cell(1, pointKindsNum);
for tempi=1:pointKindsNum
    errorType_eachOrder_eachLength.interference{tempi} = zeros(1, tempi);
    errorTransitionMatrix_eachLength.interference{tempi} = cell(1, tempi);
    errorTransitionMatrix_eachLength_validCount.interference{tempi} = cell(1, tempi);
    
    errorType_eachOrder_eachLength_negative.interference{tempi} = zeros(1, tempi);
    errorTransitionMatrix_eachLength_negative.interference{tempi} = cell(1, tempi);
    errorTransitionMatrix_eachLength_validCount_negative.interference{tempi} = cell(1, tempi);
    for temp_order=1:tempi
        errorTransitionMatrix_eachLength.interference{tempi}{temp_order} = zeros(numFrames, numFrames);
        errorTransitionMatrix_eachLength_validCount.interference{tempi}{temp_order} = zeros(numFrames, numFrames);
        
        errorTransitionMatrix_eachLength_negative.interference{tempi}{temp_order} = zeros(numFrames, numFrames);
        errorTransitionMatrix_eachLength_validCount_negative.interference{tempi}{temp_order} = zeros(numFrames, numFrames);
    end
    
    for tempj=1:numSeq(tempi)
        stimuli_seq = target_seqSet{tempi}{tempj};
        target_seqSet_inOne{stimuli_index};
        stimuli_index = sum(numSeq(1:(tempi-1))) + tempj;
        if valid_index(stimuli_index) == 0
            continue
        end
        
        response_index_candidate = find(errorType_indexMatrix.interference(stimuli_index, :)==1);
        %current_responseMatrix(stimuli_index, response_index_candidate);
        for interference_index=1:length(response_index_candidate)
            response_index = response_index_candidate(interference_index);
            response_seq = target_seqSet_inOne{response_index};
            %interference_order = find(ismember(stimuli_seq, response_seq)==0);
            interference_order = [];
            for temp_order=1:tempi
                if stimuli_seq(temp_order) ~= response_seq(temp_order)
                    interference_order = [interference_order temp_order];
                end
            end
            
            stimuli_index;
            response_index;
            for tempk=1:length(interference_order)
                interference_order(tempk);
                errorType_eachOrder_eachLength.interference{tempi}(interference_order(tempk)) = ...
                    errorType_eachOrder_eachLength.interference{tempi}(interference_order(tempk)) + ...
                    current_responseMatrix(stimuli_index, response_index);
                
                errorTransitionMatrix_eachLength_validCount.interference{tempi}{interference_order(tempk)}(stimuli_seq(interference_order(tempk)), response_seq(interference_order(tempk))) = ...
                    errorTransitionMatrix_eachLength_validCount.interference{tempi}{interference_order(tempk)}(stimuli_seq(interference_order(tempk)), response_seq(interference_order(tempk))) + 1;
                
                errorTransitionMatrix_eachLength.interference{tempi}{interference_order(tempk)}(stimuli_seq(interference_order(tempk)), response_seq(interference_order(tempk))) = ...
                    errorTransitionMatrix_eachLength.interference{tempi}{interference_order(tempk)}(stimuli_seq(interference_order(tempk)), response_seq(interference_order(tempk))) + ...
                    current_responseMatrix(stimuli_index, response_index);
                
                if valid_boostLessThan0_index(stimuli_index) == 1
                    errorType_eachOrder_eachLength_negative.interference{tempi}(interference_order(tempk)) = ...
                        errorType_eachOrder_eachLength_negative.interference{tempi}(interference_order(tempk)) + ...
                        current_responseMatrix(stimuli_index, response_index);
                    
                    errorTransitionMatrix_eachLength_validCount_negative.interference{tempi}{interference_order(tempk)}(stimuli_seq(interference_order(tempk)), response_seq(interference_order(tempk))) = ...
                        errorTransitionMatrix_eachLength_validCount_negative.interference{tempi}{interference_order(tempk)}(stimuli_seq(interference_order(tempk)), response_seq(interference_order(tempk))) + 1;
                    
                    errorTransitionMatrix_eachLength_negative.interference{tempi}{interference_order(tempk)}(stimuli_seq(interference_order(tempk)), response_seq(interference_order(tempk))) = ...
                        errorTransitionMatrix_eachLength_negative.interference{tempi}{interference_order(tempk)}(stimuli_seq(interference_order(tempk)), response_seq(interference_order(tempk))) + ...
                        current_responseMatrix(stimuli_index, response_index);
                end
                
            end
        end
    end
    if tempi == 1
        temp_range = 1:numSeq(1);
    else
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    end
    errorType_eachOrder_eachLength.interference{tempi} = ...
        errorType_eachOrder_eachLength.interference{tempi} / sum(valid_index(temp_range));
    
    errorTransitionMatrix_eachLength.interference{tempi};
    for tempj=1:tempi
        errorTransitionMatrix_eachLength.interference{tempi}{tempj} = ...
            errorTransitionMatrix_eachLength.interference{tempi}{tempj} ./ ...
            errorTransitionMatrix_eachLength_validCount.interference{tempi}{tempj};
    end
    
    temp_range;
    valid_boostLessThan0_range = find(valid_boostLessThan0_index==1);
    temp_validRange_negative = valid_boostLessThan0_range(ismember(valid_boostLessThan0_range,temp_range));
    
    errorType_eachOrder_eachLength_negative.interference{tempi} = ...
        errorType_eachOrder_eachLength_negative.interference{tempi} ./ length(temp_validRange_negative);
    
    for tempj=1:tempi
        errorTransitionMatrix_eachLength_negative.interference{tempi}{tempj} = ...
            errorTransitionMatrix_eachLength_negative.interference{tempi}{tempj} ./ ...
            errorTransitionMatrix_eachLength_validCount_negative.interference{tempi}{tempj};
    end
    
end

for temp_conditionIndex=1:3
    if temp_conditionIndex == 1
        temp_str = 'deleteOne';
    elseif temp_conditionIndex == 2
        temp_str = 'addOne';
    elseif temp_conditionIndex == 3
        temp_str = 'interference';
    end
    for tempi=1:pointKindsNum
        if isempty(errorType_eachOrder_eachLength.(temp_str){tempi})
            errorType_eachOrder_eachLength.(temp_str){tempi} = nan;
        end
        if isempty(errorType_eachOrder_eachLength_negative.(temp_str){tempi})
            errorType_eachOrder_eachLength_negative.(temp_str){tempi} = nan;
        end
    end
end

%% Figure 18
if ifFigure18 == 1
    fig18 = figure(18);
    set(gcf,'Position',[0 50 1650 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    
    temp_color = [0.00 0.00 0.00;
        0.28 0.28 0.28;
        0.56 0.56 0.56;
        0.84 0.84 0.84];
    a = 1;
    length_rangeHead = seq_length_rangeHead;
    length_rangeTail = seq_length_rangeTail+1;
    
    for temp_conditionIndex=1:3
        nexttile
        
        if temp_conditionIndex == 1
            temp_str = 'deleteOne';
            temp_titleAppend = 'delete one';
        elseif temp_conditionIndex == 2
            temp_str = 'addOne';
            temp_titleAppend = 'add one';
        elseif temp_conditionIndex == 3
            temp_str = 'interference';
            temp_titleAppend = 'interference';
        end
        
        for target_seqLength=1:pointKindsNum
            %if isempty(errorType_eachOrder_eachLength.(temp_str){target_seqLength})
            %    errorType_eachOrder_eachLength.(temp_str){target_seqLength} = nan;
            %end
            temp_rangeTail = length(errorType_eachOrder_eachLength.(temp_str){target_seqLength});
            plot(1:temp_rangeTail,errorType_eachOrder_eachLength.(temp_str){target_seqLength},...
                '-o','Color',temp_color(target_seqLength,:),'LineWidth',1,'MarkerSize',5);
            hold on
        end
        
        legend('length=1','length=2','length=3','length=4',...
            'Location','northeast','fontsize',9)
        
        % ylim([min(temp_mean_collaped)-0.1 1]);
        set(gca, 'FontSize', 20)
        % set(gca,'YTick',[-0.3:0.1:0.3]);%设置要显示坐标刻度的范围
        set(gca,'XLim',[length_rangeHead-1 length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XTick',[length_rangeHead:1:length_rangeTail]);%设置要显示坐标刻度的范围
        set(gca,'XTickLabel',[length_rangeHead:1:length_rangeTail]);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        xlabel('Order', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Proportion', 'FontSize', 18, 'FontWeight', 'bold');
        title({'\fontsize{20}{\bf Error types with each order}';temp_titleAppend});
    end
    currentFigName = ['errorTypesWithEachOrder', '_'];
    % to generate a unique file name for saving figure
    fileName_fig18 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(18) == 1
        exportgraphics(fig18, fileName_fig18, 'BackgroundColor', 'none');
    end
end

errorType_eachOrder_eachLength;
errorType_dominantOrder_eachLength = errorType_eachOrder_eachLength;
errorType_dominantOrder_eachLength_negative = errorType_eachOrder_eachLength_negative;

for temp_conditionIndex=1:3
    if temp_conditionIndex == 1
        temp_str = 'deleteOne';
    elseif temp_conditionIndex == 2
        temp_str = 'addOne';
    elseif temp_conditionIndex == 3
        temp_str = 'interference';
    end
    errorType_eachOrder_eachLength.(temp_str);
    errorType_dominantOrder_eachLength.(temp_str) = nan(1, length(errorType_eachOrder_eachLength.(temp_str)));
    errorType_dominantOrder_eachLength_negative.(temp_str) = nan(1, length(errorType_eachOrder_eachLength_negative.(temp_str)));
    for tempi=1:pointKindsNum
        temp_errorType = errorType_eachOrder_eachLength.(temp_str){tempi};
        if ~isnan(temp_errorType)
            [temp_max, temp_maxIndex] = max(temp_errorType);
            errorType_dominantOrder_eachLength.(temp_str)(tempi) = temp_maxIndex;
        end
        temp_errorType_negative = errorType_eachOrder_eachLength_negative.(temp_str){tempi};
        if ~isnan(temp_errorType_negative)
            [temp_max, temp_maxIndex] = max(temp_errorType_negative); %#ok<*ASGLU>
            errorType_dominantOrder_eachLength_negative.(temp_str)(tempi) = temp_maxIndex;
        end
    end
end
errorType_dominantOrder_eachLength;
errorTransitionMatrix_eachLength;

errorTransitionMatrix_dominantOrder_eachLength = errorType_eachOrder_eachLength;
errorTransitionMatrix_dominantOrder_eachLength_negative = errorType_eachOrder_eachLength_negative;
for temp_conditionIndex=1:3
    if temp_conditionIndex == 1
        temp_str = 'deleteOne';
    elseif temp_conditionIndex == 2
        temp_str = 'addOne';
    elseif temp_conditionIndex == 3
        temp_str = 'interference';
    end
    errorType_dominantOrder_eachLength.(temp_str);
    for tempi=1:pointKindsNum
        temp_dominantOrder = errorType_dominantOrder_eachLength.(temp_str)(tempi);
        if isnan(temp_dominantOrder)
            continue
        end
        if temp_conditionIndex < 3
            errorTransitionMatrix_dominantOrder_eachLength.(temp_str){tempi} = ...
                errorTransitionMatrix_eachLength.(temp_str){tempi}(temp_dominantOrder, :);
        else
            errorTransitionMatrix_dominantOrder_eachLength.(temp_str){tempi} = ...
                errorTransitionMatrix_eachLength.(temp_str){tempi}{temp_dominantOrder};
        end
    end
    for tempi=1:pointKindsNum
        temp_dominantOrder_negative = errorType_dominantOrder_eachLength_negative.(temp_str)(tempi);
        if isnan(temp_dominantOrder_negative)
            continue
        end
        if temp_conditionIndex < 3
            errorTransitionMatrix_dominantOrder_eachLength_negative.(temp_str){tempi} = ...
                errorTransitionMatrix_eachLength_negative.(temp_str){tempi}(temp_dominantOrder_negative, :);
        else
            errorTransitionMatrix_dominantOrder_eachLength_negative.(temp_str){tempi} = ...
                errorTransitionMatrix_eachLength_negative.(temp_str){tempi}{temp_dominantOrder_negative};
        end
    end
end
errorTransitionMatrix_dominantOrder_eachLength;

%% Figure 19
if ifFigure19 == 1
    fig19 = figure(19);
    set(gcf,'Position',[0 50 1100 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
    
    temp_color = [0.0 0.0 0.0;
        0.4 0.4 0.4;
        0.7 0.7 0.7;
        0.9 0.9 0.9];
    a = 1;
    length_rangeHead = 1;
    length_rangeTail = numFrames;
    
    for temp_conditionIndex=1:2
        nexttile
        
        if temp_conditionIndex == 1
            temp_str = 'deleteOne';
            temp_titleAppend = 'delete one';
        elseif temp_conditionIndex == 2
            temp_str = 'addOne';
            temp_titleAppend = 'add one';
        end
        
        for target_seqLength=1:pointKindsNum
            temp_rangeTail = length(errorTransitionMatrix_dominantOrder_eachLength.(temp_str){target_seqLength});
            plot(1:temp_rangeTail,errorTransitionMatrix_dominantOrder_eachLength.(temp_str){target_seqLength},...
                '-o','Color',temp_color(target_seqLength,:),'LineWidth',1.5,'MarkerSize',5);
            hold on
        end
        
        legend('length=1','length=2','length=3','length=4',...
            'Location','northwest','fontsize',10)
        
        % ylim([min(temp_mean_collaped)-0.1 1]);
        set(gca, 'FontSize', 20)
        % set(gca,'YTick',[-0.3:0.1:0.3]);%设置要显示坐标刻度的范围
        set(gca,'XLim',[length_rangeHead-1 length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XTick',[length_rangeHead:1:length_rangeTail]);%设置要显示坐标刻度的范围
        set(gca,'XTickLabel',[length_rangeHead:1:length_rangeTail]);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        xlabel('Location', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Proportion', 'FontSize', 18, 'FontWeight', 'bold');
        title({'\fontsize{20}{\bf Error type of dominant order}';temp_titleAppend});
    end
    currentFigName = ['errorTypeDeleteAddOne_dominantOrder_eachLength', '_'];
    % to generate a unique file name for saving figure
    fileName_fig19 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(19) == 1
        exportgraphics(fig19, fileName_fig19, 'BackgroundColor', 'none');
    end
end


% errorTransitionMatrix_eachLength_eachLoc_deleteOne = zeros(pointKindsNum,numFrames);
% for tempi=2:pointKindsNum
%     temp_1 = errorTransitionMatrix_eachLength.deleteOne{tempi}*numSeq(tempi);
%     temp_1(isnan(temp_1)) = 0;
%     %errorTransitionMatrix_eachLength_eachLoc_deleteOne(tempi,:) = sum(temp_1,1)./numSeq(tempi);
%     errorTransitionMatrix_eachLength_eachLoc_deleteOne(tempi,:) = sum(temp_1,1);
% end
% % errorTransitionMatrix_eachLoc_deleteOne = mean(errorTransitionMatrix_eachLength_eachLoc_deleteOne,1);
% errorTransitionMatrix_eachLoc_deleteOne = sum(errorTransitionMatrix_eachLength_eachLoc_deleteOne,1)./sum(numSeq);
%
% errorTransitionMatrix_eachLength_eachLoc_addOne = zeros(pointKindsNum,numFrames);
% for tempi=1:pointKindsNum
%     temp_1 = errorTransitionMatrix_eachLength.addOne{tempi}*numSeq(tempi);
%     temp_1(isnan(temp_1)) = 0;
%     %errorTransitionMatrix_eachLength_eachLoc_addOne(tempi,:) = sum(temp_1,1)./numSeq(tempi);
%     errorTransitionMatrix_eachLength_eachLoc_addOne(tempi,:) = sum(temp_1,1);
% end
% % errorTransitionMatrix_eachLoc_addOne = mean(errorTransitionMatrix_eachLength_eachLoc_addOne,1);
% errorTransitionMatrix_eachLoc_addOne = sum(errorTransitionMatrix_eachLength_eachLoc_addOne,1)./sum(numSeq);


errorTransitionMatrix_eachLength_raw = errorTransitionMatrix_eachLength;
errorTransitionMatrix_eachLength_new = errorTransitionMatrix_eachLength;

% errorTransitionMatrix_eachLoc_deleteOne
errorTransitionMatrix_eachLength_eachLoc_deleteOne = zeros(pointKindsNum,numFrames);
errorTransitionMatrix_eachLength_eachLoc_deleteOne_weighted = zeros(pointKindsNum,numFrames);
for tempi=2:pointKindsNum
    errorTransitionMatrix_eachLength_new.deleteOne{tempi} = ...
        errorTransitionMatrix_eachLength.deleteOne{tempi} .* ...
        errorTransitionMatrix_eachLength_validCount.deleteOne{tempi};
    errorTransitionMatrix_eachLength_new.deleteOne{tempi}(isnan(errorTransitionMatrix_eachLength_new.deleteOne{tempi})) = 0;
    errorTransitionMatrix_eachLength_new.deleteOne{tempi} = sum(errorTransitionMatrix_eachLength_new.deleteOne{tempi},1);
    errorTransitionMatrix_eachLength_new.deleteOne{tempi} = ...
        errorTransitionMatrix_eachLength_new.deleteOne{tempi} ./ sum(errorTransitionMatrix_eachLength_validCount.deleteOne{tempi},1);
    
    temp_weight = numSeq(tempi)/sum(numSeq);
    errorTransitionMatrix_eachLength_eachLoc_deleteOne(tempi,:) = errorTransitionMatrix_eachLength_new.deleteOne{tempi};
    errorTransitionMatrix_eachLength_eachLoc_deleteOne_weighted(tempi,:) = errorTransitionMatrix_eachLength_new.deleteOne{tempi}*temp_weight;
end
errorTransitionMatrix_eachLoc_deleteOne = sum(errorTransitionMatrix_eachLength_eachLoc_deleteOne_weighted,1);

% errorTransitionMatrix_eachLoc_addOne
errorTransitionMatrix_eachLength_eachLoc_addOne = zeros(pointKindsNum,numFrames);
errorTransitionMatrix_eachLength_eachLoc_addOne_weighted = zeros(pointKindsNum,numFrames);
for tempi=1:pointKindsNum
    errorTransitionMatrix_eachLength_new.addOne{tempi} = ...
        errorTransitionMatrix_eachLength.addOne{tempi} .* ...
        errorTransitionMatrix_eachLength_validCount.addOne{tempi};
    errorTransitionMatrix_eachLength_new.addOne{tempi}(isnan(errorTransitionMatrix_eachLength_new.addOne{tempi})) = 0;
    errorTransitionMatrix_eachLength_new.addOne{tempi} = sum(errorTransitionMatrix_eachLength_new.addOne{tempi},1);
    errorTransitionMatrix_eachLength_new.addOne{tempi} = ...
        errorTransitionMatrix_eachLength_new.addOne{tempi} ./ sum(errorTransitionMatrix_eachLength_validCount.addOne{tempi},1);
    
    temp_weight = numSeq(tempi)/sum(numSeq);
    errorTransitionMatrix_eachLength_eachLoc_addOne(tempi,:) = errorTransitionMatrix_eachLength_new.addOne{tempi};
    errorTransitionMatrix_eachLength_eachLoc_addOne_weighted(tempi,:) = errorTransitionMatrix_eachLength_new.addOne{tempi}*temp_weight;
end
errorTransitionMatrix_eachLoc_addOne = sum(errorTransitionMatrix_eachLength_eachLoc_addOne_weighted,1);

% errorTransitionMatrix_eachLoc_interference
errorTransitionMatrix_eachLength_eachLoc_interference = cell(1,pointKindsNum);
errorTransitionMatrix_eachLength_eachLoc_interference_weighted = cell(1,pointKindsNum);
errorTransitionMatrix_eachLoc_interference = zeros(numFrames,numFrames);
for tempi=1:pointKindsNum
    temp_interference = zeros(numFrames,numFrames);
    temp_interferenceCount = zeros(numFrames,numFrames);
    for tempj=1:tempi
        errorTransitionMatrix_eachLength_new.interference{tempi}{tempj} = ...
            errorTransitionMatrix_eachLength_new.interference{tempi}{tempj} .* ...
            errorTransitionMatrix_eachLength_validCount.interference{tempi}{tempj};
        errorTransitionMatrix_eachLength_new.interference{tempi}{tempj}(isnan(errorTransitionMatrix_eachLength_new.interference{tempi}{tempj})) = 0;
        temp_interference = temp_interference + errorTransitionMatrix_eachLength_new.interference{tempi}{tempj};
        temp_interferenceCount = temp_interferenceCount + errorTransitionMatrix_eachLength_validCount.interference{tempi}{tempj};
    end
    errorTransitionMatrix_eachLength_new.interference{tempi} = temp_interference./temp_interferenceCount;
    errorTransitionMatrix_eachLength_new.interference{tempi}(isnan(errorTransitionMatrix_eachLength_new.interference{tempi})) = 0;
    
    temp_weight = numSeq(tempi)/sum(numSeq);
    
    errorTransitionMatrix_eachLength_eachLoc_interference{tempi} = errorTransitionMatrix_eachLength_new.interference{tempi};
    errorTransitionMatrix_eachLength_eachLoc_interference_weighted{tempi} = errorTransitionMatrix_eachLength_new.interference{tempi}*temp_weight;
    errorTransitionMatrix_eachLoc_interference = errorTransitionMatrix_eachLoc_interference + errorTransitionMatrix_eachLength_eachLoc_interference_weighted{tempi};
end

errorTransitionMatrix_eachLoc = struct;
errorTransitionMatrix_eachLoc.deleteOne = errorTransitionMatrix_eachLoc_deleteOne;
errorTransitionMatrix_eachLoc.addOne = errorTransitionMatrix_eachLoc_addOne;
errorTransitionMatrix_eachLoc.interference = errorTransitionMatrix_eachLoc_interference;


%% Figure 191
if ifFigure19 == 1
    fig191 = figure(191);
    %set(gcf,'Position',[0 50 1100 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
    set(gcf,'Position',[0 50 1500 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    
    temp_color = [0.0 0.0 0.0;
        0.4 0.4 0.4;
        0.7 0.7 0.7;
        0.9 0.9 0.9];
    a = 1;
    length_rangeHead = 1;
    length_rangeTail = numFrames;
    
    for temp_conditionIndex=1:3
        nexttile
        
        if temp_conditionIndex == 1
            temp_str = 'deleteOne';
            temp_titleAppend = 'delete one';
        elseif temp_conditionIndex == 2
            temp_str = 'addOne';
            temp_titleAppend = 'add one';
        elseif temp_conditionIndex == 3
            temp_str = 'interference';
            temp_titleAppend = 'interference';
        end
        
        if temp_conditionIndex == 1 || temp_conditionIndex == 2
            
            plot(1:numFrames,errorTransitionMatrix_eachLoc.(temp_str),...
                '-o','Color',[0 0 0],'LineWidth',1.5,'MarkerSize',5);
            hold on
            
            set(gca, 'FontSize', 20)
            set(gca,'XLim',[length_rangeHead-1 length_rangeTail+1]);%X轴的数据显示范围
            set(gca,'XTick',[length_rangeHead:1:length_rangeTail]);%设置要显示坐标刻度的范围
            set(gca,'XTickLabel',[length_rangeHead:1:length_rangeTail]);%给坐标加标签
            set(gca,'box','off');% 取消右、上边框
            xlabel('Location', 'FontSize', 18, 'FontWeight', 'bold');
            ylabel('Proportion', 'FontSize', 18, 'FontWeight', 'bold');
            title(sprintf('Error type of %s',temp_titleAppend),'fontsize',20);
            
        elseif temp_conditionIndex == 3
            x = errorTransitionMatrix_eachLoc.interference;
            
            temp_max = max(x);
            temp_totalMax = max(temp_max);
            
            temp_plot = x;
            % Plot heat map
            imagesc(temp_plot, [0 temp_totalMax]);
            my_gray = gray;
            my_gray = my_gray(end:-1:1,:);
            colormap(my_gray);
            hold on
            
            
            axis equal
            xlim([0.5 numFrames+0.5]);
            ylim([0.5 numFrames+0.5]);
            
            set(gca, 'FontSize', 20)
            set(gca,'YDir','reverse');% Reverse the direction of y axis
            set(gca,'box','off');% 取消右、上边框
            set(gca,'XTick',[1:1:numFrames]);%设置要显示坐标刻度的范围
            set(gca,'YTick',[1:1:numFrames]);%设置要显示坐标刻度的范围
            
            xlabel('Response location','FontSize',18,'FontWeight','bold');
            ylabel('Stimuli location','FontSize',18,'FontWeight','bold');
            title(sprintf('Error type of %s',temp_titleAppend),'fontsize',20);
            
            c = colorbar;
            c.Label.String = 'proportion';
            
        end
    end
end

errorTransitionMatrix_dominantOrder_eachLength.interference;
if ifFigure20 == 1
    fig20 = figure(20);
    set(gcf,'Position',[0 50 1500 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,pointKindsNum,'TileSpacing','Compact','Padding','Compact');
    
    temp_max = zeros(1, pointKindsNum);
    for tempi=1:pointKindsNum
        temp_max(tempi) = max(max(errorTransitionMatrix_dominantOrder_eachLength.interference{tempi}));
    end
    temp_totalMax = max(temp_max);
    
    for tempi=1:pointKindsNum
        nexttile
        temp_plot = errorTransitionMatrix_dominantOrder_eachLength.interference{tempi};
        % Plot heat map
        imagesc(temp_plot, [0 temp_totalMax]);
        %imagesc(temp_plot);
        my_gray = gray;
        my_gray = my_gray(end:-1:1,:);
        colormap(my_gray);
        hold on
        
        
        
        axis equal
        xlim([0.5 numFrames+0.5]);
        ylim([0.5 numFrames+0.5]);
        
        set(gca, 'FontSize', 16)
        set(gca,'YDir','reverse');% Reverse the direction of y axis
        set(gca,'box','off');% 取消右、上边框
        set(gca,'XTick',[1:1:numFrames]);%设置要显示坐标刻度的范围
        set(gca,'YTick',[1:1:numFrames]);%设置要显示坐标刻度的范围
        
        xlabel('Response location','FontSize',14,'FontWeight','bold');
        ylabel('Stimuli location','FontSize',14,'FontWeight','bold');
        %title('\fontsize{16}{Inteference}');
        %title(['Inteference, length', num2str(tempi)],'FontSize',16);
        title({'\fontsize{15}{\bf Error type of dominant order}';['inteference, length', num2str(tempi)]});
        
        %         %c = colorbar('southoutside');
        %         c = colorbar;
        %         c.Label.String = 'proportion';
        %         c.Label.FontSize = 12;
    end
    %c = colorbar('southoutside');
    c = colorbar;
    c.Label.String = 'proportion';
    %c.Label.FontSize = 12;
    
    currentFigName = ['errorType_interference_dominantOrder_eachLength', '_'];
    % to generate a unique file name for saving figure
    fileName_fig20 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(20) == 1
        exportgraphics(fig20, fileName_fig20, 'BackgroundColor', 'none');
    end
end

%% Figure 202
fig202 = figure(202);
% set(gcf,'Position',[0 50 400 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set(gcf,'Position',[0 50 1100 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

errorTransitionMatrix_dominantOrder_eachLength.interference;
errorTransitionMatrix_dominantOrder_lengthAvg_interference = ...
    zeros(numSeq(1), numSeq(1));

temp_sum = errorTransitionMatrix_dominantOrder_lengthAvg_interference;
for tempi=1:pointKindsNum
    temp_data = errorTransitionMatrix_dominantOrder_eachLength.interference{tempi};
    temp_data(isnan(temp_data)) = 0;
    
    
    temp_sum = temp_sum + temp_data;
end
errorTransitionMatrix_dominantOrder_lengthAvg_interference = temp_sum / pointKindsNum;

for tempi=1:2
    nexttile;
    
    if tempi == 1
        x = errorTransitionMatrix_dominantOrder_lengthAvg_interference;
    elseif tempi == 2
        x = errorTransitionMatrix_eachLoc.interference;
    end
    
    temp_max = max(x);
    temp_totalMax = max(temp_max);
    
    temp_plot = x;
    % Plot heat map
    imagesc(temp_plot, [0 temp_totalMax]);
    my_gray = gray;
    my_gray = my_gray(end:-1:1,:);
    colormap(my_gray);
    hold on
    
    
    axis equal
    xlim([0.5 numFrames+0.5]);
    ylim([0.5 numFrames+0.5]);
    
    set(gca, 'FontSize', 16)
    set(gca,'YDir','reverse');% Reverse the direction of y axis
    set(gca,'box','off');% 取消右、上边框
    set(gca,'XTick',[1:1:numFrames]);%设置要显示坐标刻度的范围
    set(gca,'YTick',[1:1:numFrames]);%设置要显示坐标刻度的范围
    
    xlabel('Response location','FontSize',14,'FontWeight','bold');
    ylabel('Stimuli location','FontSize',14,'FontWeight','bold');
    if tempi==1
        title({'\fontsize{15}{\bf Error type of dominant order}';['inteference, length average']});
    elseif tempi==2
        title({'\fontsize{15}{\bf Error type of inteference}'});
    end
    
    c = colorbar;
    c.Label.String = 'proportion';
    
end

currentFigName = ['errorType_interference_dominantOrder_lengthAvg', '_'];
% to generate a unique file name for saving figure
fileName_fig202 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if false
    exportgraphics(fig202, fileName_fig202, 'BackgroundColor', 'none');
end


%% Correlation between offloading rate and error rate
fig21 = figure(21);

offloadingProb_inOne;
gError_noChoice_collapsed_inOne = 1 - gAcc_noChoice_collapsed_inOne;
boostMoreThan0_index = choiceMinusNoChoice_inOne'>=0;
boostLessThan0_index = choiceMinusNoChoice_inOne'<=0;
valid_boostMoreThan0_index;
valid_boostLessThan0_index;

boost_max = max(choiceMinusNoChoice_inOne(valid_index));
boost_min = min(choiceMinusNoChoice_inOne(valid_index));
boost_int = round(choiceMinusNoChoice_inOne'*100)/100;

x = gError_noChoice_collapsed_inOne(valid_boostMoreThan0_index);
y = offloadingProb_inOne(valid_boostMoreThan0_index);
n = 1;
[p_mapping_positive,S] = polyfit(x,y,n);
r2_positive = 1 - (S.normr/norm(y - mean(y)))^2;
r_positive = corr(x',y');
x_fit_positive = 0:0.1:1;
y_fit_positive = polyval(p_mapping_positive,x_fit_positive);

x = gError_noChoice_collapsed_inOne(valid_boostLessThan0_index);
y = offloadingProb_inOne(valid_boostLessThan0_index);
n = 1;
[p_mapping_negative,S] = polyfit(x,y,n);
r2_negative = 1 - (S.normr/norm(y - mean(y)))^2;
r_negative = corr(x',y');
x_fit_negative = 0:0.1:1;
y_fit_negative = polyval(p_mapping_negative,x_fit_negative);

plot(x_fit_positive, y_fit_positive, '-', 'LineWidth', 2, 'Color', [0.75 0.25 0.25 0.7]);
hold on
plot(x_fit_negative, y_fit_negative, '-', 'LineWidth', 2, 'Color', [0.25 0.25 0.75 0.7]);
hold on
temp_color = [0.75 0.25 0.25;
    0.25 0.25 0.75];

% re-order the top 2 index of temp_range_reOrder, so that legend can plot
% positive and negative correctly.
for tempi=1:sum(numSeq)
    if valid_boostMoreThan0_index(tempi) == 1
        tempIndex_positive = tempi;
        break
    end
end
for tempi=1:sum(numSeq)
    if valid_boostLessThan0_index(tempi) == 1
        tempIndex_negative = tempi;
        break
    end
end
tempIndex_positive;
tempIndex_negative;
temp_range_reOrder = 1:sum(numSeq);

temp_cache = temp_range_reOrder(1);
temp_range_reOrder(1) = tempIndex_positive;
temp_range_reOrder(tempIndex_positive) = temp_cache;

temp_cache = temp_range_reOrder(2);
temp_range_reOrder(2) = tempIndex_negative;
temp_range_reOrder(tempIndex_negative) = temp_cache;

for tempi=1:sum(numSeq)
    tempi_2 = temp_range_reOrder(tempi);
    
    if valid_boostMoreThan0_index(tempi_2) == 1
        scatter(gError_noChoice_collapsed_inOne(tempi_2), offloadingProb_inOne(tempi_2), ...
            abs(boost_int(tempi_2))*300+1, 'filled', 'MarkerFaceColor', temp_color(1, :), ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);
        hold on
        
    elseif valid_boostLessThan0_index(tempi_2) == 1
        scatter(gError_noChoice_collapsed_inOne(tempi_2), offloadingProb_inOne(tempi_2), ...
            abs(boost_int(tempi_2))*300+1, 'filled', 'MarkerFaceColor', temp_color(2, :), ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);
        hold on
    end
end

legend('Positive fit','Negative fit','Positive boost','Negative boost','Location','southeast','fontsize',8)
text(0.02,0.92,sprintf('slope(positive) = %.2f\nslope(negative) = %.2f',p_mapping_positive(1),p_mapping_negative(1)), 'fontsize',13,'FontWeight','bold');

% legend('Positive boost','Negative boost','Location','southeast','fontsize',9)
% text(0.05,0.9,sprintf('r(positive) = %.3f\nr(negative) = %.3f',r_positive,r_negative), 'fontsize',13,'FontWeight','bold');


xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Error rate', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');
temp_str1 = sprintf('Correlation, offloading rate & error rate');
temp_str2 = sprintf('Boost=[%.2f,%.2f]',boost_min,boost_max);
[t,s] = title(temp_str1,temp_str2);
t.FontSize = 18;
t.FontWeight = 'bold';
s.FontSize = 14;

currentFigName = ['correlation_rProbandError_positiveAndNegative', '_'];
% to generate a unique file name for saving figure
fileName_fig21 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(21) == 1
    exportgraphics(fig21, fileName_fig21, 'BackgroundColor', 'none');
end


%% Figure 212, Correlatioin between angle (slope) and boost
fig212 = figure(212);
set(gcf,'Position',[0 50 950 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

slope_inOne = offloadingProb_inOne./(gError_noChoice_collapsed_inOne+eps);
angle_inOne = (atan(slope_inOne)/pi)*180;

% [rho,pval] = corr(angle_inOne(valid_index)',choiceMinusNoChoice_inOne(valid_index));
% [rho,pval] = corr(slope_inOne(valid_index)',choiceMinusNoChoice_inOne(valid_index));

tempBoolIndex = ~isnan(choiceMinusNoChoice_inOne);
% tempBoolIndex = valid_index;

[rho_angle,pval_angle] = corr(angle_inOne(tempBoolIndex)',choiceMinusNoChoice_inOne(tempBoolIndex));
[rho_slope,pval_slope] = corr(slope_inOne(tempBoolIndex)',choiceMinusNoChoice_inOne(tempBoolIndex));


for temptempi=1:2
    ax = nexttile;
    
    if temptempi == 1
        x = angle_inOne(tempBoolIndex)';
    elseif temptempi == 2
        x = slope_inOne(tempBoolIndex)';
    end
    y = choiceMinusNoChoice_inOne(tempBoolIndex);
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    
    n = 1;
    [p_mapping,S] = polyfit(x,y,n);
    r2 = 1 - (S.normr/norm(y - mean(y)))^2;
    [r, p_corr] = corr(x,y);
    
    x_fit = min(x):0.1:max(x);
    y_fit = polyval(p_mapping,x_fit);
    temp_color = [0.25 0.25 0.25];
    
    scatter(x, y, ...
        25, 'filled', 'MarkerFaceColor', temp_color, ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
    
    xlim([x_min-(x_max-x_min)*0.05 x_max+(x_max-x_min)*0.05]);
    ylim([y_min-(y_max-y_min)*0.05 y_max+(y_max-y_min)*0.05]);
    
    set(gca, 'FontSize', 20)
    set(gca,'box','off');% 取消右、上边框
    ylabel('Boost effects', 'FontSize', 18, 'FontWeight', 'bold');
    if temptempi == 1
        xlabel('Angle of offloading vs error rate', 'FontSize', 18, 'FontWeight', 'bold');
        title({'\fontsize{15}Correlation, Angle & Boost';[monkey_name,', ',sprintf('r=%.3f, p=%.3f',r,p_corr)]});
    elseif temptempi == 2
        xlabel('Slope of offloading vs error rate', 'FontSize', 18, 'FontWeight', 'bold');
        title({'\fontsize{15}Correlation, Slope & Boost';[monkey_name,', ',sprintf('r=%.3f, p=%.3f',r,p_corr)]});
    end
    
end

%% Figure 22, Error types with each order, only negative boost sequences
if ifFigure22 == 1
    fig22 = figure(22);
    set(gcf,'Position',[0 50 1650 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    
    temp_color = [0.00 0.00 0.00;
        0.28 0.28 0.28;
        0.56 0.56 0.56;
        0.84 0.84 0.84];
    a = 1;
    length_rangeHead = seq_length_rangeHead;
    length_rangeTail = seq_length_rangeTail+1;
    
    for temp_conditionIndex=1:3
        nexttile
        
        if temp_conditionIndex == 1
            temp_str = 'deleteOne';
            temp_titleAppend = 'delete one, negative boost';
        elseif temp_conditionIndex == 2
            temp_str = 'addOne';
            temp_titleAppend = 'add one, negative boost';
        elseif temp_conditionIndex == 3
            temp_str = 'interference';
            temp_titleAppend = 'interference, negative boost';
        end
        
        for target_seqLength=1:pointKindsNum
            temp_rangeTail = length(errorType_eachOrder_eachLength_negative.(temp_str){target_seqLength});
            plot(1:temp_rangeTail,errorType_eachOrder_eachLength_negative.(temp_str){target_seqLength},...
                '-o','Color',temp_color(target_seqLength,:),'LineWidth',1,'MarkerSize',5);
            hold on
        end
        
        legend('length=1','length=2','length=3','length=4',...
            'Location','northeast','fontsize',9)
        
        set(gca, 'FontSize', 20)
        set(gca,'XLim',[length_rangeHead-1 length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XTick',[length_rangeHead:1:length_rangeTail]);%设置要显示坐标刻度的范围
        set(gca,'XTickLabel',[length_rangeHead:1:length_rangeTail]);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        xlabel('Order', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Proportion', 'FontSize', 18, 'FontWeight', 'bold');
        title({'\fontsize{20}{\bf Error types with each order}';temp_titleAppend});
    end
    currentFigName = ['errorTypesWithEachOrder_negative', '_'];
    % to generate a unique file name for saving figure
    fileName_fig22 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(22) == 1
        exportgraphics(fig22, fileName_fig22, 'BackgroundColor', 'none');
    end
end

%% Figure 23, Error type of dominant order, negative boost sequences
if ifFigure23 == 1
    fig23 = figure(23);
    set(gcf,'Position',[0 50 1100 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
    
    temp_color = [0.0 0.0 0.0;
        0.4 0.4 0.4;
        0.7 0.7 0.7;
        0.9 0.9 0.9];
    a = 1;
    length_rangeHead = 1;
    length_rangeTail = numFrames;
    
    for temp_conditionIndex=1:2
        nexttile
        
        if temp_conditionIndex == 1
            temp_str = 'deleteOne';
            temp_titleAppend = 'delete one, negative boost';
        elseif temp_conditionIndex == 2
            temp_str = 'addOne';
            temp_titleAppend = 'add one, negative boost';
        end
        
        for target_seqLength=1:pointKindsNum
            temp_rangeTail = length(errorTransitionMatrix_dominantOrder_eachLength_negative.(temp_str){target_seqLength});
            plot(1:temp_rangeTail,errorTransitionMatrix_dominantOrder_eachLength_negative.(temp_str){target_seqLength},...
                '-o','Color',temp_color(target_seqLength,:),'LineWidth',1.5,'MarkerSize',5);
            hold on
        end
        
        legend('length=1','length=2','length=3','length=4',...
            'Location','northwest','fontsize',10)
        
        set(gca, 'FontSize', 20)
        set(gca,'XLim',[length_rangeHead-1 length_rangeTail+1]);%X轴的数据显示范围
        set(gca,'XTick',[length_rangeHead:1:length_rangeTail]);%设置要显示坐标刻度的范围
        set(gca,'XTickLabel',[length_rangeHead:1:length_rangeTail]);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        xlabel('Location', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Proportion', 'FontSize', 18, 'FontWeight', 'bold');
        title({'\fontsize{20}{\bf Error type of dominant order}';temp_titleAppend});
    end
    currentFigName = ['errorTypeDeleteAddOne_dominantOrder_eachLength_negative', '_'];
    % to generate a unique file name for saving figure
    fileName_fig23 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(23) == 1
        exportgraphics(fig23, fileName_fig23, 'BackgroundColor', 'none');
    end
end

%% Figure 24, Error type of dominant order, negative boost sequences
if ifFigure24 == 1
    fig24 = figure(24);
    set(gcf,'Position',[0 50 1500 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,pointKindsNum,'TileSpacing','Compact','Padding','Compact');
    
    temp_max = zeros(1, pointKindsNum);
    for tempi=1:pointKindsNum
        temp_max(tempi) = max(max(errorTransitionMatrix_dominantOrder_eachLength_negative.interference{tempi}));
    end
    temp_totalMax = max(temp_max);
    
    for tempi=1:pointKindsNum
        nexttile
        temp_plot = errorTransitionMatrix_dominantOrder_eachLength_negative.interference{tempi};
        % Plot heat map
        imagesc(temp_plot, [0 temp_totalMax]);
        my_gray = gray;
        my_gray = my_gray(end:-1:1,:);
        colormap(my_gray);
        hold on
        
        
        
        axis equal
        xlim([0.5 numFrames+0.5]);
        ylim([0.5 numFrames+0.5]);
        
        set(gca, 'FontSize', 16)
        set(gca,'YDir','reverse');% Reverse the direction of y axis
        set(gca,'box','off');% 取消右、上边框
        set(gca,'XTick',[1:1:numFrames]);%设置要显示坐标刻度的范围
        set(gca,'YTick',[1:1:numFrames]);%设置要显示坐标刻度的范围
        
        xlabel('Response location','FontSize',14,'FontWeight','bold');
        ylabel('Stimuli location','FontSize',14,'FontWeight','bold');
        title({'\fontsize{15}{\bf Error type of dominant order}';['inteference, negative boost, length', num2str(tempi)]});
        
    end
    c = colorbar;
    c.Label.String = 'proportion';
    
    currentFigName = ['errorType_interference_dominantOrder_eachLength_negative', '_'];
    % to generate a unique file name for saving figure
    fileName_fig24 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(24) == 1
        exportgraphics(fig24, fileName_fig24, 'BackgroundColor', 'none');
    end
end

%% Compute summarized accuracy in each location of positive and negative seqs (GLM-like)
gAcc_noChoice_collapsed_inOne;
valid_boostMoreThan0_index;
valid_boostLessThan0_index;
target_seqSet_inOne;
location_positive_count = zeros(1, numFrames);
location_positive_sumAcc = zeros(1, numFrames);
location_negative_count = zeros(1, numFrames);
location_negative_sumAcc = zeros(1, numFrames);
for tempi=1:sum(numSeq)
    if valid_boostMoreThan0_index(tempi) == 0
        continue
    end
    target_seqSet_inOne{tempi};
    location_positive_count(target_seqSet_inOne{tempi}) = ...
        location_positive_count(target_seqSet_inOne{tempi}) + 1;
    location_positive_sumAcc(target_seqSet_inOne{tempi}) = ...
        location_positive_sumAcc(target_seqSet_inOne{tempi}) + gAcc_noChoice_collapsed_inOne(tempi);
end
for tempi=1:sum(numSeq)
    if valid_boostLessThan0_index(tempi) == 0
        continue
    end
    target_seqSet_inOne{tempi};
    location_negative_count(target_seqSet_inOne{tempi}) = ...
        location_negative_count(target_seqSet_inOne{tempi}) + 1;
    location_negative_sumAcc(target_seqSet_inOne{tempi}) = ...
        location_negative_sumAcc(target_seqSet_inOne{tempi}) + gAcc_noChoice_collapsed_inOne(tempi);
end
location_positive_sumAcc_n11 = location_positive_sumAcc ./ location_positive_count;
location_negative_sumAcc_n11 = location_negative_sumAcc ./ location_negative_count;


%
gAcc_noChoice_collapsed;
valid_boostMoreThan0_index;
valid_boostLessThan0_index;
valid_boostMoreThan0_index_inLength = cell(1, pointKindsNum);
valid_boostLessThan0_index_inLength = cell(1, pointKindsNum);
for tempi=1:pointKindsNum
    if tempi == 1
        temp_range = 1:numSeq(1);
    else
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    end
    temp_range;
    valid_boostMoreThan0_index_inLength{tempi} = valid_boostMoreThan0_index(temp_range);
    valid_boostLessThan0_index_inLength{tempi} = valid_boostLessThan0_index(temp_range);
    
end

%% Figure25
fig25 = figure(25);
set(gcf,'Position',[0 50 1650 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout('flow','TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

backgroundColor = [1 1 1];

p_noChoice_PONE = zeros(1, pointKindsNum);
p_choice_PONE = zeros(1, pointKindsNum);
for tempi=2:pointKindsNum
    ax = nexttile;
    set(gca,'color',backgroundColor);
    
    temp_gAcc_choice_positive = gAcc_choice_collapsed{tempi}(valid_boostMoreThan0_index_inLength{tempi})';
    temp_gAcc_choice_negative = gAcc_choice_collapsed{tempi}(valid_boostLessThan0_index_inLength{tempi})';
    [h_choice_PONE,p_choice_PONE(tempi)]=ttest2(temp_gAcc_choice_positive, temp_gAcc_choice_negative);
    
    temp_gAcc_noChoice_positive = gAcc_noChoice_collapsed{tempi}(valid_boostMoreThan0_index_inLength{tempi})';
    temp_gAcc_noChoice_negative = gAcc_noChoice_collapsed{tempi}(valid_boostLessThan0_index_inLength{tempi})';
    [h_noChoice_PONE,p_noChoice_PONE(tempi)]=ttest2(temp_gAcc_noChoice_positive, temp_gAcc_noChoice_negative);
    
    
    temp_mean_choice = [mean(temp_gAcc_choice_positive) mean(temp_gAcc_choice_negative)];
    temp_mean_noChoice = [mean(temp_gAcc_noChoice_positive) mean(temp_gAcc_noChoice_negative)];
    y_bar = [temp_mean_choice; temp_mean_noChoice]';
    
    bar_width = 0.8;
    b = bar([0 1],y_bar,bar_width,'grouped',...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    b(1).FaceColor = [0.4660 0.6740 0.1880];
    b(2).FaceColor = [0.45 0.45 0.45];
    
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(y_bar);
    % Get the x coordinate of the bars
    x_bar = nan(ngroups, nbars);
    for i = 1:nbars
        x_bar(:,i) = b(i).XEndPoints;
    end
    
    temp_SEM_choice = [std(temp_gAcc_choice_positive)./sqrt(length(temp_gAcc_choice_positive)) ...
        std(temp_gAcc_choice_negative)./sqrt(length(temp_gAcc_choice_negative))];
    temp_SEM_noChoice = [std(temp_gAcc_noChoice_positive)./sqrt(length(temp_gAcc_noChoice_positive)) ...
        std(temp_gAcc_noChoice_negative)./sqrt(length(temp_gAcc_noChoice_negative))];
    errorbar(x_bar, y_bar,[temp_SEM_choice; temp_SEM_noChoice]', '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
    hold on
    
    tempTxt_choice = sprintf('');
    if p_choice_PONE(tempi) < 0.001/(pointKindsNum-1)
        tempTxt_choice = sprintf('***');
    elseif p_choice_PONE(tempi) < 0.01/(pointKindsNum-1)
        tempTxt_choice = sprintf('**');
    elseif p_choice_PONE(tempi) < 0.05/(pointKindsNum-1)
        tempTxt_choice = sprintf('*');
    end
    
    tempTxt_noChoice = sprintf('');
    if p_noChoice_PONE(tempi) < 0.001/(pointKindsNum-1)
        tempTxt_noChoice = sprintf('***');
    elseif p_noChoice_PONE(tempi) < 0.01/(pointKindsNum-1)
        tempTxt_noChoice = sprintf('**');
    elseif p_noChoice_PONE(tempi) < 0.05/(pointKindsNum-1)
        tempTxt_noChoice = sprintf('*');
    end
    
    x_bar;
    y_bar;
    plot([x_bar(1,1)+0.6*bar_width/4 x_bar(2,1)],[y_bar(1,1) y_bar(1,1)],'Color',b(1).FaceColor,'LineWidth',2);
    hold on
    plot([x_bar(1,2)+0.6*bar_width/4 x_bar(2,2)],[y_bar(1,2) y_bar(1,2)],'Color',b(2).FaceColor,'LineWidth',2);
    hold on
    
    tempTxt_choice;
    tempTxt_noChoice;
    text(0.5,y_bar(1,1)+0.02,tempTxt_choice, 'fontsize',15,'FontWeight','bold');
    text(0.6,y_bar(1,2)+0.02,tempTxt_noChoice, 'fontsize',15,'FontWeight','bold');
    
    legend('choice','no choice',...
        'Location','northeast','fontsize',10)
    
    temp_category1 = 'Positive';
    temp_category2 = 'Negative';
    set(gca,'XTick',0:1);
    set(gca,'XTickLabel', [temp_category1; temp_category2]);%给坐标加标签
    ylim([0 1]);
    set(gca, 'FontSize', 16)
    set(gca,'box','off');% 取消右、上边框
    ylabel('Accuracy', 'FontSize', 16, 'FontWeight', 'bold');
    
    
    
    temp_str1 = ['Length=',num2str(tempi),',monkey ',monkey_name];
    temp_str2 = sprintf('p(choice)=%.3f, p(noChoice)=%.3f',p_choice_PONE(tempi),p_noChoice_PONE(tempi));
    [t,s] = title(temp_str1,temp_str2);
    t.FontSize = 14;
    t.FontWeight = 'bold';
    s.FontSize = 12;
    
    currentFigName = ['boostAccCompare_positiveAndNegative', '_'];
    % to generate a unique file name for saving figure
    fileName_fig25 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(25) == 1
        exportgraphics(fig25, fileName_fig25, 'BackgroundColor', 'none');
    end
end

%% Figure252
fig252 = figure(252);
set(gcf,'Position',[0 50 1650 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout('flow','TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

backgroundColor = [1 1 1];

p_noChoice_PONE = zeros(1, pointKindsNum);
p_choice_PONE = zeros(1, pointKindsNum);
for tempi=2:pointKindsNum
    ax = nexttile;
    set(gca,'color',backgroundColor);
    
    temp_gAcc_choice_positive = gAcc_choice_collapsed{tempi}(valid_boostMoreThan0_index_inLength{tempi})';
    temp_gAcc_choice_negative = gAcc_choice_collapsed{tempi}(valid_boostLessThan0_index_inLength{tempi})';
    [h_choice_PONE,p_choice_PONE(tempi)]=ttest2(temp_gAcc_choice_positive, temp_gAcc_choice_negative);
    
    temp_rProb_positive = offloadingProb_collapsed{tempi}(valid_boostMoreThan0_index_inLength{tempi})';
    temp_rProb_negative = offloadingProb_collapsed{tempi}(valid_boostLessThan0_index_inLength{tempi})';
    [h_rProb_PONE,p_rProb_PONE(tempi)]=ttest2(temp_rProb_positive, temp_rProb_negative);
    
    
    temp_mean_choice = [mean(temp_gAcc_choice_positive) mean(temp_gAcc_choice_negative)];
    temp_mean_rProb = [mean(temp_rProb_positive) mean(temp_rProb_negative)];
    y_bar = [temp_mean_choice; temp_mean_rProb]';
    
    bar_width = 0.8;
    b = bar([0 1],y_bar,bar_width,'grouped',...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    b(1).FaceColor = [0.4660 0.6740 0.1880];
    b(2).FaceColor = [0.6350 0.0780 0.1840];
    
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(y_bar);
    % Get the x coordinate of the bars
    x_bar = nan(ngroups, nbars);
    for i = 1:nbars
        x_bar(:,i) = b(i).XEndPoints;
    end
    
    temp_SEM_choice = [std(temp_gAcc_choice_positive)./sqrt(length(temp_gAcc_choice_positive)) ...
        std(temp_gAcc_choice_negative)./sqrt(length(temp_gAcc_choice_negative))];
    temp_SEM_rProb = [std(temp_rProb_positive)./sqrt(length(temp_rProb_positive)) ...
        std(temp_rProb_negative)./sqrt(length(temp_rProb_negative))];
    errorbar(x_bar, y_bar,[temp_SEM_choice; temp_SEM_rProb]', '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
    hold on
    
    tempTxt_choice = sprintf('');
    if p_choice_PONE(tempi) < 0.001/(pointKindsNum-1)
        tempTxt_choice = sprintf('***');
    elseif p_choice_PONE(tempi) < 0.01/(pointKindsNum-1)
        tempTxt_choice = sprintf('**');
    elseif p_choice_PONE(tempi) < 0.05/(pointKindsNum-1)
        tempTxt_choice = sprintf('*');
    end
    
    tempTxt_rProb = sprintf('');
    if p_rProb_PONE(tempi) < 0.001/(pointKindsNum-1)
        tempTxt_rProb = sprintf('***');
    elseif p_rProb_PONE(tempi) < 0.01/(pointKindsNum-1)
        tempTxt_rProb = sprintf('**');
    elseif p_rProb_PONE(tempi) < 0.05/(pointKindsNum-1)
        tempTxt_rProb = sprintf('*');
    end
    
    x_bar;
    y_bar;
    if strcmp(tempTxt_choice,sprintf('')) == 0
        plot([x_bar(1,1)+0.6*bar_width/4 x_bar(2,1)],[y_bar(1,1) y_bar(1,1)],'Color',b(1).FaceColor,'LineWidth',2);
        hold on
    end
    if strcmp(tempTxt_rProb,sprintf('')) == 0
        plot([x_bar(1,2)+0.6*bar_width/4 x_bar(2,2)],[y_bar(2,2) y_bar(2,2)],'Color',b(2).FaceColor,'LineWidth',2);
        hold on
    end
    
    tempTxt_choice;
    tempTxt_noChoice;
    %     text(0.5,y_bar(1,1)+0.02,tempTxt_choice, 'fontsize',15,'FontWeight','bold');
    %     text(0.6,y_bar(2,2)+0.02,tempTxt_rProb, 'fontsize',15,'FontWeight','bold');
    text(0.2,y_bar(1,1)-0.04,tempTxt_choice, 'fontsize',15,'FontWeight','bold');
    text(0.8,y_bar(2,2)+0.02,tempTxt_rProb, 'fontsize',15,'FontWeight','bold');
    
    legend('internal memory','offloading',...
        'Location','northwest','fontsize',10)
    
    temp_category1 = 'Positive';
    temp_category2 = 'Negative';
    set(gca,'XTick',0:1);
    set(gca,'XTickLabel', [temp_category1; temp_category2]);%给坐标加标签
    ylim([0 1]);
    set(gca, 'FontSize', 16)
    set(gca,'box','off');% 取消右、上边框
    ylabel('Accuracy', 'FontSize', 16, 'FontWeight', 'bold');
    
    
    
    temp_str1 = ['Length=',num2str(tempi),',monkey ',monkey_name];
    temp_str2 = sprintf('p(memory)=%.3f, p(offload)=%.3f',p_choice_PONE(tempi),p_rProb_PONE(tempi));
    [t,s] = title(temp_str1,temp_str2);
    t.FontSize = 14;
    t.FontWeight = 'bold';
    s.FontSize = 12;
    
    currentFigName = ['boostAccRProbCompare_positiveAndNegative', '_'];
    % to generate a unique file name for saving figure
    fileName_fig252 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if false
        exportgraphics(fig252, fileName_fig252, 'BackgroundColor', 'none');
    end
end

%% Figure253
fig253 = figure(253);
set(gcf,'Position',[0 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
backgroundColor = [1 1 1];

p_choice = 0;
p_rProb = 0;

set(gca,'color',backgroundColor);

temp_gAcc_choice_positive = gAcc_choice_collapsed_inOne(valid_boostMoreThan0_index);
temp_gAcc_choice_negative = gAcc_choice_collapsed_inOne(valid_boostLessThan0_index);
[h_choice,p_choice]=ttest2(temp_gAcc_choice_positive, temp_gAcc_choice_negative);

temp_rProb_positive = offloadingProb_inOne(valid_boostMoreThan0_index);
temp_rProb_negative = offloadingProb_inOne(valid_boostLessThan0_index);
[h_rProb,p_rProb]=ttest2(temp_rProb_positive, temp_rProb_negative);


temp_mean_choice = [mean(temp_gAcc_choice_positive) mean(temp_gAcc_choice_negative)];
temp_mean_rProb = [mean(temp_rProb_positive) mean(temp_rProb_negative)];
y_bar = [temp_mean_choice; temp_mean_rProb]';

bar_width = 0.8;
b = bar([0 1],y_bar,bar_width,'grouped',...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
hold on
b(1).FaceColor = [0.4660 0.6740 0.1880];
b(2).FaceColor = [0.6350 0.0780 0.1840];

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(y_bar);
% Get the x coordinate of the bars
x_bar = nan(ngroups, nbars);
for i = 1:nbars
    x_bar(:,i) = b(i).XEndPoints;
end

temp_SEM_choice = [std(temp_gAcc_choice_positive)./sqrt(length(temp_gAcc_choice_positive)) ...
    std(temp_gAcc_choice_negative)./sqrt(length(temp_gAcc_choice_negative))];
temp_SEM_rProb = [std(temp_rProb_positive)./sqrt(length(temp_rProb_positive)) ...
    std(temp_rProb_negative)./sqrt(length(temp_rProb_negative))];
errorbar(x_bar, y_bar,[temp_SEM_choice; temp_SEM_rProb]', '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
hold on

tempTxt_choice = sprintf('');
if p_choice < 0.001
    tempTxt_choice = sprintf('***');
elseif p_choice < 0.01
    tempTxt_choice = sprintf('**');
elseif p_choice < 0.05
    tempTxt_choice = sprintf('*');
end

tempTxt_rProb = sprintf('');
if p_rProb < 0.001
    tempTxt_rProb = sprintf('***');
elseif p_rProb < 0.01
    tempTxt_rProb = sprintf('**');
elseif p_rProb < 0.05
    tempTxt_rProb = sprintf('*');
end

x_bar;
y_bar;
if strcmp(tempTxt_choice,sprintf('')) == 0
    plot([x_bar(1,1)+0.6*bar_width/4 x_bar(2,1)],[y_bar(1,1) y_bar(1,1)],'Color',b(1).FaceColor,'LineWidth',2);
    hold on
end
if strcmp(tempTxt_rProb,sprintf('')) == 0
    plot([x_bar(1,2)+0.6*bar_width/4 x_bar(2,2)],[y_bar(2,2) y_bar(2,2)],'Color',b(2).FaceColor,'LineWidth',2);
    hold on
end

tempTxt_choice;
tempTxt_noChoice;
%     text(0.5,y_bar(1,1)+0.02,tempTxt_choice, 'fontsize',15,'FontWeight','bold');
%     text(0.6,y_bar(2,2)+0.02,tempTxt_rProb, 'fontsize',15,'FontWeight','bold');
text(0.2,y_bar(1,1)-0.04,tempTxt_choice, 'fontsize',15,'FontWeight','bold');
text(0.8,y_bar(2,2)+0.02,tempTxt_rProb, 'fontsize',15,'FontWeight','bold');

legend('internal memory','offloading',...
    'Location','northwest','fontsize',10)

temp_category1 = 'Positive';
temp_category2 = 'Negative';
set(gca,'XTick',0:1);
set(gca,'XTickLabel', [temp_category1; temp_category2]);%给坐标加标签
ylim([0 1]);
set(gca, 'FontSize', 16)
set(gca,'box','off');% 取消右、上边框
ylabel('Accuracy', 'FontSize', 16, 'FontWeight', 'bold');



temp_str1 = ['Monkey ',monkey_name];
temp_str2 = sprintf('p(memory)=%.3f, p(offload)=%.3f',p_choice,p_rProb);
[t,s] = title(temp_str1,temp_str2);
t.FontSize = 14;
t.FontWeight = 'bold';
s.FontSize = 12;

currentFigName = ['boostAccRProbCompare_positiveAndNegative_lengthMerge', '_'];
% to generate a unique file name for saving figure
fileName_fig253 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if false
    exportgraphics(fig253, fileName_fig253, 'BackgroundColor', 'none');
end


%% Figure254
fig254 = figure(254);
set(gcf,'Position',[0 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
backgroundColor = [1 1 1];

p_choice = 0;
p_noChoice = 0;
p_positive = 0;
p_negative = 0;

set(gca,'color',backgroundColor);

temp_gAcc_choice_positive = gAcc_choice_collapsed_inOne(valid_boostMoreThan0_index);
temp_gAcc_choice_negative = gAcc_choice_collapsed_inOne(valid_boostLessThan0_index);
[h_choice,p_choice]=ttest2(temp_gAcc_choice_positive, temp_gAcc_choice_negative);

temp_gAcc_noChoice_positive = gAcc_noChoice_collapsed_inOne(valid_boostMoreThan0_index);
temp_gAcc_noChoice_negative = gAcc_noChoice_collapsed_inOne(valid_boostLessThan0_index);
[h_noChoice,p_noChoice]=ttest2(temp_gAcc_noChoice_positive, temp_gAcc_noChoice_negative);

[h_positive,p_positive]=ttest(temp_gAcc_choice_positive, temp_gAcc_noChoice_positive);
[h_negative,p_negative]=ttest(temp_gAcc_choice_negative, temp_gAcc_noChoice_negative);


temp_mean_choice = [mean(temp_gAcc_choice_positive) mean(temp_gAcc_choice_negative)];
temp_mean_noChoice = [mean(temp_gAcc_noChoice_positive) mean(temp_gAcc_noChoice_negative)];
y_bar = [temp_mean_choice; temp_mean_noChoice]';

bar_width = 0.8;
b = bar([0 1],y_bar,bar_width,'grouped',...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
hold on
b(1).FaceColor = [0.4660 0.6740 0.1880];
b(2).FaceColor = [0.45 0.45 0.45];

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(y_bar);
% Get the x coordinate of the bars
x_bar = nan(ngroups, nbars);
for i = 1:nbars
    x_bar(:,i) = b(i).XEndPoints;
end

temp_SEM_choice = [std(temp_gAcc_choice_positive)./sqrt(length(temp_gAcc_choice_positive)) ...
    std(temp_gAcc_choice_negative)./sqrt(length(temp_gAcc_choice_negative))];
temp_SEM_noChoice = [std(temp_gAcc_noChoice_positive)./sqrt(length(temp_gAcc_noChoice_positive)) ...
    std(temp_gAcc_noChoice_negative)./sqrt(length(temp_gAcc_noChoice_negative))];
errorbar(x_bar, y_bar,[temp_SEM_choice; temp_SEM_noChoice]', '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
hold on

tempTxt_choice = sprintf('');
if p_choice < 0.001
    tempTxt_choice = sprintf('***');
elseif p_choice < 0.01
    tempTxt_choice = sprintf('**');
elseif p_choice < 0.05
    tempTxt_choice = sprintf('*');
end

tempTxt_noChoice = sprintf('');
if p_noChoice < 0.001
    tempTxt_noChoice = sprintf('***');
elseif p_noChoice < 0.01
    tempTxt_noChoice = sprintf('**');
elseif p_noChoice < 0.05
    tempTxt_noChoice = sprintf('*');
end

tempTxt_positive = sprintf('');
if p_positive < 0.001
    tempTxt_positive = sprintf('***');
elseif p_positive < 0.01
    tempTxt_positive = sprintf('**');
elseif p_positive < 0.05
    tempTxt_positive = sprintf('*');
end

tempTxt_negative = sprintf('');
if p_negative < 0.001
    tempTxt_negative = sprintf('***');
elseif p_negative < 0.01
    tempTxt_negative = sprintf('**');
elseif p_negative < 0.05
    tempTxt_negative = sprintf('*');
end

x_bar;
y_bar;
if strcmp(tempTxt_choice,sprintf('')) == 0
    plot([x_bar(1,1)+0.6*bar_width/4 x_bar(2,1)],[y_bar(1,1) y_bar(1,1)],'Color',b(1).FaceColor,'LineWidth',2);
    hold on
end
if strcmp(tempTxt_noChoice,sprintf('')) == 0
    plot([x_bar(1,2)+0.6*bar_width/4 x_bar(2,2)],[y_bar(1,2) y_bar(1,2)],'Color',b(2).FaceColor,'LineWidth',2);
    hold on
end


tempTxt_choice;
tempTxt_noChoice;
text(0.5,y_bar(1,1)+0.02,tempTxt_choice, 'fontsize',15,'FontWeight','bold');
text(0.6,y_bar(1,2)+0.02,tempTxt_noChoice, 'fontsize',15,'FontWeight','bold');
% text(0.2,y_bar(1,1)-0.04,tempTxt_choice, 'fontsize',15,'FontWeight','bold');
% text(0.8,y_bar(2,2)+0.02,tempTxt_noChoice, 'fontsize',15,'FontWeight','bold');

tempTxt_positive;
tempTxt_negative;
text(0-0.05,max(y_bar(1,:))+0.07,tempTxt_positive, 'fontsize',15,'FontWeight','bold');
text(1-0.05,max(y_bar(2,:))+0.07,tempTxt_negative, 'fontsize',15,'FontWeight','bold');

legend('choice memory','no choice memory',...
    'Location','northeast','fontsize',10)

temp_category1 = 'Positive';
temp_category2 = 'Negative';
set(gca,'XTick',0:1);
set(gca,'XTickLabel', [temp_category1; temp_category2]);%给坐标加标签
ylim([0 1]);
set(gca, 'FontSize', 16)
set(gca,'box','off');% 取消右、上边框
ylabel('Accuracy', 'FontSize', 16, 'FontWeight', 'bold');



temp_str1 = ['Monkey ',monkey_name];
temp_str2 = sprintf('p(choice)=%.3f, p(noChoice)=%.3f',p_choice,p_noChoice);
[t,s] = title(temp_str1,temp_str2);
t.FontSize = 14;
t.FontWeight = 'bold';
s.FontSize = 12;

currentFigName = ['boostAccCompare_positiveAndNegative_lengthMerge', '_'];
% to generate a unique file name for saving figure
fileName_fig254 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if false
    exportgraphics(fig254, fileName_fig254, 'BackgroundColor', 'none');
end


%% Correlation between error rate (choice) and error rate (noChoice)
fig26 = figure(26);

gAcc_choice_collapsed_inOne;
gAcc_noChoice_collapsed_inOne;
valid_boostMoreThan0_index;
valid_boostLessThan0_index;

boost_max = max(choiceMinusNoChoice_inOne(valid_index));
boost_min = min(choiceMinusNoChoice_inOne(valid_index));
boost_int = round(choiceMinusNoChoice_inOne'*100)/100;

x = gAcc_noChoice_collapsed_inOne(valid_boostMoreThan0_index);
y = gAcc_choice_collapsed_inOne(valid_boostMoreThan0_index);
n = 1;
[p_mapping_positive,S] = polyfit(x,y,n);
r2_positive = 1 - (S.normr/norm(y - mean(y)))^2;
r_positive = corr(x',y');
x_fit_positive = 0:0.1:1;
y_fit_positive = polyval(p_mapping_positive,x_fit_positive);

x = gAcc_noChoice_collapsed_inOne(valid_boostLessThan0_index);
y = gAcc_choice_collapsed_inOne(valid_boostLessThan0_index);
n = 1;
[p_mapping_negative,S] = polyfit(x,y,n);
r2_negative = 1 - (S.normr/norm(y - mean(y)))^2;
r_negative = corr(x',y');
x_fit_negative = 0:0.1:1;
y_fit_negative = polyval(p_mapping_negative,x_fit_negative);

plot([0 1], [0 1], '-', 'LineWidth', 2, 'Color', [0.45 0.45 0.45 0.6]);
hold on

% plot(x_fit_positive, y_fit_positive, '-', 'LineWidth', 2, 'Color', [0.75 0.25 0.25 0.7]);
% hold on
% plot(x_fit_negative, y_fit_negative, '-', 'LineWidth', 2, 'Color', [0.25 0.25 0.75 0.7]);
% hold on
temp_color = [0.75 0.25 0.25;
    0.25 0.25 0.75];

% re-order the top 2 index of temp_range_reOrder, so that legend can plot
% positive and negative correctly.
for tempi=1:sum(numSeq)
    if valid_boostMoreThan0_index(tempi) == 1
        tempIndex_positive = tempi;
        break
    end
end
for tempi=1:sum(numSeq)
    if valid_boostLessThan0_index(tempi) == 1
        tempIndex_negative = tempi;
        break
    end
end
tempIndex_positive;
tempIndex_negative;
temp_range_reOrder = 1:sum(numSeq);

temp_cache = temp_range_reOrder(1);
temp_range_reOrder(1) = tempIndex_positive;
temp_range_reOrder(tempIndex_positive) = temp_cache;

temp_cache = temp_range_reOrder(2);
temp_range_reOrder(2) = tempIndex_negative;
temp_range_reOrder(tempIndex_negative) = temp_cache;

for tempi=1:sum(numSeq)
    tempi_2 = temp_range_reOrder(tempi);
    
    if valid_boostMoreThan0_index(tempi_2) == 1
        scatter(gAcc_noChoice_collapsed_inOne(tempi_2), gAcc_choice_collapsed_inOne(tempi_2), ...
            abs(boost_int(tempi_2))*300+1, 'filled', 'MarkerFaceColor', temp_color(1, :), ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);
        hold on
        
    elseif valid_boostLessThan0_index(tempi_2) == 1
        scatter(gAcc_noChoice_collapsed_inOne(tempi_2), gAcc_choice_collapsed_inOne(tempi_2), ...
            abs(boost_int(tempi_2))*300+1, 'filled', 'MarkerFaceColor', temp_color(2, :), ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);
        hold on
    end
end
legend('Diagonal','Positive boost','Negative boost','Location','southeast','fontsize',8)
% legend('Positive fit','Negative fit','Positive boost','Negative boost','Location','southeast','fontsize',8)
% text(0.02,0.92,sprintf('slope(positive) = %.2f\nslope(negative) = %.2f',p_mapping_positive(1),p_mapping_negative(1)), 'fontsize',13,'FontWeight','bold');


xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Accuracy (noChoice)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Accuracy (choice)', 'FontSize', 18, 'FontWeight', 'bold');
temp_str1 = sprintf('Correlation of accuracy, noChoice & choice');
temp_str2 = sprintf('Boost=[%.2f,%.2f]',boost_min,boost_max);
[t,s] = title(temp_str1,temp_str2);
t.FontSize = 18;
t.FontWeight = 'bold';
s.FontSize = 14;

currentFigName = ['correlation_accuracy_choiceAndNoChoice_positiveAndNegative', '_'];
% to generate a unique file name for saving figure
fileName_fig26 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(26) == 1
    exportgraphics(fig26, fileName_fig26, 'BackgroundColor', 'none');
end

if false
    %temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\sequenceInfo_Z';
    temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\sequenceInfo_D'; %#ok<*UNRCH>
    save(temp_str, 'valid_boostMoreThan0_index', 'valid_boostLessThan0_index', ...
        'target_seqSet_inOne', 'valid_index');
end

%% Correlation between error rate (choice) and error rate (noChoice)
fig262 = figure(262);

gAcc_choice_collapsed_inOne;
gAcc_noChoice_collapsed_inOne;
valid_boostMoreThan0_index;
valid_boostLessThan0_index;

boost_max = max(choiceMinusNoChoice_inOne(valid_index));
boost_min = min(choiceMinusNoChoice_inOne(valid_index));
boost_int = round(choiceMinusNoChoice_inOne'*100)/100;

x = gAcc_noChoice_collapsed_inOne(valid_boostMoreThan0_index);
y = gAcc_choice_collapsed_inOne(valid_boostMoreThan0_index);
n = 1;
[p_mapping_positive,S] = polyfit(x,y,n);
r2_positive = 1 - (S.normr/norm(y - mean(y)))^2;
r_positive = corr(x',y');
x_fit_positive = 0:0.1:1;
y_fit_positive = polyval(p_mapping_positive,x_fit_positive);

x = gAcc_noChoice_collapsed_inOne(valid_boostLessThan0_index);
y = gAcc_choice_collapsed_inOne(valid_boostLessThan0_index);
n = 1;
[p_mapping_negative,S] = polyfit(x,y,n);
r2_negative = 1 - (S.normr/norm(y - mean(y)))^2;
r_negative = corr(x',y');
x_fit_negative = 0:0.1:1;
y_fit_negative = polyval(p_mapping_negative,x_fit_negative);

plot([0 1], [0 1], '-', 'LineWidth', 2, 'Color', [0.45 0.45 0.45 0.6]);
hold on

temp_color = [0.75 0.25 0.25;
    0.25 0.25 0.75];

% re-order the top 2 index of temp_range_reOrder, so that legend can plot
% positive and negative correctly.
for tempi=1:sum(numSeq)
    if valid_boostMoreThan0_index(tempi) == 1
        tempIndex_positive = tempi;
        break
    end
end
for tempi=1:sum(numSeq)
    if valid_boostLessThan0_index(tempi) == 1
        tempIndex_negative = tempi;
        break
    end
end
tempIndex_positive;
tempIndex_negative;
temp_range_reOrder = 1:sum(numSeq);

temp_cache = temp_range_reOrder(1);
temp_range_reOrder(1) = tempIndex_positive;
temp_range_reOrder(tempIndex_positive) = temp_cache;

temp_cache = temp_range_reOrder(2);
temp_range_reOrder(2) = tempIndex_negative;
temp_range_reOrder(tempIndex_negative) = temp_cache;

temp_length = zeros(1, sum(numSeq));
for tempi=1:pointKindsNum
    if tempi == 1
        temp_range = 1:numSeq(1);
    else
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    end
    temp_range;
    temp_length(temp_range) = tempi;
end
% temp_size = temp_length*15 - 10;
% temp_size = temp_length*30 - 20;
% temp_size = (temp_length.^2)*10 - 5;
temp_size = (temp_length.^3)*2 + 3;


for tempi=1:sum(numSeq)
    tempi_2 = temp_range_reOrder(tempi);
    if valid_boostMoreThan0_index(tempi_2) == 1
        scatter(gAcc_noChoice_collapsed_inOne(tempi_2), gAcc_choice_collapsed_inOne(tempi_2), ...
            temp_size(tempi_2), 'filled', 'MarkerFaceColor', temp_color(1, :), ...
            'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);%0.6
        hold on
        
    elseif valid_boostLessThan0_index(tempi_2) == 1
        scatter(gAcc_noChoice_collapsed_inOne(tempi_2), gAcc_choice_collapsed_inOne(tempi_2), ...
            temp_size(tempi_2), 'filled', 'MarkerFaceColor', temp_color(2, :), ...
            'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);%0.6
        hold on
    end
end
legend('Diagonal','Positive boost','Negative boost','Location','southeast','fontsize',8)


xlim([0 1]);
ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Accuracy (noChoice)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Accuracy (choice)', 'FontSize', 18, 'FontWeight', 'bold');
temp_str1 = sprintf('Correlation of accuracy, noChoice & choice');
% temp_str2 = sprintf('Boost=[%.2f,%.2f]',boost_min,boost_max);
temp_str2 = sprintf('Size represent length');
[t,s] = title(temp_str1,temp_str2);
t.FontSize = 18;
t.FontWeight = 'bold';
s.FontSize = 14;
% t = title(temp_str1);
% t.FontSize = 18;
% t.FontWeight = 'bold';

currentFigName = ['correlation_accuracy_choiceAndNoChoice_positiveAndNegative', '_'];
% to generate a unique file name for saving figure
fileName_fig262 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if false
    exportgraphics(fig262, fileName_fig262, 'BackgroundColor', 'none');
end

% if false
%     %temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\sequenceInfo_Z';
%     temp_str = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\sequenceInfo_D'; %#ok<*UNRCH>
%     save(temp_str, 'valid_boostMoreThan0_index', 'valid_boostLessThan0_index', ...
%         'target_seqSet_inOne', 'valid_index');
% end


%% Compute and plot endingHold.choiceMemory, compare all seqs and negative boost seqs
endingHold;
endingHold.noChoiceMemory.all;
endingHold.noChoiceMemory.all_mean = cell(1, pointKindsNum);
endingHold.noChoiceMemory.all_mean_inOne = [];
endingHold.choiceMemory.all;
endingHold.choiceMemory.all_mean = cell(1, pointKindsNum);
endingHold.choiceMemory.all_mean_inOne = [];
for tempi=1:pointKindsNum
    endingHold.noChoiceMemory.all_mean{tempi} = nan(1, numSeq(tempi));
    endingHold.choiceMemory.all_mean{tempi} = nan(1, numSeq(tempi));
    for tempj=1:numSeq(tempi)
        endingHold.choiceMemory.all_mean{tempi}(tempj) = mean(endingHold.choiceMemory.all{tempi}{tempj});
        endingHold.noChoiceMemory.all_mean{tempi}(tempj) = mean(endingHold.noChoiceMemory.all{tempi}{tempj});
    end
end
for tempi=1:pointKindsNum
    endingHold.noChoiceMemory.all_mean_inOne = [endingHold.noChoiceMemory.all_mean_inOne ...
        endingHold.noChoiceMemory.all_mean{tempi}];
    endingHold.choiceMemory.all_mean_inOne = [endingHold.choiceMemory.all_mean_inOne ...
        endingHold.choiceMemory.all_mean{tempi}];
end


endingHold.choiceMemory.error_mean = cell(1, pointKindsNum);
endingHold.choiceMemory.error_mean_inOne = [];
for tempi=1:pointKindsNum
    endingHold.choiceMemory.error_mean{tempi} = nan(1, numSeq(tempi));
    for tempj=1:numSeq(tempi)
        endingHold.choiceMemory.error_mean{tempi}(tempj) = mean(endingHold.choiceMemory.error{tempi}{tempj});
    end
end
for tempi=1:pointKindsNum
    endingHold.choiceMemory.error_mean_inOne = [endingHold.choiceMemory.error_mean_inOne ...
        endingHold.choiceMemory.error_mean{tempi}];
end

endingHold.choiceMemory.correct_mean = cell(1, pointKindsNum);
endingHold.choiceMemory.correct_mean_inOne = [];
for tempi=1:pointKindsNum
    endingHold.choiceMemory.correct_mean{tempi} = nan(1, numSeq(tempi));
    for tempj=1:numSeq(tempi)
        endingHold.choiceMemory.correct_mean{tempi}(tempj) = mean(endingHold.choiceMemory.correct{tempi}{tempj});
    end
end
for tempi=1:pointKindsNum
    endingHold.choiceMemory.correct_mean_inOne = [endingHold.choiceMemory.correct_mean_inOne ...
        endingHold.choiceMemory.correct_mean{tempi}];
end


temp_valid_boostMoreThan0_index = (choiceMinusNoChoice_inOne'>0.1) & valid_index;
temp_valid_boostLessThan0_index = (choiceMinusNoChoice_inOne'<-0.1) & valid_index;

endingHold.choiceMemory.all_allSeq = endingHold.choiceMemory.all_mean_inOne;
endingHold.choiceMemory.all_positiveSeq = endingHold.choiceMemory.all_allSeq(temp_valid_boostMoreThan0_index);
endingHold.choiceMemory.all_negativeSeq = endingHold.choiceMemory.all_allSeq(temp_valid_boostLessThan0_index);
[h_endingHold_all_positiveSeq_negative,p_endingHold_all_positiveSeq_negative]=...
    ttest2(endingHold.choiceMemory.all_positiveSeq, endingHold.choiceMemory.all_negativeSeq);
% [h_endingHold_all_positiveSeq_negative,p_endingHold_all_positiveSeq_negative]=...
%     ttest2(endingHold.choiceMemory.all_positiveSeq, endingHold.choiceMemory.all_negativeSeq,'Tail', 'left');
p_endingHold_all_positiveSeq_negative;
mean(endingHold.choiceMemory.all_positiveSeq);
mean(endingHold.choiceMemory.all_negativeSeq);


endingHold.choiceMemory.error_allSeq = endingHold.choiceMemory.error_mean_inOne;
endingHold.choiceMemory.error_positiveSeq = endingHold.choiceMemory.error_allSeq(temp_valid_boostMoreThan0_index);
endingHold.choiceMemory.error_negativeSeq = endingHold.choiceMemory.error_allSeq(temp_valid_boostLessThan0_index);
[h_endingHold_error_positiveSeq_negative,p_endingHold_error_positiveSeq_negative]=...
    ttest2(endingHold.choiceMemory.error_positiveSeq, endingHold.choiceMemory.error_negativeSeq,'Tail', 'right');
% [h_endingHold_error_positiveSeq_negative,p_endingHold_error_positiveSeq_negative]=...
%     ttest2(endingHold.choiceMemory.error_positiveSeq, endingHold.choiceMemory.error_negativeSeq);
p_endingHold_error_positiveSeq_negative;
mean(endingHold.choiceMemory.error_positiveSeq);
mean(endingHold.choiceMemory.error_negativeSeq);

endingHold.choiceMemory.correct_allSeq = endingHold.choiceMemory.correct_mean_inOne;
endingHold.choiceMemory.correct_positiveSeq = endingHold.choiceMemory.correct_allSeq(temp_valid_boostMoreThan0_index);
endingHold.choiceMemory.correct_negativeSeq = endingHold.choiceMemory.correct_allSeq(temp_valid_boostLessThan0_index);
% [h_endingHold_correct_positiveSeq_negative,p_endingHold_correct_positiveSeq_negative]=...
%     ttest2(endingHold.choiceMemory.correct_positiveSeq, ...
%     endingHold.choiceMemory.correct_negativeSeq(~isnan(endingHold.choiceMemory.correct_negativeSeq)),'Tail', 'right');
[h_endingHold_correct_positiveSeq_negative,p_endingHold_correct_positiveSeq_negative]=...
    ttest2(endingHold.choiceMemory.correct_positiveSeq, ...
    endingHold.choiceMemory.correct_negativeSeq(~isnan(endingHold.choiceMemory.correct_negativeSeq)));

endingHold.choiceMemory.delta_allSeq = endingHold.choiceMemory.correct_allSeq - endingHold.choiceMemory.error_allSeq;
endingHold.choiceMemory.delta_positiveSeq = ...
    endingHold.choiceMemory.correct_positiveSeq - ...
    endingHold.choiceMemory.error_positiveSeq;
endingHold.choiceMemory.delta_negativeSeq = ...
    endingHold.choiceMemory.correct_negativeSeq - ...
    endingHold.choiceMemory.error_negativeSeq;
endingHold.choiceMemory.delta_negativeSeq_valid = ...
    endingHold.choiceMemory.delta_negativeSeq(~isnan(endingHold.choiceMemory.delta_negativeSeq));
% [h_endingHold_delta_positiveSeq_negative,p_endingHold_delta_positiveSeq_negative]=...
%     ttest2(endingHold.choiceMemory.delta_positiveSeq, endingHold.choiceMemory.delta_negativeSeq_valid);
[h_endingHold_delta_positiveSeq_negative,p_endingHold_delta_positiveSeq_negative]=...
    ttest2(endingHold.choiceMemory.delta_positiveSeq, endingHold.choiceMemory.delta_negativeSeq_valid,'Tail', 'left');
p_endingHold_delta_positiveSeq_negative;

mean(endingHold.choiceMemory.delta_positiveSeq);
mean(endingHold.choiceMemory.delta_negativeSeq_valid);


%% Figure27
if false
    fig27 = figure(27);
    set(gcf,'Position',[0 50 1650 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout('flow','TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    
    % temp_p = [p_endingHold_all_positiveSeq_negative,...
    %     p_endingHold_error_positiveSeq_negative,...
    %     p_endingHold_delta_positiveSeq_negative];
    temp_p = [p_endingHold_correct_positiveSeq_negative,...
        p_endingHold_error_positiveSeq_negative,...
        p_endingHold_delta_positiveSeq_negative];
    for tempi=1:3
        ax = nexttile;
        
        if tempi == 1
            %         temp_str1 = 'all_positiveSeq';
            %         temp_str2 = 'all_negativeSeq';
            %         temp_str3 = 'all response';
            temp_str1 = 'correct_positiveSeq';
            temp_str2 = 'correct_negativeSeq';
            temp_str3 = 'correct response';
        elseif tempi == 2
            temp_str1 = 'error_positiveSeq';
            temp_str2 = 'error_negativeSeq';
            temp_str3 = 'error response';
        elseif tempi == 3
            temp_str1 = 'delta_positiveSeq';
            temp_str2 = 'delta_negativeSeq_valid';
            temp_str3 = 'correct - error';
        end
        temp_1 = endingHold.choiceMemory.(temp_str1);
        temp_2 = endingHold.choiceMemory.(temp_str2);
        temp_2 = temp_2(~isnan(temp_2));
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        temp_p(tempi);
        
        bar([0 1], [mean(temp_1) mean(temp_2)], ...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        errorbar([0 1], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
        tempTxt = sprintf('');
        if temp_p(tempi) < 0.001
            tempTxt = sprintf('***');
        elseif temp_p(tempi) < 0.01
            tempTxt = sprintf('**');
        elseif temp_p(tempi) < 0.05
            tempTxt = sprintf('*');
        end
        text(0.5,max([mean(temp_1) mean(temp_2)])+0.1,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        if tempi == 3
            ylim_min = min([mean(temp_1) mean(temp_2)])-0.3;
            ylim_max = max([mean(temp_1) mean(temp_2)])+0.3;
            ylim([ylim_min ylim_max]);
            %set(gca,'YTick',[ylim_min:0.2:ylim_max]);%设置要显示坐标刻度的范围
        else
            ylim([0 max([mean(temp_1) mean(temp_2)])+0.2]);
            set(gca,'YTick',[0:0.2:1.3]);%设置要显示坐标刻度的范围
        end
        set(gca, 'FontSize', 20)
        set(gca,'XTickLabel', ["Positive seqs"; "Negative seqs"],'FontSize', 12);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        ylabel('Reaction time (s)', 'FontSize', 18, 'FontWeight', 'bold');
        title({'\fontsize{15}EndingHold RT, choice memory';[monkey_name,', ',temp_str3]});
    end
    currentFigName = ['endingHold_positive_negative', '_'];
    % to generate a unique file name for saving figure
    fileName_fig27 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if ifSave_fig(27) == 1
        exportgraphics(fig27, fileName_fig27, 'BackgroundColor', 'none');
    end
    
    
    fig272 = figure(272);
    set(gcf,'Position',[0 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    backgroundColor = [1 1 1];
    
    p_correct = p_endingHold_correct_positiveSeq_negative;
    p_error = p_endingHold_error_positiveSeq_negative;
    
    set(gca,'color',backgroundColor);
    
    
    temp_correct_positiveSeq = endingHold.choiceMemory.correct_positiveSeq;
    temp_error_positiveSeq = endingHold.choiceMemory.error_positiveSeq;
    
    temp_correct_negativeSeq = endingHold.choiceMemory.correct_negativeSeq;
    temp_correct_negativeSeq = temp_correct_negativeSeq(~isnan(temp_correct_negativeSeq));
    temp_error_negativeSeq = endingHold.choiceMemory.error_negativeSeq;
    temp_error_negativeSeq = temp_error_negativeSeq(~isnan(temp_error_negativeSeq));
    
    temp_mean_positiveSeq = [mean(temp_correct_positiveSeq) mean(temp_error_positiveSeq)];
    temp_mean_negativeSeq = [mean(temp_correct_negativeSeq) mean(temp_error_negativeSeq)];
    
    y_bar = [temp_mean_positiveSeq; temp_mean_negativeSeq]';
    
    temp_correct_positiveSeq_SEM = std(temp_correct_positiveSeq)/sqrt(length(temp_correct_positiveSeq));
    temp_error_positiveSeq_SEM = std(temp_error_positiveSeq)/sqrt(length(temp_error_positiveSeq));
    
    temp_correct_negativeSeq_SEM = std(temp_correct_negativeSeq)/sqrt(length(temp_correct_negativeSeq));
    temp_error_negativeSeq_SEM = std(temp_error_negativeSeq)/sqrt(length(temp_error_negativeSeq));
    
    temp_SEM_positiveSeq = [temp_correct_positiveSeq_SEM temp_error_positiveSeq_SEM];
    temp_SEM_negativeSeq = [temp_correct_negativeSeq_SEM temp_error_negativeSeq_SEM];
    
    bar_width = 0.8;
    b = bar([0 1],y_bar,bar_width,'grouped',...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    b(1).FaceColor = [0.75 0.25 0.25];
    b(2).FaceColor = [0.25 0.25 0.75];
    b(1).FaceAlpha = 0.4;
    b(2).FaceAlpha = 0.4;
    
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(y_bar);
    % Get the x coordinate of the bars
    x_bar = nan(ngroups, nbars);
    for i = 1:nbars
        x_bar(:,i) = b(i).XEndPoints;
    end
    
    errorbar(x_bar, y_bar,[temp_SEM_positiveSeq; temp_SEM_negativeSeq]', '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
    hold on
    
    tempTxt_correct = sprintf('');
    if p_correct < 0.001
        tempTxt_correct = sprintf('***');
    elseif p_correct < 0.01
        tempTxt_correct = sprintf('**');
    elseif p_correct < 0.05
        tempTxt_correct = sprintf('*');
    end
    
    tempTxt_error = sprintf('');
    if p_error < 0.001
        tempTxt_error = sprintf('***');
    elseif p_error < 0.01
        tempTxt_error = sprintf('**');
    elseif p_error < 0.05
        tempTxt_error = sprintf('*');
    end
    
    text(0,max(y_bar(1,:))+0.07,tempTxt_correct, 'fontsize',15,'FontWeight','bold','HorizontalAlignment','center');
    text(1,max(y_bar(2,:))+0.07,tempTxt_error, 'fontsize',15,'FontWeight','bold','HorizontalAlignment','center');
    
    legend('Positive','Negative',...
        'Location','northwest','fontsize',10)
    
    temp_category1 = 'Correct';
    temp_category2 = 'Error  ';
    set(gca,'XTick',0:1);
    set(gca,'XTickLabel', [temp_category1; temp_category2]);%给坐标加标签
    ylim([0 ceil((max(max(y_bar))+0.1)*10)/10]);
    temp_ylim = ylim;
    set(gca,'YTick',[0:0.2:temp_ylim(2)]);%设置要显示坐标刻度的范围
    set(gca, 'FontSize', 16)
    set(gca,'box','off');% 取消右、上边框
    ylabel('Reaction time (s)', 'FontSize', 16, 'FontWeight', 'bold');
    
    temp_str1 = 'EndingHold RT, choice memory';
    temp_str2_1 = ['Monkey ',monkey_name];
    temp_str2_2 = sprintf(', p(correct)=%.3f, p(error)=%.3f',p_correct,p_error);
    temp_str2 = [temp_str2_1, temp_str2_2];
    [t,s] = title(temp_str1,temp_str2);
    t.FontSize = 14;
    t.FontWeight = 'bold';
    s.FontSize = 12;
    
    currentFigName = ['endingHold_positive_negative_inOne', '_'];
    % to generate a unique file name for saving figure
    fileName_fig272 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if false
        exportgraphics(fig272, fileName_fig272, 'BackgroundColor', 'none');
    end
    
    
    % endingHold_noChoiceMemory_correct_seqLevel
    temp1 = endingHold.noChoiceMemory.correct;
    temp1 = [];
    for tempi=1:length(numSeq)
        temp1 = [temp1 endingHold.noChoiceMemory.correct{tempi}];
    end
    temp2 = [];
    for tempi=1:sum(numSeq)
        temp2 = [temp2 mean(temp1{tempi},'omitnan')];
    end
    endingHold_noChoiceMemory_correct_seqLevel = temp2;
    
    % endingHold_noChoiceMemory_error_seqLevel
    temp1 = endingHold.noChoiceMemory.error;
    temp1 = [];
    for tempi=1:length(numSeq)
        temp1 = [temp1 endingHold.noChoiceMemory.error{tempi}];
    end
    temp2 = [];
    for tempi=1:sum(numSeq)
        temp2 = [temp2 mean(temp1{tempi},'omitnan')];
    end
    endingHold_noChoiceMemory_error_seqLevel = temp2;
    
    % endingHold_choiceMemory_correct_seqLevel
    temp1 = endingHold.choiceMemory.correct;
    temp1 = [];
    for tempi=1:length(numSeq)
        temp1 = [temp1 endingHold.choiceMemory.correct{tempi}];
    end
    temp2 = [];
    for tempi=1:sum(numSeq)
        temp2 = [temp2 mean(temp1{tempi},'omitnan')];
    end
    endingHold_choiceMemory_correct_seqLevel = temp2;
    
    % endingHold_choiceMemory_error_seqLevel
    temp1 = endingHold.choiceMemory.error;
    temp1 = [];
    for tempi=1:length(numSeq)
        temp1 = [temp1 endingHold.choiceMemory.error{tempi}];
    end
    temp2 = [];
    for tempi=1:sum(numSeq)
        temp2 = [temp2 mean(temp1{tempi},'omitnan')];
    end
    endingHold_choiceMemory_error_seqLevel = temp2;
    
    
    endingHold_choiceMemory_correct_seqLevel;
    endingHold_noChoiceMemory_correct_seqLevel;
    endingHold_choiceMemory_error_seqLevel;
    endingHold_noChoiceMemory_error_seqLevel;
    
    [h_endingHold_correct_choice_noChocie_seqLevel,p_endingHold_correct_choice_noChocie_seqLevel]=...
        ttest(endingHold_choiceMemory_correct_seqLevel(valid_index), endingHold_noChoiceMemory_correct_seqLevel(valid_index));
    
    [h_endingHold_error_choice_noChocie_seqLevel,p_endingHold_error_choice_noChocie_seqLevel]=...
        ttest(endingHold_choiceMemory_error_seqLevel(valid_index), endingHold_noChoiceMemory_error_seqLevel(valid_index));
    
    p_endingHold_correct_choice_noChocie_seqLevel;
    p_endingHold_error_choice_noChocie_seqLevel;
    
    fig273 = figure(273);
    set(gcf,'Position',[0 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    backgroundColor = [1 1 1];
    
    endingHold;
    
    p_correct = p_endingHold_correct_choice_noChocie_seqLevel;
    p_error = p_endingHold_error_choice_noChocie_seqLevel;
    
    set(gca,'color',backgroundColor);
    
    
    endingHold_choiceMemory_correct_seqLevel;
    endingHold_noChoiceMemory_correct_seqLevel;
    endingHold_choiceMemory_error_seqLevel;
    endingHold_noChoiceMemory_error_seqLevel;
    
    temp_mean_choiceMemory = [mean(endingHold_choiceMemory_correct_seqLevel(valid_index),'omitnan') mean(endingHold_choiceMemory_error_seqLevel(valid_index),'omitnan')];
    temp_mean_noChoiceMemory = [mean(endingHold_noChoiceMemory_correct_seqLevel(valid_index),'omitnan') mean(endingHold_noChoiceMemory_error_seqLevel(valid_index),'omitnan')];
    
    % y_bar = [mean(endingHold_choiceMemory_correct_seqLevel,'omitnan') mean(endingHold_noChoiceMemory_correct_seqLevel,'omitnan')];
    y_bar = [temp_mean_choiceMemory; temp_mean_noChoiceMemory]';
    
    temp_choiceMemory_correct_SEM = std(endingHold_choiceMemory_correct_seqLevel(valid_index),'omitnan')/sqrt(sum(~isnan(endingHold_choiceMemory_correct_seqLevel(valid_index))));
    temp_noChoiceMemory_correct_SEM = std(endingHold_noChoiceMemory_correct_seqLevel(valid_index),'omitnan')/sqrt(sum(~isnan(endingHold_noChoiceMemory_correct_seqLevel(valid_index))));
    temp_choiceMemory_error_SEM = std(endingHold_choiceMemory_error_seqLevel(valid_index),'omitnan')/sqrt(sum(~isnan(endingHold_choiceMemory_error_seqLevel(valid_index))));
    temp_noChoiceMemory_error_SEM = std(endingHold_noChoiceMemory_error_seqLevel(valid_index),'omitnan')/sqrt(sum(~isnan(endingHold_noChoiceMemory_error_seqLevel(valid_index))));
    
    temp_SEM_choiceMemory = [temp_choiceMemory_correct_SEM temp_choiceMemory_error_SEM];
    temp_SEM_noChoiceMemory = [temp_noChoiceMemory_correct_SEM temp_noChoiceMemory_error_SEM];
    
    
    bar_width = 0.8;
    b = bar([0 1],y_bar,bar_width,'grouped',...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    b(1).FaceColor = [0.4660 0.6740 0.1880];
    b(2).FaceColor = [0.45 0.45 0.45];
    
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(y_bar);
    % Get the x coordinate of the bars
    x_bar = nan(ngroups, nbars);
    for i = 1:nbars
        x_bar(:,i) = b(i).XEndPoints;
    end
    
    errorbar(x_bar, y_bar,[temp_SEM_choiceMemory; temp_SEM_noChoiceMemory]', '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
    hold on
    
    tempTxt_correct = sprintf('');
    if p_correct < 0.001
        tempTxt_correct = sprintf('***');
    elseif p_correct < 0.01
        tempTxt_correct = sprintf('**');
    elseif p_correct < 0.05
        tempTxt_correct = sprintf('*');
    end
    
    tempTxt_error = sprintf('');
    if p_error < 0.001
        tempTxt_error = sprintf('***');
    elseif p_error < 0.01
        tempTxt_error = sprintf('**');
    elseif p_error < 0.05
        tempTxt_error = sprintf('*');
    end
    
    text(0,max(y_bar(1,:))+0.07,tempTxt_correct, 'fontsize',15,'FontWeight','bold','HorizontalAlignment','center');
    text(1,max(y_bar(2,:))+0.07,tempTxt_error, 'fontsize',15,'FontWeight','bold','HorizontalAlignment','center');
    
    legend('choice','noChoice',...
        'Location','northwest','fontsize',10)
    
    temp_category1 = 'Correct';
    temp_category2 = 'Error  ';
    set(gca,'XTick',0:1);
    set(gca,'XTickLabel', [temp_category1; temp_category2]);%给坐标加标签
    ylim([0 ceil((max(max(y_bar))+0.1)*10)/10]);
    temp_ylim = ylim;
    set(gca,'YTick',[0:0.2:temp_ylim(2)]);%设置要显示坐标刻度的范围
    set(gca, 'FontSize', 16)
    set(gca,'box','off');% 取消右、上边框
    ylabel('Reaction time (s)', 'FontSize', 16, 'FontWeight', 'bold');
    
    temp_str1 = 'EndingHold RT';
    temp_str2_1 = ['Monkey ',monkey_name];
    temp_str2_2 = sprintf(', p(correct)=%.3f, p(error)=%.3f',p_correct,p_error);
    temp_str2 = [temp_str2_1, temp_str2_2];
    [t,s] = title(temp_str1,temp_str2);
    t.FontSize = 14;
    t.FontWeight = 'bold';
    s.FontSize = 12;
    
end

%% Figure 28
fig28 = figure(28);

% gError_noChoice_collapsed_inOne = 1-gAcc_noChoice_collapsed_inOne;
% x = gError_noChoice_collapsed_inOne(valid_index);
% x = offloadingProb_inOne(valid_index);
% y = endingHold.noChoiceMemory.all_mean_inOne(valid_index);
x = offloadingProb_inOne;
y = endingHold.noChoiceMemory.all_mean_inOne;
n = 1;
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
temp_color = [0.25 0.25 0.25];

scatter(x, y, ...
    25, 'filled', 'MarkerFaceColor', temp_color, ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
hold on

plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
% text(0.1,max(y),sprintf('r2 = %.3f',r2), 'fontsize',16,'FontWeight','bold');
text(0.1,max(y)-0.1,sprintf('r = %.3f, %s',r,tempTxt), 'fontsize',14,'FontWeight','bold');

xlim([0 1]);
ylim([0 max(y)+0.1]);
set(gca, 'FontSize', 20)
set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
% set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
set(gca,'box','off');% 取消右、上边框
xlabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Reaction time (s)', 'FontSize', 18, 'FontWeight', 'bold');
title('Correlation, Ending hold RT & offloading rate', 'FontSize', 18, 'FontWeight', 'bold');


%% Save data of paper fig3
if false == 1
    currentName = ['fig3_data', '_'];
    fileName_figData = getFileNameMonkey_MAT_append(subject_name, currentName, monkey_name);
    
    fig3_panelB_data = struct;
    fig3_panelB_data.offloadingProb_inOne = offloadingProb_inOne;
    fig3_panelB_data.endingHold_noChoiceMemory = endingHold.noChoiceMemory.all_mean_inOne;
    fig3_panelB_data.r = r;
    fig3_panelB_data.p_corr = p_corr;
    fig3_panelB_data.numSeq = numSeq;
    save(fileName_figData, 'fig3_panelB_data','-append');
end

%% Figure 29
fig29 = figure(29);
set(gcf,'Position',[0 50 1650 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout('flow','TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

endingHold.choiceMemory.all_allSeq;
endingHold.choiceMemory.error_allSeq;
endingHold.choiceMemory.delta_allSeq;
choiceMinusNoChoice_inOne';
valid_index;

y_max = max([endingHold.choiceMemory.all_allSeq(valid_index),...
    endingHold.choiceMemory.error_allSeq(valid_index),endingHold.choiceMemory.delta_allSeq(valid_index)]);
y_min = min([endingHold.choiceMemory.all_allSeq(valid_index),...
    endingHold.choiceMemory.error_allSeq(valid_index),endingHold.choiceMemory.delta_allSeq(valid_index)]);
x_max = max(choiceMinusNoChoice_inOne(valid_index));
x_min = min(choiceMinusNoChoice_inOne(valid_index));

temp_valid_boostMoreThan0_index = (choiceMinusNoChoice_inOne'>0.1) & valid_index;
temp_valid_boostLessThan0_index = (choiceMinusNoChoice_inOne'<-0.1) & valid_index;

fig29_valid_index = valid_index;
% fig29_valid_index = temp_valid_boostMoreThan0_index ...
%     | temp_valid_boostLessThan0_index;

for tempi=1:3
    ax = nexttile;
    
    x = choiceMinusNoChoice_inOne(fig29_valid_index)';
    if tempi == 1
        y = endingHold.choiceMemory.all_allSeq(fig29_valid_index);
        temp_str = 'all response';
    elseif tempi == 2
        y = endingHold.choiceMemory.error_allSeq(fig29_valid_index);
        temp_str = 'error response';
    elseif tempi == 3
        y = endingHold.choiceMemory.delta_allSeq(fig29_valid_index);
        temp_str = 'correct - error';
    end
    temp_valid_index = ~isnan(y);
    y = y(temp_valid_index);
    x = x(temp_valid_index);
    
    n = 1;
    [p_mapping,S] = polyfit(x,y,n);
    r2 = 1 - (S.normr/norm(y - mean(y)))^2;
    [r, p_corr] = corr(x',y');
    
    x_fit = min(x):0.1:max(x);
    y_fit = polyval(p_mapping,x_fit);
    temp_color = [0.25 0.25 0.25];
    
    scatter(x, y, ...
        25, 'filled', 'MarkerFaceColor', temp_color, ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
    text(x_min-0.05,y_min+0.35,sprintf('r2 = %.3f',r2), 'fontsize',14,'FontWeight','bold');
    text(x_min-0.05,y_min+0.1,sprintf('r=%.3f, p=%.3f',r,p_corr), 'fontsize',14,'FontWeight','bold');
    
    xlim([x_min-0.1 x_max+0.1]);
    ylim([y_min-0.1 y_max+0.1]);
    
    set(gca, 'FontSize', 20)
    %set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    % set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    set(gca,'box','off');% 取消右、上边框
    xlabel('Boost effect', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Reaction time (s)', 'FontSize', 18, 'FontWeight', 'bold');
    title({'\fontsize{15}Correlation, Ending hold RT & error rate';[monkey_name,', ',temp_str]});
end
currentFigName = ['endingHold_errorRate_correlation', '_'];
% to generate a unique file name for saving figure
fileName_fig29 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if ifSave_fig(29) == 1
    exportgraphics(fig29, fileName_fig29, 'BackgroundColor', 'none');
end


if if_computePupil == 1
    
    %% Figure 30
    fig30 = figure(30);
    
    % gError_noChoice_collapsed_inOne = 1-gAcc_noChoice_collapsed_inOne;
    % x = gError_noChoice_coll apsed_inOne(valid_index);
    % y = pupilSize_mean_seqSet(valid_index);
    
    % temp_num = length(pupilSize_mean_seqSet);
    % x = gError_noChoice_collapsed_inOne(1:temp_num);
    
    % x = gError_noChoice_collapsed_inOne;
    % y = pupilSize_mean_seqSet;
    
    % x = offloadingProb_inOne(valid_index);
    % y = pupilSize_mean_seqSet(valid_index);
    x = offloadingProb_inOne;
    y = pupilSize_mean_seqSet;
    x = x(1:length(y));
    n = 1;
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
    temp_color = [0.25 0.25 0.25];
    
    scatter(x, y, ...
        25, 'filled', 'MarkerFaceColor', temp_color, ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    hold on
    
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', [0.35 0.35 0.35 0.7]);
    % text(0.1,min(y)+0.1,sprintf('r2 = %.3f',r2), 'fontsize',16,'FontWeight','bold');
    % text(0.1,min(y)+0.03,sprintf('r2 = %.3f',r2), 'fontsize',16,'FontWeight','bold');
    % text(0.1,min(y),sprintf('r = %.3f',r), 'fontsize',16,'FontWeight','bold');
    text(0.1,min(y)+0.05*(max(y)-min(y)),sprintf('r = %.3f, %s',r,tempTxt), 'fontsize',14,'FontWeight','bold');
    
    xlim([0 1]);
    % ylim([0 max(y)+0.1]);
    set(gca, 'FontSize', 20)
    set(gca,'XTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    % set(gca,'YTick',[0:0.2:1]);%设置要显示坐标刻度的范围
    set(gca,'box','off');% 取消右、上边框
    xlabel('Offloading rate', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Pupil diameter (mm)', 'FontSize', 18, 'FontWeight', 'bold');
    temp_str1 = 'Correlation, Pupil diameter & offloading rate';
    % temp_str2 = 'Only during last target, [filter, minus baseline, mm]';
    % temp_str2 = 'During last item + delay1 450 ms, [filter, minus baseline, mm]';
    temp_str2 = 'During last item + delay1 300 ms, [filter, minus baseline, mm]';
    [t,s] = title(temp_str1,temp_str2);
    t.FontSize = 18;
    t.FontWeight = 'bold';
    s.FontSize = 14;
    
    %% Save data of paper fig3
    if false == 1
        currentName = ['fig3_data', '_'];
        fileName_figData = getFileNameMonkey_MAT_append(subject_name, currentName, monkey_name);
        
        fig3_panelD_data = struct;
        fig3_panelD_data.offloadingProb_inOne = offloadingProb_inOne;
        fig3_panelD_data.pupilSize_mean_seqSet = pupilSize_mean_seqSet;
        fig3_panelD_data.r = r;
        fig3_panelD_data.p_corr = p_corr;
        fig3_panelD_data.numSeq = numSeq;
        save(fileName_figData, 'fig3_panelD_data','-append');
    end
    
    % %% Figure 31
    % pupilSize_mean_seqSet;
    % gAcc_noChoice_collapsed_inOne;
    % % temp_lowAcc_index = gAcc_noChoice_collapsed_inOne < 0.5;
    % % temp_highAcc_index = ~temp_lowAcc_index;
    % temp_lowAcc_index = gAcc_noChoice_collapsed_inOne < 0.3;
    % temp_highAcc_index = gAcc_noChoice_collapsed_inOne > 0.7;
    %
    % temp_mean = [mean(pupilSize_mean_seqSet(temp_lowAcc_index)) ...
    %     mean(pupilSize_mean_seqSet(temp_highAcc_index))];
    % temp_sem = zeros(1, 2);
    % temp_sem(1) = std(pupilSize_mean_seqSet(temp_lowAcc_index))./sqrt(sum(temp_lowAcc_index));
    % temp_sem(2) = std(pupilSize_mean_seqSet(temp_highAcc_index))./sqrt(sum(temp_highAcc_index));
    %
    % [h_pupilSize_lowAcc_highAcc,p_pupilSize_lowAcc_highAcc]=...
    %     ttest2(pupilSize_mean_seqSet(temp_lowAcc_index), pupilSize_mean_seqSet(temp_highAcc_index));
    %
    %
    % fig31 = figure(31);
    %
    % bar([0 1], temp_mean, 'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
    % hold on
    % errorbar([0 1],temp_mean, temp_sem, '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
    % tempTxt = sprintf('');
    % if p_pupilSize_lowAcc_highAcc < 0.001
    %     tempTxt = sprintf('***');
    % elseif p_pupilSize_lowAcc_highAcc < 0.01
    %     tempTxt = sprintf('**');
    % elseif p_pupilSize_lowAcc_highAcc < 0.05
    %     tempTxt = sprintf('*');
    % end
    % text(0.5,max(pupilSize_mean_seqSet)+0.1,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
    %     'HorizontalAlignment','center');
    %
    % % ylim([0 max(pupilSize_mean_seqSet)+0.2]);
    % % set(gca,'YTick',[0:0.2:1.3]);%设置要显示坐标刻度的范围
    % set(gca, 'FontSize', 20)
    % set(gca,'XTickLabel', ["Low accuracy"; "High accuracy"], 'FontSize', 14);%给坐标加标签
    % set(gca,'box','off');% 取消右、上边框
    % ylabel('Pupil diameter (mm)', 'FontSize', 18, 'FontWeight', 'bold');
    % title('Pupil diameter & accuracy', 'FontSize', 18, 'FontWeight', 'bold');
    
    choiceMinusNoChoice_inOne;
    valid_index;
    
    positive_seq_valid_index = (choiceMinusNoChoice_inOne'>=0 & valid_index);
    negative_seq_valid_index = (choiceMinusNoChoice_inOne'<0 & valid_index);
    
    positive_negative = struct;
    positive_negative.positive_seq_valid_index = positive_seq_valid_index;
    positive_negative.negative_seq_valid_index = negative_seq_valid_index;
    positive_negative.valid_index = valid_index;
    positive_negative.choiceMinusNoChoice_inOne = choiceMinusNoChoice_inOne;
    positive_negative.target_seqSet_inOne = target_seqSet_inOne;
    
    
    positive_negative_pupilSize;
    positive_negative_seq_pupilSize;
    positive_negative_sum = struct;
    % positive_negative_sum.positive_correct = positive_negative_seq_pupilSize.positive_correct_mean;
    % positive_negative_sum.positive_error = positive_negative_seq_pupilSize.positive_error_mean;
    % positive_negative_sum.positive_boost = choiceMinusNoChoice_inOne(positive_seq_valid_index)';
    % positive_negative_sum.positive_seq = target_seqSet_inOne(positive_seq_valid_index)';
    % positive_negative_sum.negative_correct = positive_negative_seq_pupilSize.negative_correct_mean;
    % positive_negative_sum.negative_error = positive_negative_seq_pupilSize.negative_error_mean;
    % positive_negative_sum.negative_boost = choiceMinusNoChoice_inOne(negative_seq_valid_index)';
    % positive_negative_sum.negative_seq = target_seqSet_inOne(negative_seq_valid_index)';
    
    positive_negative_sum.positive_correct = positive_negative_choiceMemory_seq_pupilSize_ROI.positive_correct_mean;
    positive_negative_sum.positive_error = positive_negative_choiceMemory_seq_pupilSize_ROI.positive_error_mean;
    positive_negative_sum.positive_boost = choiceMinusNoChoice_inOne(positive_seq_valid_index)';
    positive_negative_sum.positive_seq = target_seqSet_inOne(positive_seq_valid_index)';
    positive_negative_sum.negative_correct = positive_negative_choiceMemory_seq_pupilSize_ROI.negative_correct_mean;
    positive_negative_sum.negative_error = positive_negative_choiceMemory_seq_pupilSize_ROI.negative_error_mean;
    positive_negative_sum.negative_boost = choiceMinusNoChoice_inOne(negative_seq_valid_index)';
    positive_negative_sum.negative_seq = target_seqSet_inOne(negative_seq_valid_index)';
    
    positive_negative_sum_partial = positive_negative_sum;
    
    
    positive_seq_valid_index_partial = (choiceMinusNoChoice_inOne'>=0 & valid_index);
    negative_seq_valid_index_partial = (choiceMinusNoChoice_inOne'<-0 & valid_index);
    
    find1 = find(positive_seq_valid_index==1);
    find2 = find(positive_seq_valid_index_partial==1);
    positive_partial_index = ismember(find1,find2);
    
    find1 = find(negative_seq_valid_index==1);
    find2 = find(negative_seq_valid_index_partial==1);
    negative_partial_index = ismember(find1,find2);
    
    % positive_negative_sum_partial.positive_correct = positive_negative_seq_pupilSize.positive_correct_mean(positive_partial_index);
    % positive_negative_sum_partial.positive_error = positive_negative_seq_pupilSize.positive_error_mean(positive_partial_index);
    % positive_negative_sum_partial.positive_boost = choiceMinusNoChoice_inOne(positive_seq_valid_index_partial)';
    % positive_negative_sum_partial.positive_seq = target_seqSet_inOne(positive_seq_valid_index_partial)';
    % positive_negative_sum_partial.negative_correct = positive_negative_seq_pupilSize.negative_correct_mean(negative_partial_index);
    % positive_negative_sum_partial.negative_error = positive_negative_seq_pupilSize.negative_error_mean(negative_partial_index);
    % positive_negative_sum_partial.negative_boost = choiceMinusNoChoice_inOne(negative_seq_valid_index_partial)';
    % positive_negative_sum_partial.negative_seq = target_seqSet_inOne(negative_seq_valid_index_partial)';
    
    positive_negative_sum_partial.positive_correct = positive_negative_choiceMemory_seq_pupilSize_ROI.positive_correct_mean(positive_partial_index);
    positive_negative_sum_partial.positive_error = positive_negative_choiceMemory_seq_pupilSize_ROI.positive_error_mean(positive_partial_index);
    positive_negative_sum_partial.positive_boost = choiceMinusNoChoice_inOne(positive_seq_valid_index_partial)';
    positive_negative_sum_partial.positive_seq = target_seqSet_inOne(positive_seq_valid_index_partial)';
    positive_negative_sum_partial.negative_correct = positive_negative_choiceMemory_seq_pupilSize_ROI.negative_correct_mean(negative_partial_index);
    positive_negative_sum_partial.negative_error = positive_negative_choiceMemory_seq_pupilSize_ROI.negative_error_mean(negative_partial_index);
    positive_negative_sum_partial.negative_boost = choiceMinusNoChoice_inOne(negative_seq_valid_index_partial)';
    positive_negative_sum_partial.negative_seq = target_seqSet_inOne(negative_seq_valid_index_partial)';
    
    [h_pupilSize_error_positiveSeq_negative,p_pupilSize_error_positiveSeq_negative]=...
        ttest2(positive_negative_sum_partial.positive_error, positive_negative_sum_partial.negative_error,'Tail', 'left');
    % [h_pupilSize_error_positiveSeq_negative,p_pupilSize_error_positiveSeq_negative]=...
    %     ttest2(positive_negative_sum_partial.positive_error, positive_negative_sum_partial.negative_error);
    p_pupilSize_error_positiveSeq_negative;
    
    [h_pupilSize_correct_positiveSeq_negative,p_pupilSize_correct_positiveSeq_negative]=...
        ttest2(positive_negative_sum_partial.positive_correct, positive_negative_sum_partial.negative_correct);
    
    
    %% Figure32
    fig32 = figure(32);
    set(gcf,'Position',[0 50 1100 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout('flow','TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    
    temp_p = [p_pupilSize_correct_positiveSeq_negative,...
        p_pupilSize_error_positiveSeq_negative];
    for tempi=1:2
        ax = nexttile;
        
        if tempi == 1
            temp_str1 = 'positive_correct';
            temp_str2 = 'negative_correct';
            temp_str3 = 'correct response';
        elseif tempi == 2
            temp_str1 = 'positive_error';
            temp_str2 = 'negative_error';
            temp_str3 = 'error response';
        elseif tempi == 3
            temp_str1 = 'positive_delta';
            temp_str2 = 'negative_delta';
            temp_str3 = 'correct - error';
        end
        temp_str4 = '_mean';
        temp_str5 = '_sem';
        
        
        temp_1 = positive_negative_sum_partial.(temp_str1);
        temp_1 = temp_1(~isnan(temp_1));
        temp_2 = positive_negative_sum_partial.(temp_str2);
        temp_2 = temp_2(~isnan(temp_2));
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        %     temp_1 = positive_negative_pupilSize.([temp_str1 temp_str4]);
        %     temp_2 = positive_negative_pupilSize.([temp_str2 temp_str4]);
        %     temp_2 = temp_2(~isnan(temp_2));
        %
        %     temp1_SEM = positive_negative_pupilSize.([temp_str1 temp_str5]);
        %     temp2_SEM = positive_negative_pupilSize.([temp_str2 temp_str5]);
        
        
        
        temp_p(tempi);
        
        bar([0 1], [mean(temp_1) mean(temp_2)], ...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
        hold on
        errorbar([0 1], [mean(temp_1) mean(temp_2)],[temp1_SEM temp2_SEM], '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
        tempTxt = sprintf('');
        if temp_p(tempi) < 0.001
            tempTxt = sprintf('***');
        elseif temp_p(tempi) < 0.01
            tempTxt = sprintf('**');
        elseif temp_p(tempi) < 0.05
            tempTxt = sprintf('*');
        end
        text(0.5,min([mean(temp_1) mean(temp_2)])-0.05,tempTxt,'Color','black','FontSize',20,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        %     if tempi == 3
        %         ylim_min = min([mean(temp_1) mean(temp_2)])-0.2;
        %         ylim_max = max([mean(temp_1) mean(temp_2)])+0.2;
        %         ylim([ylim_min ylim_max]);
        %         %set(gca,'YTick',[ylim_min:0.2:ylim_max]);%设置要显示坐标刻度的范围
        %     else
        %         ylim([0 max([mean(temp_1) mean(temp_2)])+0.2]);
        %         set(gca,'YTick',[0:0.2:1.3]);%设置要显示坐标刻度的范围
        %     end
        ylim_min = min([mean(temp_1) mean(temp_2)])-0.2;
        ylim_max = max([mean(temp_1) mean(temp_2)])+0.2;
        ylim([ylim_min ylim_max]);
        set(gca, 'FontSize', 20)
        set(gca,'XTickLabel', ["Positive seqs"; "Negative seqs"],'FontSize', 12);%给坐标加标签
        set(gca,'box','off');% 取消右、上边框
        ylabel('Pupil size (mm)', 'FontSize', 18, 'FontWeight', 'bold');
        title({'\fontsize{15}Pupil size';[monkey_name,', ',temp_str3,', delay1 300ms']});
    end
    currentFigName = ['pupilSize_positive_negative', '_'];
    % to generate a unique file name for saving figure
    fileName_fig32 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if false
        exportgraphics(fig32, fileName_fig32, 'BackgroundColor', 'none');
    end
    
    %% Figure 322
    fig322 = figure(322);
    set(gcf,'Position',[0 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    backgroundColor = [1 1 1];
    
    p_correct = p_pupilSize_correct_positiveSeq_negative;
    p_error = p_pupilSize_error_positiveSeq_negative;
    
    set(gca,'color',backgroundColor);
    
    
    temp_correct_positiveSeq = positive_negative_sum_partial.positive_correct;
    temp_correct_positiveSeq = temp_correct_positiveSeq(~isnan(temp_correct_positiveSeq));
    temp_error_positiveSeq = positive_negative_sum_partial.positive_error;
    temp_error_positiveSeq = temp_error_positiveSeq(~isnan(temp_error_positiveSeq));
    
    temp_correct_negativeSeq = positive_negative_sum_partial.negative_correct;
    temp_correct_negativeSeq = temp_correct_negativeSeq(~isnan(temp_correct_negativeSeq));
    temp_error_negativeSeq = positive_negative_sum_partial.negative_error;
    temp_error_negativeSeq = temp_error_negativeSeq(~isnan(temp_error_negativeSeq));
    
    temp_mean_positiveSeq = [mean(temp_correct_positiveSeq) mean(temp_error_positiveSeq)];
    temp_mean_negativeSeq = [mean(temp_correct_negativeSeq) mean(temp_error_negativeSeq)];
    
    y_bar = [temp_mean_positiveSeq; temp_mean_negativeSeq]';
    
    temp_correct_positiveSeq_SEM = std(temp_correct_positiveSeq)/sqrt(length(temp_correct_positiveSeq));
    temp_error_positiveSeq_SEM = std(temp_error_positiveSeq)/sqrt(length(temp_error_positiveSeq));
    
    temp_correct_negativeSeq_SEM = std(temp_correct_negativeSeq)/sqrt(length(temp_correct_negativeSeq));
    temp_error_negativeSeq_SEM = std(temp_error_negativeSeq)/sqrt(length(temp_error_negativeSeq));
    
    temp_SEM_positiveSeq = [temp_correct_positiveSeq_SEM temp_error_positiveSeq_SEM];
    temp_SEM_negativeSeq = [temp_correct_negativeSeq_SEM temp_error_negativeSeq_SEM];
    
    bar_width = 0.8;
    b = bar([0 1],y_bar,bar_width,'grouped',...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1);
    hold on
    b(1).FaceColor = [0.75 0.25 0.25];
    b(2).FaceColor = [0.25 0.25 0.75];
    b(1).FaceAlpha = 0.4;
    b(2).FaceAlpha = 0.4;
    
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(y_bar);
    % Get the x coordinate of the bars
    x_bar = nan(ngroups, nbars);
    for i = 1:nbars
        x_bar(:,i) = b(i).XEndPoints;
    end
    
    errorbar(x_bar, y_bar,[temp_SEM_positiveSeq; temp_SEM_negativeSeq]', '.', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'CapSize', 15);
    hold on
    
    tempTxt_correct = sprintf('');
    if p_correct < 0.001
        tempTxt_correct = sprintf('***');
    elseif p_correct < 0.01
        tempTxt_correct = sprintf('**');
    elseif p_correct < 0.05
        tempTxt_correct = sprintf('*');
    end
    
    tempTxt_error = sprintf('');
    if p_error < 0.001
        tempTxt_error = sprintf('***');
    elseif p_error < 0.01
        tempTxt_error = sprintf('**');
    elseif p_error < 0.05
        tempTxt_error = sprintf('*');
    end
    
    text(0,min(y_bar(1,:))-0.1,tempTxt_correct, 'fontsize',15,'FontWeight','bold','HorizontalAlignment','center');
    text(1,min(y_bar(2,:))-0.1,tempTxt_error, 'fontsize',15,'FontWeight','bold','HorizontalAlignment','center');
    
    legend('Positive','Negative',...
        'Location','southwest','fontsize',10)
    
    temp_category1 = 'Correct';
    temp_category2 = 'Error  ';
    set(gca,'XTick',0:1);
    set(gca,'XTickLabel', [temp_category1; temp_category2]);%给坐标加标签
    ylim([ceil((min(min(y_bar))-0.4)*10)/10 0]);
    temp_ylim = ylim;
    set(gca,'YTick',sort([0:-0.2:temp_ylim(1)]));%设置要显示坐标刻度的范围
    set(gca, 'FontSize', 16)
    set(gca,'box','off');% 取消右、上边框
    ylabel('Pupil size (mm)', 'FontSize', 16, 'FontWeight', 'bold');
    
    temp_str1 = 'Pupil size, delay1 300ms';
    temp_str2_1 = ['Monkey ',monkey_name];
    temp_str2_2 = sprintf(', p(correct)=%.3f, p(error)=%.3f',p_correct,p_error);
    temp_str2 = [temp_str2_1, temp_str2_2];
    [t,s] = title(temp_str1,temp_str2);
    t.FontSize = 14;
    t.FontWeight = 'bold';
    s.FontSize = 12;
    
    currentFigName = ['pupilSize_positive_negative_inOne', '_'];
    % to generate a unique file name for saving figure
    fileName_fig322 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
    % to judge whether save figure or not
    if false
        exportgraphics(fig322, fileName_fig322, 'BackgroundColor', 'none');
    end
    
end

%% Save data of paper fig4
if false == 1
    currentName = ['fig4_data', '_'];
    fileName_figData = getFileNameMonkey_MAT(subject_name, currentName, monkey_name);
    
    fig4_panelA_data = struct;
    fig4_panelB_data = struct;
    fig4_panelC_data = struct;
    fig4_panelD_data = struct;
    
    fig4_panelA_data.gAcc_noChoice_collapsed_inOne = gAcc_noChoice_collapsed_inOne;
    fig4_panelA_data.gAcc_choice_collapsed_inOne = gAcc_choice_collapsed_inOne;
    fig4_panelA_data.valid_boostMoreThan0_index = valid_boostMoreThan0_index;
    fig4_panelA_data.valid_boostLessThan0_index = valid_boostLessThan0_index;
    fig4_panelA_data.numSeq = numSeq;
    fig4_panelA_data.temp_range_reOrder = temp_range_reOrder;
    fig4_panelA_data.temp_size = temp_size;
    
    %     fig4_panelB_data.temp_mean_choice = temp_mean_choice;
    %     fig4_panelB_data.temp_mean_rProb = temp_mean_rProb;
    %     fig4_panelB_data.temp_SEM_choice = temp_SEM_choice;
    %     fig4_panelB_data.temp_SEM_rProb = temp_SEM_rProb;
    %     fig4_panelB_data.p_choice = p_choice;
    %     fig4_panelB_data.p_rProb = p_rProb;
    
    fig4_panelB_data.temp_mean_choice = temp_mean_choice;
    fig4_panelB_data.temp_mean_noChoice = temp_mean_noChoice;
    fig4_panelB_data.temp_SEM_choice = temp_SEM_choice;
    fig4_panelB_data.temp_SEM_noChoice = temp_SEM_noChoice;
    fig4_panelB_data.p_choice = p_choice;
    fig4_panelB_data.p_noChoice = p_noChoice;
    fig4_panelB_data.p_positive = p_positive;
    fig4_panelB_data.p_negative = p_negative;
    
    fig4_panelC_data.endingHold_choiceMemory = endingHold.choiceMemory;
    fig4_panelC_data.p_correct = p_endingHold_correct_positiveSeq_negative;
    fig4_panelC_data.p_error = p_endingHold_error_positiveSeq_negative;
    
    fig4_panelD_data.positive_negative_sum_partial = positive_negative_sum_partial;
    fig4_panelD_data.p_correct = p_pupilSize_correct_positiveSeq_negative;
    fig4_panelD_data.p_error = p_pupilSize_error_positiveSeq_negative;
    
    save(fileName_figData, 'fig4_panelA_data','fig4_panelB_data',...
        'fig4_panelC_data','fig4_panelD_data');
end

%% Interference proportion

fig33 = figure(33);
set(gcf,'Position',[0 50 700 500]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point

interference_responseProportion_underSameLengthButError;

temp_data = cell(1, pointKindsNum);
temp_mean = zeros(1, pointKindsNum);
temp_SEM = zeros(1, pointKindsNum);
temp_mean_collaped = [];

for tempi=1:pointKindsNum
    if tempi == 1
        temp_range = 1:numSeq(1);
    else
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
    end
    temp_range;
    
    temp_data{tempi} = interference_responseProportion_underSameLengthButError(temp_range);
    temp_valid = ~isnan(temp_data{tempi});
    temp_mean(tempi) = mean(temp_data{tempi}(temp_valid));
    temp_SEM(tempi) = std(temp_data{tempi}(temp_valid))...
        ./sqrt( sum(temp_valid) );
end

errorbar(seq_length_rangeHead:seq_length_rangeTail,temp_mean,temp_SEM,...
    '-o','Color',[0 0 0],'LineWidth',1,'CapSize',12,'MarkerSize',5);
hold on


ylim([0 1]);
set(gca, 'FontSize', 20)
set(gca,'XLim',[seq_length_rangeHead-1 seq_length_rangeTail+1]);%X轴的数据显示范围
set(gca,'XTick',[seq_length_rangeHead:1:seq_length_rangeTail]);%设置要显示坐标刻度的范围
set(gca,'XTickLabel',[seq_length_rangeHead:1:seq_length_rangeTail]);%给坐标加标签
set(gca,'box','off');% 取消右、上边框
xlabel('Length', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Proportion', 'FontSize', 18, 'FontWeight', 'bold');
title({'\fontsize{20}{\bf Interference proportion when length correct}'});

currentFigName = ['interferenceProportion', '_'];
% to generate a unique file name for saving figure
fileName_fig33 = getFileNameMonkey(subject_name, currentFigName, monkey_name);
% to judge whether save figure or not
if false
    exportgraphics(fig33, fileName_fig33, 'BackgroundColor', 'none');
end


%% End
cd 'C:\ASDROOT\STUDY\MonkeyDataAnalysis'