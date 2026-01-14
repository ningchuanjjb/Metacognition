% Chuan's 33th script (20260114)
%This script: To load preprocessed eyelink data (in mat file type), conduct trial parcellation and save them
%% Initialization
% clear all %#ok<CLALL>
close all

%% Some flags to control the following script
if_monkey_D0_Z1 = 0;% To decide whether dealing with Ding's data or Zelku's data
ifSavePupilSize = 1;

if_load = 1;


targetPATH = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\fixation';
cd(targetPATH)


path_preprocessed = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\Data\preprocessed';
path_results = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\fixation\data';
% path_edf = 'C:\ASDROOT\STUDY\edfData\mat';
path_edf = 'D:\edfData\mat';

if if_monkey_D0_Z1 == 0
    monkey_name = 'D';
    
    %     searchName = 'preprocessed_DnumFrames8_A';
    %     searchName_mappingParam = 'from22-06-21to22-08-12_D_mappingParam_1';
    %     searchName_eye = 'D_edf';
    
    %     searchName = 'preprocessed_DnumFrames6_Train2';
    %     searchName_mappingParam = 'from23-03-13to23-04-05_D_mappingParam_1';
    %     searchName_eye = 'DTrain2_edf';
    
    searchName = 'preprocessed_DnumFrames6_twoSession2p';
    searchName_eye = 'D2p_edf';
    
    
elseif if_monkey_D0_Z1 == 1
    monkey_name = 'Z';
    
    %     searchName = 'preprocessed_ZnumFrames6_B2'; % 没有error stop
    %     searchName_mappingParam = 'fromf_2022-06-15tof_2022-07-01_Z_mappingParam_1';
    %     searchName_eye = 'Zx_edf';
    
    %     searchName = 'preprocessed_ZnumFrames6_Train1';
    %     searchName_mappingParam = 'from23-03-26to23-04-03_Z_mappingParam_1';
    %     searchName_eye = 'ZTrain1_edf';
    
    searchName = 'preprocessed_ZnumFrames6_twoSession2p';
    searchName_eye = 'Z2p_edf';
    
end

t0 = tic;

%% To load data from files
if if_load == 1
    % Load ptb data
    [total_ptbData, MAT_file_name] = loadMat_multi(searchName, path_preprocessed); %load all file ptb(basic&trial）
    fileSize = length(MAT_file_name);
    
    % Load edf data
    total_edfData = loadMat_multi(searchName_eye, path_edf);
    
end

a = 1;

% clear MAT_file*
% clear searchName_*
% clear targetPATH*


choice_g_index = cell(1, fileSize);
choice_r_index = cell(1, fileSize);
choice_g_index_collapsed = [];
choice_r_index_collapsed = [];

isCorrect_collapsed = [];

sequence = cell(1, fileSize);
sequence_collapsed = [];

eyeBaseline_pupilSize = cell(1, fileSize);
eyeSample_pupilSize = cell(1, fileSize);
eyeDelay1_pupilSize = cell(1, fileSize);

eyeBaseline_pupilSize_collapsed = [];
eyeSample_pupilSize_collapsed = [];
eyeDelay1_pupilSize_collapsed = [];

trial_num_multi = nan(1, fileSize);


%% Get data from each file
subject_name = cell(1, fileSize);
for file_index=1:fileSize
    %% Preparation
    file_name = cell2mat(MAT_file_name(file_index)); % file文件名
    temp_tail = strfind(file_name,'train'); %文件名中‘train’在第几个字符
    subject_name(file_index) = {file_name(length(searchName)+1+2: temp_tail-1)};%日期
    
    edfData = eval(['total_edfData.file',sprintf('%d',file_index),'.edf0_saved']); %edfData为当前file data
    ptbData = eval(['total_ptbData.file',sprintf('%d',file_index)]);
    % Sample rate is 1000 hz.
    %sampleRate = 1000;
    trial_num = ptbData.trial_para.trial_count;
    pointShowTime = ptbData.basic_para.pointShowTime;
    pointShowPWM = ptbData.basic_para.pointShowPWM;
    delay1 = ptbData.basic_para.fixationTime(1:trial_num);
    %delay1_min = floor(min(ptbData.basic_para.fixationTime));
    sequence{file_index} = ptbData.trial_para.currentSequence;
    sequence_collapsed = [sequence_collapsed sequence{file_index}];
    pointKindsNum = ptbData.basic_para.seq_length_rangeTail;
    
    %% Get choiceCondition trial index
    ifSelectOffloading = ptbData.trial_para.ifSelectOffloading; %1 offloading/-1 memory
    choiceCondition_flag = ptbData.trial_para.choiceCondition_flag; % 2 choice /0 nochoice
    
    choice_g_index{file_index} = (choiceCondition_flag == 2 & ifSelectOffloading == -1);
    choice_r_index{file_index} = (choiceCondition_flag == 2 & ifSelectOffloading == 1);
    choice_g_index_collapsed = [choice_g_index_collapsed choice_g_index{file_index}];
    choice_r_index_collapsed = [choice_r_index_collapsed choice_r_index{file_index}];
    
    isCorrect_collapsed = [isCorrect_collapsed,ptbData.trial_para.isCorrect];
    
    
    %% Get marker of each trial
    marker_info_raw = edfData.Messages.info;
    marker_time_raw = edfData.Messages.time;
    marker_time = zeros(1, length(marker_info_raw));
    trial_duration = zeros(1, length(marker_info_raw));
    temp_replicate_index = [];
    temp_replicate_marker = [];
    last_marker_count = 0;
    for tempi=1:length(marker_info_raw)
        k = strfind(marker_info_raw{tempi}, 'TRIALID ');% if find, k shouldn't be empty
        if k == 1
            marker_count = str2double(marker_info_raw{tempi}(length('TRIALID ')+1:end)); % get trial number
            if marker_count == last_marker_count %什么情况下=0呢
                temp_replicate_index = [temp_replicate_index tempi]; %#ok<*AGROW>
                temp_replicate_marker = [temp_replicate_marker marker_count];
            end
            marker_time(marker_count) = marker_time_raw(tempi);% 把marker time（trial onset time）与marker count对应
            last_marker_count = marker_count;
        end
    end
    
    for tempi=1:marker_count-1
        trial_duration(tempi) = marker_time(tempi+1) - marker_time(tempi); %trial duration是后一个trial onset-这个trial onset，也包含了delay和respon
    end
    trial_duration(marker_count:end) = 1234567;
    
    %clear marker_info*
    %clear marker_time_raw*
    temp_marker_time = marker_time(ptbData.basic_para.if_outlierTrial==0);
    marker_time = temp_marker_time; % delete outlierTrial
    
    
    %% Get pupil size timecourse of each trial
    periodBaseline_pupilSize = cell(1, trial_num);    
    periodSample_pupilSize = cell(1, trial_num);
    periodDelay1_pupilSize = cell(1, trial_num);
    
    for tempi=1:trial_num
        % mark here in 20220502, 22:30
        timestamp1_sampleOn = marker_time(tempi); %每个trial onset time
        timestamp2_sampleOff_true = timestamp1_sampleOn + length(sequence{file_index}{tempi})*pointShowTime - 1;% sequence = currentseq// 600*点数 -1 为刺激持续时长
        %timestamp2_sampleOff = timestamp1_sampleOn + length(sequence{file_index}{tempi})*pointShowTime - 1 ...
        %    + floor(pointShowTime*(1-pointShowPWM));
        timestamp2_sampleOff = timestamp1_sampleOn + length(sequence{file_index}{tempi})*pointShowTime - 1 ...
           + floor(pointShowTime*pointShowPWM+1);        
        %timestamp3_delay1Off = timestamp2_sampleOff_true + round(delay1(tempi)) - 1; % delay off 时刻
        timestamp3_delay1Off = timestamp2_sampleOff_true + round(delay1(tempi)); % delay off 时刻
        
        index1_sampleOn = timestamp1_sampleOn-edfData.time(1)+1;  % 实际开始时间，因为开始时刻是time（1）
        index2_sampleOff_true = timestamp2_sampleOff_true-edfData.time(1)+1; %实际结束时间
        index2_sampleOff = timestamp2_sampleOff-edfData.time(1)+1;
        index3_delay1Off = timestamp3_delay1Off-edfData.time(1)+1;
        
        lengthBaseline1 = 700;%400
        lengthBaseline2 = 400;
        periodBaseline_pupilSize{tempi} = edfData.pupilSize(index1_sampleOn-lengthBaseline1:index1_sampleOn+lengthBaseline2-1)';
        periodSample_pupilSize{tempi} = edfData.pupilSize(index1_sampleOn+lengthBaseline2:index2_sampleOff)';
        %periodDelay1_pupilSize{tempi} = edfData.pupilSize(index2_sampleOff_true:index3_delay1Off)';
        periodDelay1_pupilSize{tempi} = edfData.pupilSize(index2_sampleOff+1:index3_delay1Off)';        
    end
    
    
    eyeBaseline_pupilSize{file_index} = periodBaseline_pupilSize;    
    eyeSample_pupilSize{file_index} = periodSample_pupilSize;
    eyeDelay1_pupilSize{file_index} = periodDelay1_pupilSize;
    eyeBaseline_pupilSize_collapsed = [eyeBaseline_pupilSize_collapsed,periodBaseline_pupilSize];    
    eyeSample_pupilSize_collapsed = [eyeSample_pupilSize_collapsed,periodSample_pupilSize];
    eyeDelay1_pupilSize_collapsed = [eyeDelay1_pupilSize_collapsed,periodDelay1_pupilSize];
    
    trial_num_multi(file_index) = trial_num;
    
end
% clear total_*
% clear ptbData edfData
% clear period*

eyeBaseline_pupilSize_collapsed;
eyeSample_pupilSize_collapsed;
eyeDelay1_pupilSize_collapsed;
sequence_collapsed;
choice_g_index_collapsed;
choice_r_index_collapsed;
isCorrect_collapsed;
pointKindsNum;

pupil_extracted = struct;
pupil_extracted.eyeBaseline_pupilSize_collapsed = eyeBaseline_pupilSize_collapsed;
pupil_extracted.eyeSample_pupilSize_collapsed = eyeSample_pupilSize_collapsed;
pupil_extracted.eyeDelay1_pupilSize_collapsed = eyeDelay1_pupilSize_collapsed;
pupil_extracted.sequence_collapsed = sequence_collapsed;
pupil_extracted.choice_g_index_collapsed = choice_g_index_collapsed;
pupil_extracted.choice_r_index_collapsed = choice_r_index_collapsed;
pupil_extracted.isCorrect_collapsed = isCorrect_collapsed;
pupil_extracted.pointKindsNum = pointKindsNum;
pupil_extracted.pointShowTime = pointShowTime;
pupil_extracted.pointShowPWM = pointShowPWM;
pupil_extracted.trial_num_multi = trial_num_multi;

%% Save
if ifSavePupilSize == 1
    %save in a mat file
    currentName = ['_','pupilSize_forAnalysis'];
    %currentName2 = [monkey_name,currentName];
    currentName2 = [searchName_eye,currentName];
    
    
    %currentName2 = ['ZTrain1_pupilSize_forAnalysis'];
    %currentName2 = ['DTrain2_pupilSize_forAnalysis'];
    
    fileName_pupilSize = [path_results,'\',currentName2,'.mat'];
    save(fileName_pupilSize, 'pupil_extracted');    
    
end

fprintf('Pupil size extraction time = %.1f secs.\n',toc(t0));

%% End
cd 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\fixation'