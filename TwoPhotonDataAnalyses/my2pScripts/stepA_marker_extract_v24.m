% Chuan's 27th script (20251214)
% This script: To conduct marker extraction from raw data.
%% Initialization
clear
close all
home;% To scroll down in command window

if_plot = 0;
if_fig4_ImageCorr = 0;

currentSession_multi = string;


% currentSession_multi = [currentSession_multi; '113Recording_20231205A_Zelku_noImage'];
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
% currentSession_multi = [currentSession_multi; '113Recording_20240315A_Zelku_Site06RA_sameFOV0312'];
% currentSession_multi = [currentSession_multi; '113Recording_20240319A_Zelku_Site09E'];
% currentSession_multi = [currentSession_multi; '113Recording_20240320A_Zelku_Site09E_sameFOV0319'];
% currentSession_multi = [currentSession_multi; '113Recording_20240322A_Zelku_Site07B'];
% currentSession_multi = [currentSession_multi; '113Recording_20240323A_Zelku_Site07B_sameFOV0322'];
% currentSession_multi = [currentSession_multi; '113Recording_20240329A_Zelku_Site05B'];
% currentSession_multi = [currentSession_multi; '113Recording_20240330A_Zelku_Site05B_sameFOV0329'];
% currentSession_multi = [currentSession_multi; '113Recording_20240402A_Zelku_Site14B'];
% currentSession_multi = [currentSession_multi; '113Recording_20240403A_Zelku_Site14B_sameFOV0402'];
% currentSession_multi = [currentSession_multi; '113Recording_20240410A_Zelku_Site17B'];
currentSession_multi = [currentSession_multi; '113Recording_20240411A_Zelku_Site17B_sameFOV0410'];


% currentSession_multi = [currentSession_multi; 'test0608'];

% currentSession_multi = [currentSession_multi; '113Recording_20230426A_Ding_Site16'];
% currentSession_multi = [currentSession_multi; '113Recording_20230427A_Ding_Site16_sameFOV0426'];
% currentSession_multi = [currentSession_multi; '113Recording_20230502A_Ding_Site13'];
% currentSession_multi = [currentSession_multi; '113Recording_20230503A_Ding_Site13_sameFOV0502'];
% currentSession_multi = [currentSession_multi; '113Recording_20230504A_Ding_Site02'];
% currentSession_multi = [currentSession_multi; '113Recording_20230508A_Ding_Site02_sameFOV0509'];
% currentSession_multi = [currentSession_multi; '113Recording_20230509A_Ding_Site02'];% 660000 frames, easy to crash
% 
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
% 
% currentSession_multi = [currentSession_multi; '113Recording_20230527A_Ding_Site02B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230528A_Ding_Site02B_sameFOV0527'];
% currentSession_multi = [currentSession_multi; '113Recording_20230530A_Ding_Site05C'];
% currentSession_multi = [currentSession_multi; '113Recording_20230531A_Ding_Site05C_sameFOV0530'];
% currentSession_multi = [currentSession_multi; '113Recording_20230601A_Ding_Site13B'];
% currentSession_multi = [currentSession_multi; '113Recording_20230602A_Ding_Site13B_sameFOV0601'];
% 
% currentSession_multi = [currentSession_multi; '113Recording_20230604A_Ding_Site07'];
% currentSession_multi = [currentSession_multi; '113Recording_20230605A_Ding_Site07_sameFOV0604'];
% currentSession_multi = [currentSession_multi; '113Recording_20230612A_Ding_Site14'];
% currentSession_multi = [currentSession_multi; '113Recording_20230614A_Ding_Site05D'];
% currentSession_multi = [currentSession_multi; '113Recording_20230615A_Ding_Site05D_sameFOV0614'];
% currentSession_multi = [currentSession_multi; '113Recording_20230619A_Ding_Site02C'];
% currentSession_multi = [currentSession_multi; '113Recording_20230620A_Ding_Site05E'];

currentSession_multi(1) = [];
num_FOV = length(currentSession_multi);

t0 = tic;
for tempSessionIndex=1:num_FOV
    currentSession = currentSession_multi{tempSessionIndex};
    
    targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
    cd(targetPATH)
    targetPATH_norm = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
    
    % rawData_path = 'D:\twoPhotonRawData\ToBeProcessed';
    rawData_path = 'R:\twoPhotonRawData\ToBeProcessed';
    marker_2p_name = '2P-marker';
    marker_PTB_name = 'PTB-marker';
    PTB_name = '.mat';
    markerExtractTestName_2p='fun_markerExtract_2p';
    markerExtractTestName_PTB='fun_markerExtract_PTB';
    markerDownSamplingTestName='fun_markerDownSampling';
    markerParse_trialLevelTestName='fun_markerParse_trialLevel';
    
    parellelTestName = 'parellelTest';
    
    currentSession_path = [rawData_path '\' currentSession];
    
    marker_2p_fullName = autoGetMarkerFileName(['*',marker_2p_name], currentSession_path);
    marker_PTB_fullName = autoGetMarkerFileName(['*',marker_PTB_name], currentSession_path);
    % [~,parts_fname,~] = fileparts(marker_PTB_fullName);
    % PTB_fullName = autoGetMarkerFileName([parts_fname(1:4),'*',PTB_name], currentSession_path);
    PTB_fullName = autoGetMarkerFileName(['2*',PTB_name], currentSession_path);
    
    %output_shortPath = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected';
    output_shortPath = 'D:\twoPhotonData_motionCorrected';
    temp_currentSession_path = [output_shortPath '\' currentSession];
    temp_if_max0_min1 = 0;
    
    
    MAT_file=dir([temp_currentSession_path,'\Result*']);
    if isempty(MAT_file) == 1
        datestr_now30 = datestr(now,30);
        output_path = [temp_currentSession_path,'\Result',datestr_now30];
        mkdir(output_path);
    else
        output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
    end
    
    
    %% Convert and load edf file
    % file_name_eyeMat_short = Edf2Mat_jjb_single_v1('edf',currentSession_path,currentSession_path,targetPATH);
    % file_name_eyeMat = fullfile(currentSession_path,[file_name_eyeMat_short,'.mat']);
    % load(file_name_eyeMat,'edf0_saved');
    
    %% Perform marker extraction
    markerExtractTestName_2p_v = autoGetFunName_myScripts(markerExtractTestName_2p, [targetPATH '\functions']);
    fprintf('Now runing %s.  ------> ', markerExtractTestName_2p_v);
    fun_markerExtractTest_2p = str2func(markerExtractTestName_2p_v);
    extractOutput_2p = fun_markerExtractTest_2p(marker_2p_fullName);
    
    markerExtractTestName_PTB_v = autoGetFunName_myScripts(markerExtractTestName_PTB, [targetPATH '\functions']);
    fprintf('Now runing %s. ------> ', markerExtractTestName_PTB_v);
    fun_markerExtractTest_PTB = str2func(markerExtractTestName_PTB_v);
    % extractOutput_PTB = fun_markerExtractTest_PTB(marker_PTB_fullName,PTB_fullName,edf0_saved);
    extractOutput_PTB = fun_markerExtractTest_PTB(marker_PTB_fullName,PTB_fullName);
    
    
    
    
    
    %% Downsampling PTB marker to 2p frame space
    markerDownSamplingTestName_v = autoGetFunName_myScripts(markerDownSamplingTestName, [targetPATH '\functions']);
    fprintf('Now runing %s.  ------> ', markerDownSamplingTestName_v);
    fun_markerDownSamplingTest = str2func(markerDownSamplingTestName_v);
    extractOutput_PTB.TOTAL_downSampling = fun_markerDownSamplingTest(extractOutput_PTB.marker_sampleIndex,extractOutput_2p.valid_frames_sampleIndex);
    extractOutput_PTB.TRIALID_downSampling = fun_markerDownSamplingTest(extractOutput_PTB.marker_TRIALID_sampleIndex,extractOutput_2p.valid_frames_sampleIndex);
    extractOutput_PTB.REWARD_downSampling = fun_markerDownSamplingTest(extractOutput_PTB.marker_REWARD_sampleIndex,extractOutput_2p.valid_frames_sampleIndex);
    extractOutput_PTB.validTRIALID_downSampling = fun_markerDownSamplingTest(extractOutput_PTB.marker_validTRIALID_sampleIndex,extractOutput_2p.valid_frames_sampleIndex);
    fprintf('\n');
    a = 1;
    
    %% Parse marker into trial level
    markerParse_trialLevelTestName_v = autoGetFunName_myScripts(markerParse_trialLevelTestName, [targetPATH '\functions']);
    fprintf('Now runing %s. ------> ', markerParse_trialLevelTestName_v);
    fun_markerParse_trialLevel = str2func(markerParse_trialLevelTestName_v);
    markerParse_trialLevel = fun_markerParse_trialLevel(extractOutput_PTB);
    % fprintf('\n');
    
    save([output_path,'\','markerParse_trialLevel.mat'],'markerParse_trialLevel');
    
    
    
    if if_plot == 1
        %% Plot 2p marker
        fig1 = figure('Name','Fig1, Two photon marker','NumberTitle','off');
        set(gcf,'Position',[35 35+350 700 250]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
        
        nexttile
        plot(extractOutput_2p.marker_array, '-', 'LineWidth', 1, 'Color', [0.40 0.40 0.40]);
        hold on
        plot(extractOutput_2p.validRecordingBlock_range(1)*ones(1,2),[0 1], '-', 'LineWidth', 3, 'Color', [0.75 0.25 0.25]);
        hold on
        plot(extractOutput_2p.validRecordingBlock_range(2)*ones(1,2),[0 1], '-', 'LineWidth', 3, 'Color', [0.75 0.25 0.25]);
        hold on
        plot(extractOutput_2p.validRecordingBlock_range,[1 1], '-', 'LineWidth', 3, 'Color', [0.75 0.25 0.25]);
        hold on
        ylim([0 1]);
        set(gca,'yticklabel',[])
        set(gca, 'FontSize', 11)
        set(gca,'box','off');% 取消右、上边框
        % title(sprintf('NI card sample count'), 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Image sample index (NI card), %d samples',sum(extractOutput_2p.marker_array)), 'FontSize', 12, 'FontWeight', 'bold');
        
        nexttile
        plot(extractOutput_2p.valid_frames_boolIndex, '-', 'LineWidth', 5, 'Color', [0.40 0.40 0.40]);
        hold on
        set(gca,'yticklabel',[])
        ylim([0 1]);
        set(gca, 'FontSize', 11)
        set(gca,'box','off');% 取消右、上边框
        title(sprintf('valid frames index (image), %d frames',sum(extractOutput_2p.valid_frames_boolIndex)), 'FontSize', 12, 'FontWeight', 'bold');
        
        
        %% Plot PTB marker
        fig2 = figure('Name','Fig2, Behavior marker','NumberTitle','off');
        set(gcf,'Position',[35+700 35+350 700 250]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        for tempi=1:length(extractOutput_PTB.marker_TRIALID_sampleIndex)
            plot(extractOutput_PTB.marker_TRIALID_sampleIndex(tempi)*ones(1,2),[0 1], '-', 'LineWidth', 0.01, 'Color', [0.40 0.40 0.40]);
            hold on
        end
        xlim([0 max(extractOutput_PTB.marker_TRIALID_sampleIndex)]);
        ylim([0 1]);
        set(gca,'yticklabel',[])
        set(gca, 'FontSize', 11)
        set(gca,'box','off');% 取消右、上边框
        title(sprintf('Trial start sample index (NI card), %d trials',length(extractOutput_PTB.marker_TRIALID_sampleIndex)), 'FontSize', 12, 'FontWeight', 'bold');
        
        nexttile
        for tempi=1:length(extractOutput_PTB.marker_REWARD_sampleIndex)
            plot(extractOutput_PTB.marker_REWARD_sampleIndex(tempi)*ones(1,2),[0 1], '-', 'LineWidth', 0.01, 'Color', [0.40 0.40 0.40]);
            hold on
        end
        xlim([0 max(extractOutput_PTB.marker_TRIALID_sampleIndex)]);
        ylim([0 1]);
        set(gca,'yticklabel',[])
        set(gca, 'FontSize', 1)
        set(gca,'box','off');% 取消右、上边框
        title(sprintf('Reward sample index (NI card), %d rewards',length(extractOutput_PTB.marker_REWARD_sampleIndex)), 'FontSize', 12, 'FontWeight', 'bold');
        
        
        %% Plot PTB marker in 2p frame space
        fig3 = figure('Name','Fig3, Behavior marker in 2p frame space','NumberTitle','off');
        set(gcf,'Position',[35+700 35 700 250]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
        
        nexttile
        for tempi=1:length(extractOutput_PTB.TRIALID_downSampling.valid_marker_X_frameIndex)
            plot(extractOutput_PTB.TRIALID_downSampling.valid_marker_X_frameIndex(tempi)*ones(1,2),[0 1], '-', 'LineWidth', 0.01, 'Color', [0.40 0.40 0.40]);
            hold on
        end
        xlim([0 max(extractOutput_PTB.TRIALID_downSampling.valid_marker_X_frameIndex)]);
        ylim([0 1]);
        set(gca,'yticklabel',[])
        set(gca, 'FontSize', 11)
        set(gca,'box','off');% 取消右、上边框
        title(sprintf('Trial start frame index (image), %d trials',length(extractOutput_PTB.TRIALID_downSampling.valid_marker_X_frameIndex)), 'FontSize', 12, 'FontWeight', 'bold');
        
        nexttile
        for tempi=1:length(extractOutput_PTB.REWARD_downSampling.valid_marker_X_frameIndex)
            plot(extractOutput_PTB.REWARD_downSampling.valid_marker_X_frameIndex(tempi)*ones(1,2),[0 1], '-', 'LineWidth', 0.01, 'Color', [0.40 0.40 0.40]);
            hold on
        end
        xlim([0 max(extractOutput_PTB.TRIALID_downSampling.valid_marker_X_frameIndex)]);
        ylim([0 1]);
        set(gca,'yticklabel',[])
        set(gca, 'FontSize', 11)
        set(gca,'box','off');% 取消右、上边框
        title(sprintf('Reward frame index (image), %d rewards',length(extractOutput_PTB.REWARD_downSampling.valid_marker_X_frameIndex)), 'FontSize', 12, 'FontWeight', 'bold');
        
        if if_fig4_ImageCorr == 1
            parellelTestName_v = autoGetFunName_myScripts(parellelTestName, [targetPATH_norm '\functions']);
            fprintf('Now runing %s. ', parellelTestName_v);
            fun_parellelTest = str2func(parellelTestName_v);
            
            %output_shortPath = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected';
            %temp_currentSession_path = [output_shortPath '\' currentSession];
            %temp_if_max0_min1 = 0;
            %output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
            
            temp_if_imageAverage = 0;
            temp_finalOutputType = 'suite2pbin';
            temp_if_plot = 0;
            % [cY, cM] = fun_parellelTest(currentSession,temp_if_imageAverage,temp_if_plot,temp_if_max0_min1);
            [cY, cM] = fun_parellelTest(output_path,currentSession,temp_if_imageAverage,temp_finalOutputType,temp_if_plot);
            
            %% Co-plot PTB marker (in 2p frame space) and image results
            fig4 = figure('Name','Fig4, Behavior marker and image results','NumberTitle','off');
            set(gcf,'Position',[35 35 700 250]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            
            
            plot(cY, '-', 'LineWidth', 1, 'Color', [0.40 0.40 0.40]);
            hold on
            plot(cM, '-', 'LineWidth', 1, 'Color', [0.25 0.85 0.85]);
            hold on
            
            [min_cM,max_cM] = bounds(cM);
            for tempi=1:length(extractOutput_PTB.REWARD_downSampling.valid_marker_X_frameIndex)
                %plot(extractOutput_PTB.REWARD_downSampling.valid_marker_X_frameIndex(tempi)*ones(1,2),[min_cM max_cM], '-', 'LineWidth', 0.01, 'Color', [0.75 0.25 0.25]);
                plot(extractOutput_PTB.REWARD_downSampling.valid_marker_X_frameIndex(tempi)*ones(1,2),[max_cM/1.5 max_cM], '-', 'LineWidth', 0.01, 'Color', [0.75 0.25 0.25]);
                hold on
            end
            
            hl = legend('Raw','Motion correction','Reward',...
                'Location','northwest','fontsize',8);
            
            
            xlim([0 max(extractOutput_PTB.TRIALID_downSampling.valid_marker_X_frameIndex)]);
            ylim([0 1]);
            set(gca,'yticklabel',[])
            set(gca, 'FontSize', 11)
            set(gca,'box','off');% 取消右、上边框
            xlabel(sprintf('Frame number, T=%d',length(cY)), 'FontSize', 14, 'FontWeight', 'bold');
            % title(sprintf('Reward frame index (image), %d rewards',length(extractOutput_PTB.REWARD_downSampling.valid_marker_X_frameIndex)), 'FontSize', 12, 'FontWeight', 'bold');
            
        end
        
    end
    a = 1;
    
end

fprintf('Time of the stepA_marker_extract is %.1f secs.\n',toc(t0));
fprintf('num_FOV=%d.\n',num_FOV);

cd(targetPATH)
%% End