% trials after correct
% Period 1: [ITI(end-1000ms:end), prefixation(1:700), fixation(1:400)]
% Period 2: reward(1:2800ms)

%% Initialization
clear
close all
% home;% To scroll down in command window

% currentSession = '113Recording_20230123A_Ding_Site09B_sameFOV0122';
% currentSession = '113Recording_20230423A_Ding_Site02';
% currentSession = '113Recording_20230424A_Ding_Site02_sameFOV0423';
% currentSession = '113Recording_20230502A_Ding_Site13';
% currentSession = '113Recording_20230503A_Ding_Site13_sameFOV0502';
% currentSession = '113Recording_20230504A_Ding_Site02';
% currentSession = '113Recording_20230508A_Ding_Site02';
% currentSession = '113Recording_20230509A_Ding_Site02_sameFOV0508';
currentSession = '113Recording_20230510A_Ding_Site05';
% currentSession = '113Recording_20230510B_Ding_Site05';


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)

% rawData_path = 'D:\twoPhotonRawData\ToBeProcessed';
rawData_path = 'R:\twoPhotonRawData\ToBeProcessed';
output_shortPath = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected';

temp_currentSession_path = [output_shortPath '\' currentSession];
currentSession_path = [rawData_path '\' currentSession];

PTB_name = '.mat';
PTB_fullName = autoGetMarkerFileName(['2*',PTB_name], currentSession_path);

temp_if_max0_min1 = 0;
output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
path_plane = [output_path,'\plane0'];

load([output_path,'\','markerParse_trialLevel.mat'],'markerParse_trialLevel');


fileName_Fall = 'Fall.mat';
fileName_iscell = 'iscell.npy';
fullFileName_Fall = [path_plane,'\',fileName_Fall];
fullFileName_iscell = [path_plane,'\',fileName_iscell];
load(fullFileName_Fall,'F_dff');
%load(fullFileName_Fall);
load(fullFileName_Fall,'F_dff_raw');
iscell = readNPY(fullFileName_iscell);

if isempty(markerParse_trialLevel{end}.currentTrialTotalMarkers_frameIndex) == 1
    markerParse_trialLevel(end) = [];
end

F_dff_notcell = F_dff_raw(iscell(:,1)==0,:);
clear F_dff_raw

markerParse_trialLevel;
iscell;
F_dff;
F_dff_notcell;

load(PTB_fullName,'basic_para','trial_para');

%% Further computation

get_F_dff_highCorrClusterName = 'get_F_dff_mean_highCorrCluster_notcell';
get_F_dff_highCorrClusterName_v = autoGetFunName_myScripts(get_F_dff_highCorrClusterName, [targetPATH '\functions']);
fun_F_dff_highCorrCluster_v = str2func(get_F_dff_highCorrClusterName_v);

temp_binSize = 6;%6
corr_threshold = 0.5;%0.5

[F_dff_median_highCorrCluster_notcell,F_dff_mean_highCorrCluster_notcell,highCorrClusterIndex_notcell_suite2p] = ...
    fun_F_dff_highCorrCluster_v(F_dff_notcell,temp_binSize,corr_threshold,iscell);



correctTrialCount = 0;
for trialIndex_marker=1:(length(markerParse_trialLevel)-1)
    trialIndex_PTB = markerParse_trialLevel{trialIndex_marker}.trial_index;
    
    if trial_para.isCorrect(trialIndex_PTB)==1
        correctTrialCount = correctTrialCount + 1;
    end
end

prePeriod2_duration = 60;
period2_duration = (2800/1000)*30;% 2800ms --> x framesF_dff_afterReward_period2 = zeros(size(F_dff,1),correctTrialCount,totalPeriod_duration,'single');

totalPeriod_duration = prePeriod2_duration + period2_duration;

F_dff_afterReward_period2 = zeros(size(F_dff,1),correctTrialCount,totalPeriod_duration,'single');
F_dff_mean_highCorrCluster_afterReward = zeros(correctTrialCount,totalPeriod_duration,'single');
F_dff_median_highCorrCluster_afterReward = zeros(correctTrialCount,totalPeriod_duration,'single');


rewardFrameIndex = zeros(1, correctTrialCount);
correctTrialCount2 = 0;
for trialIndex_marker=1:(length(markerParse_trialLevel)-1)
    trialIndex_PTB = markerParse_trialLevel{trialIndex_marker}.trial_index;
    
    if trial_para.isCorrect(trialIndex_PTB)==1
        correctTrialCount2 = correctTrialCount2 + 1;
        
        temp_rewardFrameIndex = markerParse_trialLevel{trialIndex_marker}.currentTrialTotalMarkers_frameIndex(end);
        a = 1;                %#ok<*NASGU>
        %F_dff_afterReward_period2(:,correctTrialCount2,:) = F_dff(:,temp_rewardFrameIndex:(temp_rewardFrameIndex+totalPeriod_duration-1));
        F_dff_afterReward_period2(:,correctTrialCount2,:) = F_dff(:,temp_rewardFrameIndex-prePeriod2_duration:(temp_rewardFrameIndex+period2_duration-1));        
        a = 1;
        rewardFrameIndex(correctTrialCount2) = temp_rewardFrameIndex;
        
        %F_dff_mean_highCorrCluster_afterReward(correctTrialCount2,:) = F_dff_mean_highCorrCluster_notcell(temp_rewardFrameIndex:(temp_rewardFrameIndex+totalPeriod_duration-1));
        %F_dff_median_highCorrCluster_afterReward(correctTrialCount2,:) = F_dff_median_highCorrCluster_notcell(temp_rewardFrameIndex:(temp_rewardFrameIndex+totalPeriod_duration-1));        
        
        F_dff_mean_highCorrCluster_afterReward(correctTrialCount2,:) = F_dff_mean_highCorrCluster_notcell(temp_rewardFrameIndex-prePeriod2_duration:(temp_rewardFrameIndex+period2_duration-1));
        F_dff_median_highCorrCluster_afterReward(correctTrialCount2,:) = F_dff_median_highCorrCluster_notcell(temp_rewardFrameIndex-prePeriod2_duration:(temp_rewardFrameIndex+period2_duration-1));        
                
    end
end
a = 1;
correctTrialCount;
F_dff_afterReward = F_dff_afterReward_period2;
% F_dff_afterReward_roiMean = squeeze(mean(F_dff_afterReward,1));
% F_dff_afterReward_roiMedian = squeeze(median(F_dff_afterReward,1));
F_dff_afterReward_roiMean = squeeze(mean(F_dff_afterReward(1:50,:,:),1));
F_dff_afterReward_roiMedian = squeeze(median(F_dff_afterReward(1:50,:,:),1));

roiNum = size(F_dff,1);

corrWindow = 7;%5-->7

filterWindow = 2; 
b = (1/filterWindow)*ones(1,filterWindow);
a = 1;

vibrate_trialXframe = zeros(correctTrialCount,totalPeriod_duration);
vibrate_trialXframe2 = zeros(correctTrialCount,totalPeriod_duration);

% tempi = 10;
for tempi=1:correctTrialCount
    temp_rewardFrameIndex = rewardFrameIndex(tempi);
    temp_F = squeeze(F_dff_afterReward(:,tempi,:))';
    temp_F_smoothed = filter(b,a,temp_F);
    
    temp_F_roiMedian = F_dff_afterReward_roiMedian(tempi,:)';
    temp_F_roiMean = F_dff_afterReward_roiMean(tempi,:)';
    
        
    temp_F_highCorrCluster = F_dff_median_highCorrCluster_afterReward(tempi,:)';

    
    rho = zeros(totalPeriod_duration,roiNum);
    for tempj=1:totalPeriod_duration-corrWindow+1
        %temp2_F = temp_F(tempj:tempj+corrWindow-1,:);
        temp2_F = temp_F_smoothed(tempj:tempj+corrWindow-1,:);
        temp2_F_roiMedian = temp_F_roiMedian(tempj:tempj+corrWindow-1);
        temp2_F_roiMean = temp_F_roiMean(tempj:tempj+corrWindow-1);
        temp2_F_highCorrCluster = temp_F_highCorrCluster(tempj:tempj+corrWindow-1);
        
        
        %rho(tempj,:) = corr(temp2_F,temp2_F_roiMedian);
        %rho(:,tempj) = corr(temp2_F,temp2_F_roiMean);               
        rho(tempj,:) = corr(temp2_F,temp2_F_highCorrCluster);
    end
    
    %boolIndex1 = rho<-0.5;
    %boolIndex1 = rho<-0.8;    
    %boolIndex1 = rho>0.8;
    %boolIndex1 = rho>0.5;  
    boolIndex1 = abs(rho)>0.85;%0.8


    %boolIndex2 = temp_F'<0;        
    boolIndex2 = temp_F_smoothed<0;    
    
    
    temp_F_negative_roiMedian = zeros(totalPeriod_duration,1);
    temp_F_negative_CI = zeros(totalPeriod_duration,1);
    for tempj=1:totalPeriod_duration
        current_F_negative = temp_F_smoothed(tempj,boolIndex2(tempj,:));
        temp_F_negative_roiMedian(tempj) = median(current_F_negative);
        
        %temp_F_negative_CI(tempj) = 0-3*std(current_F_negative);
        temp_F_negative_CI(tempj) = 0-2.5*std(current_F_negative);
        %temp_F_negative_CI(tempj) = 0-2*std(current_F_negative);
    end
    temp_F_negative_roiMedian;
    temp_F_negative_CI;
    
    boolIndex3 = temp_F_smoothed<temp_F_negative_CI;
    
    
    %bb=sum(boolIndex1 & boolIndex2,2);
    %bb=sum(boolIndex1 & boolIndex3,2);
    %bb=sum(boolIndex1,1);
    %bb=sum(boolIndex3,2);
    
    vibrate_trialXframe(tempi,:) = sum(boolIndex1,2);
    vibrate_trialXframe2(tempi,:) = sum(boolIndex3,2);
    %vibrate_trialXframe(tempi,:) = sum(boolIndex1 & boolIndex3,2);
    %vibrate_trialXframe(tempi,:) = sum(boolIndex1 | boolIndex3,2);
    
    temp_rewardFrameIndex;
    rewardFrameIndex;
    a = 1;
end
a = 1;

vibrate_trialSum = sum(vibrate_trialXframe,2);
[M, I] = max(vibrate_trialSum); %#ok<*ASGLU>
vibrate_frameMean = mean(vibrate_trialXframe,1);
vibrate_frameprctileA = prctile(vibrate_trialXframe,95,'Method','approximate');
vibrate_frameprctileB = prctile(vibrate_trialXframe,75,'Method','approximate');

vibrate_trialSum2 = sum(vibrate_trialXframe2,2);
[M, I] = max(vibrate_trialSum2);
vibrate_frameMean2 = mean(vibrate_trialXframe2,1);
vibrate_frameprctileA2 = prctile(vibrate_trialXframe2,95,'Method','approximate');
vibrate_frameprctileB2 = prctile(vibrate_trialXframe2,75,'Method','approximate');

%% Plot
% fig1 = figure('Name','Fig1, Two photon marker','NumberTitle','off');
fig1 = figure(1);
set(gcf,'Position',[400 100 600 750]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

nexttile

plot(mean(F_dff_median_highCorrCluster_afterReward,1)');
hold on
plot(prctile(F_dff_median_highCorrCluster_afterReward,75,1,'Method','approximate')');
hold on
plot(prctile(F_dff_median_highCorrCluster_afterReward,95,1,'Method','approximate')');
hold on
tempy = ylim;
min_y = tempy(1);
max_y = tempy(2);
plot([1 1]*prePeriod2_duration,[min_y max_y],'k','LineWidth',1);

hl = legend('Mean','75 percentile','95 percentile',...
    'Location','northeast','fontsize',8);

set(gca,'XTick',0:10:totalPeriod_duration);

title(sprintf('Ref roi cluster dff'), 'FontSize', 12, 'FontWeight', 'bold');


nexttile

plot(vibrate_frameMean');
hold on
plot(vibrate_frameprctileB');
hold on
plot(vibrate_frameprctileA');
hold on
tempy = ylim;
min_y = tempy(1);
max_y = tempy(2);
plot([1 1]*prePeriod2_duration,[min_y max_y],'k','LineWidth',1);

hl = legend('Mean','75 percentile','95 percentile',...
    'Location','northeast','fontsize',8);
        
set(gca,'XTick',0:10:totalPeriod_duration);

title(sprintf('High correlation ROI number with ref'), 'FontSize', 12, 'FontWeight', 'bold');

nexttile

plot(vibrate_frameMean2');
hold on
plot(vibrate_frameprctileB2');
hold on
plot(vibrate_frameprctileA2');
hold on
tempy = ylim;
min_y = tempy(1);
max_y = tempy(2);
plot([1 1]*prePeriod2_duration,[min_y max_y],'k','LineWidth',1);

hl = legend('Mean','75 percentile','95 percentile',...
    'Location','northeast','fontsize',8);

set(gca,'XTick',0:10:totalPeriod_duration);

xlabel(sprintf('Frame number, after reward'), 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('dff negative ROI number'), 'FontSize', 12, 'FontWeight', 'bold');
a = 1;

%% End