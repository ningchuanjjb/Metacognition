% Chuan's 30th script (20251214)
% This script: patch of 'stepF3_test_decodingDataSimplified_AB.m'
%% Initialization
% clear
% clear decodingDataSimplified_A decodingDataSimplified_B
close all

if_monkey_D0_Z1 = 1;

if_save = 0;



% % % 0122A-0123A (Z)
% currentSessionIndex_A = 7;
% currentSessionIndex_B = 8;

% % 0124A-0126A (Z)
% currentSessionIndex_A = 9;
% currentSessionIndex_B = 10;

% % 0129A-0131A (Z)
% currentSessionIndex_A = 11;
% currentSessionIndex_B = 12;

% % 0202A-0203A (Z)
% currentSessionIndex_A = 13;
% currentSessionIndex_B = 14;

% % 0207A-0208A (Z)
% currentSessionIndex_A = 15;
% currentSessionIndex_B = 16;

% % 0210A-0211A (Z)
% currentSessionIndex_A = 17;
% currentSessionIndex_B = 18;

% % 0216A-0218A (Z)
% currentSessionIndex_A = 19;
% currentSessionIndex_B = 20;

% % 0220A-0221A (Z)
% currentSessionIndex_A = 21;
% currentSessionIndex_B = 22;

% % 0226A-0227A (Z)
% currentSessionIndex_A = 23;
% currentSessionIndex_B = 24;

% % 0229A-0301A (Z)
% currentSessionIndex_A = 25;
% currentSessionIndex_B = 26;

% % 0304A-0305A (Z)
% currentSessionIndex_A = 27;
% currentSessionIndex_B = 28;

% % 0307A-0308A (Z)
% currentSessionIndex_A = 29;
% currentSessionIndex_B = 30;

% % 0312A-0315A (Z)
currentSessionIndex_A = 31;
currentSessionIndex_B = 32;

% % 0319A-0320A (Z)
% currentSessionIndex_A = 33;
% currentSessionIndex_B = 34;

% % 0322A-0323A (Z)
% currentSessionIndex_A = 35;
% currentSessionIndex_B = 36;

% % 0329A-0330A (Z)
% currentSessionIndex_A = 37;
% currentSessionIndex_B = 38;

% % 0402A-0403A (Z)
% currentSessionIndex_A = 39;
% currentSessionIndex_B = 40;

% % 0410A-0411A (Z)
% currentSessionIndex_A = 41;
% currentSessionIndex_B = 42;


% % 0527A-0528A
% currentSessionIndex_A = 19;%19
% currentSessionIndex_B = 20;%20

% % 0512A-0513A
% currentSessionIndex_A = 11;
% currentSessionIndex_B = 12;

% % 0530A-0531A
% currentSessionIndex_A = 21;
% currentSessionIndex_B = 22;


% % 0426A-0427A
% currentSessionIndex_A = 1;
% currentSessionIndex_B = 2;

% % 0502A-0503A
% currentSessionIndex_A = 3;
% currentSessionIndex_B = 4;

% % 0504A-0509A
% currentSessionIndex_A = 5;
% currentSessionIndex_B = 7;

% % 0508A-0509A
% currentSessionIndex_A = 6;
% currentSessionIndex_B = 7;

% % 0515A-0516A
% currentSessionIndex_A = 13;
% currentSessionIndex_B = 14;

% % 0525A-0526A
% currentSessionIndex_A = 17;
% currentSessionIndex_B = 18;

% % 0601A-0602A
% currentSessionIndex_A = 23;
% currentSessionIndex_B = 24;

% % 0604A-0605A
% currentSessionIndex_A = 25;
% currentSessionIndex_B = 26;

% % 0614A-0615A
% currentSessionIndex_A = 28;
% currentSessionIndex_B = 29;

% % 0510A-0511A
% currentSessionIndex_A = 8;
% currentSessionIndex_B = 10;

% % 0510B-0511A
% currentSessionIndex_A = 9;
% currentSessionIndex_B = 10;




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


% currentSession_A = '113Recording_20230527A_Ding_Site02B';
% currentSession_B = '113Recording_20230528A_Ding_Site02B_sameFOV0527';

currentSession_A = currentSession_multi{currentSessionIndex_A};
currentSession_B = currentSession_multi{currentSessionIndex_B};

fprintf('currentSession_A = %s.\n',currentSession_A);
fprintf('currentSession_B = %s.\n',currentSession_B);


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)

output_shortPath = 'D:\twoPhotonData_motionCorrected';
temp_if_max0_min1 = 0;

%% Load
currentSession = currentSession_A;
temp_currentSession_path = [output_shortPath '\' currentSession];
input_path_A = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
if exist('decodingDataSimplified_A','var') == 0
    temp_load = load([input_path_A,'\','decodingDataSimplified.mat'],'decodingDataSimplified');
    decodingDataSimplified_A = temp_load.decodingDataSimplified;
end

currentSession = currentSession_B;
temp_currentSession_path = [output_shortPath '\' currentSession];
input_path_B = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
if exist('decodingDataSimplified_B','var') == 0
    temp_load = load([input_path_B,'\','decodingDataSimplified.mat'],'decodingDataSimplified');
    decodingDataSimplified_B = temp_load.decodingDataSimplified;
end

temp_fileName = [currentSession_A(14:22),'_and_',currentSession_B(14:22)];
input_path_AB = [output_shortPath '\twoSessions\raw\' temp_fileName];
if exist('decodingDataSimplified_AB_raw','var') == 0
    temp_load = load(input_path_AB,'decodingDataSimplified_AB');
    decodingDataSimplified_AB_raw = temp_load.decodingDataSimplified_AB;
end


%%
decodingDataSimplified_A;
decodingDataSimplified_B;
decodingDataSimplified_AB_raw;

tempMappingCellIndex = decodingDataSimplified_AB_raw.extraForMerged.tempMappingCellIndex;

roiNum_AB_shared = decodingDataSimplified_AB_raw.extraForMerged.roiNum_AB;
trialNum_A = decodingDataSimplified_AB_raw.extraForMerged.trialNum_A;
trialNum_B = decodingDataSimplified_AB_raw.extraForMerged.trialNum_B;

% F_dff_decisionPeriodB
temp1 = decodingDataSimplified_A.F_dff_decisionPeriodB(tempMappingCellIndex(:,1),:,:);
temp2 = decodingDataSimplified_B.F_dff_decisionPeriodB(tempMappingCellIndex(:,2),:,:);
temp3 = zeros(roiNum_AB_shared,trialNum_A+trialNum_B,size(decodingDataSimplified_A.F_dff_decisionPeriodB,3));
temp3(:,1:trialNum_A,:) = temp1;
temp3(:,trialNum_A+(1:trialNum_B),:) = temp2;
F_dff_decisionPeriodB = temp3;

decodingDataSimplified_AB = decodingDataSimplified_AB_raw;
decodingDataSimplified_AB.F_dff_decisionPeriodB = F_dff_decisionPeriodB;
decodingDataSimplified_AB.decisionPeriodB_interval = decodingDataSimplified_A.decisionPeriodB_interval;

%% Save
if if_save == 1
    temp_fileName = [currentSession_A(14:22),'_and_',currentSession_B(14:22)];
    output_path = [output_shortPath '\twoSessions\' temp_fileName];
    save(output_path,'decodingDataSimplified_AB');
end

%% End