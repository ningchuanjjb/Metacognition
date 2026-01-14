close all

currentSessionIndex = 20;


temp_range = FOVAllCellRange_multiFOV(currentSessionIndex,1):FOVAllCellRange_multiFOV(currentSessionIndex,2);
tempSelectiveCellBoolIndex_length1_currentFOV = tempBoolIndex1(temp_range);
tempSelectiveCellBoolIndex_length2_currentFOV = tempBoolIndex2(temp_range);
tempSelectiveCellBoolIndex_length3_currentFOV = tempBoolIndex3(temp_range);
selectiveCellBoolIndex_rProb_glm_delayx_currentFOV = selectiveCellBoolIndex_rProb_glm_delayx_multiFOV(temp_range);
cellIndex_suite2p = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range);



tempSelectiveCellBoolIndex_length1_currentFOV;
tempSelectiveCellBoolIndex_length2_currentFOV;
tempSelectiveCellBoolIndex_length3_currentFOV;
selectiveCellBoolIndex_rProb_glm_delayx_currentFOV;



currentSession_multi = string;


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

currentSession_multi(1) = [];
num_FOV = length(currentSession_multi);



targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)


currentSession = currentSession_multi{currentSessionIndex};
fprintf('currentSession = %s.\n',currentSession);
    
    
output_shortPath = 'D:\twoPhotonData_motionCorrected';
temp_currentSession_path = [output_shortPath '\' currentSession];
temp_if_max0_min1 = 0;
output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);
path_plane = [output_path,'\plane0'];
    
temp_if_max0_min1 = 0;
%template_path = autoGetFileName_general('template*.tif', output_path,temp_if_max0_min1);
template_path = autoGetFileName_general('maxProjection*.tif', output_path,temp_if_max0_min1);
template = double(loadtiff(template_path));


fileName_Fall = 'Fall.mat';
fileName_iscell = 'iscell.npy';
fullFileName_Fall = [path_plane,'\',fileName_Fall];
fullFileName_iscell = [path_plane,'\',fileName_iscell];

iscell = readNPY(fullFileName_iscell);

s = load(fullFileName_Fall,'stat');
roi_stats_raw = s.stat;
temp_cellIndex = find(iscell(:,1)==1);
roi_stats = roi_stats_raw(temp_cellIndex);


%% roi anatomy info
fig = figure('Name',' ','NumberTitle','off');
set(gcf,'Position',[400 50 800 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

nexttile

%% plot whole FOV
temp_template = template / max(template,[],'all');

%temp_template2 = (temp_template.^0.7)*255;
temp_template2 = (temp_template.^0.55)*255;

image(temp_template2);
hold on
colormap(gray);


%% plot roi
temp_bound_offset_step = 3;

for tempi=1:4
% for tempi=1:3
% for tempi=4:-1:1
  
    if tempi == 1
        temp_selectiveCellBoolIndex = tempSelectiveCellBoolIndex_length1_currentFOV;
        temp_roiColor = [0.6350 0.0780 0.1840];
    elseif tempi == 2
        temp_selectiveCellBoolIndex = tempSelectiveCellBoolIndex_length2_currentFOV;        
        temp_roiColor = [0.4660 0.6740 0.1880];        
    elseif tempi == 3
        temp_selectiveCellBoolIndex = tempSelectiveCellBoolIndex_length3_currentFOV;        
        temp_roiColor = [0.3010 0.7450 0.9330];        
    elseif tempi == 4
        temp_selectiveCellBoolIndex = selectiveCellBoolIndex_rProb_glm_delayx_currentFOV;        
        temp_roiColor = [0.9290 0.6940 0.1250];        
    end
        
    temp_selectiveCellIndex = find(temp_selectiveCellBoolIndex==1);
    temp_selectiveCellBoolIndex_suite2p = cellIndex_suite2p(temp_selectiveCellBoolIndex);    
    temp_roiNum = sum(temp_selectiveCellBoolIndex);
    
    for tempj=1:temp_roiNum
    %for tempj=1:5
        
        
        temp_cellIndex = temp_selectiveCellIndex(tempj);
        
        %temp_cellIndex = 10;     
        
        if temp_cellIndex == 10
            a = 1;
        end
        
        %temp_cellIndex_suite2p = temp_selectiveCellBoolIndex_suite2p(tempj);
        temp_cellIndex_suite2p = cellIndex_suite2p(temp_cellIndex);
        
        temp_bound_offset = 0;
        
        if tempi >= 2
            if tempSelectiveCellBoolIndex_length1_currentFOV(temp_cellIndex) == 1
                temp_bound_offset = temp_bound_offset + temp_bound_offset_step;
            end
        end
        if tempi >= 3
            if tempSelectiveCellBoolIndex_length2_currentFOV(temp_cellIndex) == 1
                temp_bound_offset = temp_bound_offset + temp_bound_offset_step;
            end
        end
        if tempi >= 4
            if tempSelectiveCellBoolIndex_length3_currentFOV(temp_cellIndex) == 1
                temp_bound_offset = temp_bound_offset + temp_bound_offset_step;
            end
        end
        
        if tempi == 4
            temp_bound_offset = max(temp_bound_offset - 1,0);
        end
        
        a = 1;
        
                                
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
            plot(boundary(:,2)+1,boundary(:,1)+1,'color',temp_roiColor, 'LineWidth', 0.5);%1.5
            hold on
        end
        
        a = 1;
    end
    
end

axis equal;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'Visible','off');



a = 1;