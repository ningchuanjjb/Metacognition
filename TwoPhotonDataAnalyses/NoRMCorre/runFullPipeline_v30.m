function runFullPipeline_v30(currentSession_external)
% clear
% close all

if(~exist('currentSession_external','var'))
    clear
%     currentSession_external = '113Recording_20230503A_Ding_Site13_sameFOV0502';
%     currentSession_external = '113Recording_20230509A_Ding_Site02';
%     currentSession_external = '113Recording_20230511A_Ding_Site05';
%     currentSession_external = '113Recording_20230522A_Ding_Site05B';
%     currentSession_external = '113Recording_20230525A_Ding_Site09B';
%     currentSession_external = '113Recording_20230528A_Ding_Site02B_sameFOV0527';
%     currentSession_external = '113Recording_20230530A_Ding_Site05C';
%     currentSession_external = '113Recording_20230620A_Ding_Site05E';
    
%     currentSession_external = '113Recording_20231205A_Zelku_noImage';
%     currentSession_external = '329Recording_20231227A_YRK';

    %currentSession_external = '113Recording_20240111A_Zelku_Site09A';
    %currentSession_external = '113Recording_20240112A_Zelku_Site06A';
    %currentSession_external = '113Recording_20240115A_Zelku_Site06A';
    %currentSession_external = '113Recording_20240117A_Zelku_Site14A';
    %currentSession_external = '113Recording_20240118A_Zelku_Site18A';
    %currentSession_external = '113Recording_20240119A_Zelku_Site17A';
    %currentSession_external = '113Recording_20240122A_Zelku_Site09B';
    %currentSession_external = '113Recording_20240123A_Zelku_Site09B_sameFOV0122';
    %currentSession_external = '113Recording_20240124A_Zelku_Site06B';
    %currentSession_external = '113Recording_20240126A_Zelku_Site06B_sameFOV0124';
    %currentSession_external = '113Recording_20240129A_Zelku_Site07A';
    %currentSession_external = '113Recording_20240131A_Zelku_Site07A_sameFOV0129';
    %currentSession_external = '113Recording_20240202A_Zelku_Site06XA';
    %currentSession_external = '113Recording_20240203A_Zelku_Site06XA_sameFOV0202';
    %currentSession_external = '113Recording_20240207A_Zelku_Site05A';
    %currentSession_external = '113Recording_20240208A_Zelku_Site05A_sameFOV0207';
    %currentSession_external = '113Recording_20240210A_Zelku_Site10A';
    %currentSession_external = '113Recording_20240211A_Zelku_Site10A_sameFOV0210';
    %currentSession_external = '113Recording_20240216A_Zelku_Site09C';
    %currentSession_external = '113Recording_20240218A_Zelku_Site09C_sameFOV0216';
    %currentSession_external = '113Recording_20240220A_Zelku_Site06XB';
    %currentSession_external = '113Recording_20240221A_Zelku_Site06XB_sameFOV0220';
    %currentSession_external = '113Recording_20240226A_Zelku_Site10B';
    %currentSession_external = '113Recording_20240227A_Zelku_Site10B_sameFOV0226';
    %currentSession_external = '113Recording_20240229A_Zelku_Site06C';
    %currentSession_external = '113Recording_20240301A_Zelku_Site06C_sameFOV0229';
    %currentSession_external = '113Recording_20240304A_Zelku_Site09D';
    %currentSession_external = '113Recording_20240305A_Zelku_Site09D_sameFOV0304';
    %currentSession_external = '113Recording_20240307A_Zelku_Site10C';
    %currentSession_external = '113Recording_20240308A_Zelku_Site10C_sameFOV0307';
    %currentSession_external = '113Recording_20240312A_Zelku_Site06RA';
    %currentSession_external = '113Recording_20240315A_Zelku_Site06RA_sameFOV0312';
    %currentSession_external = '113Recording_20240319A_Zelku_Site09E';
    %currentSession_external = '113Recording_20240320A_Zelku_Site09E_sameFOV0319';
    %currentSession_external = '113Recording_20240322A_Zelku_Site07B';
    %currentSession_external = '113Recording_20240323A_Zelku_Site07B_sameFOV0322';
    %currentSession_external = '113Recording_20240329A_Zelku_Site05B';
    %currentSession_external = '113Recording_20240330A_Zelku_Site05B_sameFOV0329';
    %currentSession_external = '113Recording_20240402A_Zelku_Site14B';
    %currentSession_external = '113Recording_20240403A_Zelku_Site14B_sameFOV0402';
    %currentSession_external = '113Recording_20240410A_Zelku_Site17B';
    currentSession_external = '113Recording_20240411A_Zelku_Site17B_sameFOV0410';
    
end
home

if_delete_dataBin = 1;%1

if_sessionA0_sessionB1 = 0;

if_NoRMCorre = 1;
if_runScript = 1;
% if_modifyCellProb = 1;
%if_mappingTwoSessionFall = 0;
%if_categorySomaDendrite = 0;
if_dff = 1;
if_modifyCellProb_dff = 1;
if_getMaxProj = 1;
if_guiScript = 1;

if if_sessionA0_sessionB1 == 1
    % for sessionB
    if_mappingTwoSessionFall = 1;
    if_dff_sessionA = 0;
    if_categorySomaDendrite_sessionA = 1;
    if_categorySomaDendrite = 1;
    
elseif if_sessionA0_sessionB1 == 0
    % for sessionA
    if_mappingTwoSessionFall = 0; % always 0
    if_dff_sessionA = 0;
    if_categorySomaDendrite_sessionA = 0; % always 0
    if_categorySomaDendrite = 0; % always 0
end

% Params of dff
options_dff = struct;
options_dff.if_dff = if_dff;
options_dff.if_dff_sessionA = if_dff_sessionA;

options_dff.df_window = 1000;
options_dff.if_save_dff = 1;
options_dff.if_plot_dff = 0;%1
options_dff.roi_ID_forPlot = 15;%from 0

options_category = struct;
options_category.if_categorySomaDendrite_sessionA = if_categorySomaDendrite_sessionA;

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)

%% Run motion correction

NoRMCorreName = 'myNoRMCorre';

NoRMCorreName_v = autoGetFunName(NoRMCorreName, targetPATH);
fprintf('Now runing %s to process %s.\n',NoRMCorreName_v,currentSession_external);
pause(0.2);%1
fun_NoRMCorre = str2func(NoRMCorreName_v);

t0 = tic;

[output_shortPath,output_shortPath_temporary,currentSession,existingTemplateSession,if_useExistingTemplate] = ...
    fun_NoRMCorre(if_NoRMCorre,currentSession_external);

temp_currentSession_path = [output_shortPath '\' currentSession];
if ~exist(temp_currentSession_path, 'dir')
    mkdir(temp_currentSession_path);
end
temp_if_max0_min1 = 0;
output_path = autoGetFileName_general('Result', temp_currentSession_path,temp_if_max0_min1);

[~,parts_fname,~] = fileparts(output_path);
Result_str = parts_fname;
output_path_temporary = [output_shortPath_temporary '\' currentSession '\' Result_str];

if if_useExistingTemplate == 1
    existingTemplateSession_path = [output_shortPath '\' existingTemplateSession];
    temp_if_max0_min1 = 0;
    Result_existingTemplateSession_path = autoGetFileName_general('Result', existingTemplateSession_path,temp_if_max0_min1);
    
    file_list=dir([Result_existingTemplateSession_path,'\template*']);
    existingTemplate_fileName = [Result_existingTemplateSession_path,'\',file_list(end).name];
    
    [~,parts_existingTemplate_fileName,~] = fileparts(existingTemplate_fileName);
    fprintf('Now using existing template from %s.\n',parts_existingTemplate_fileName);
end

% home

targetPATH2 = temp_currentSession_path;
% targetPATH2 = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230115A_Ding_Site06_sameFOV0112'
% targetPATH2 = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230122A_Ding_Site09B'
% targetPATH2 = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122'

cd(targetPATH2)


%% Run suite2p script
% Use ! to run cmd, use & to run cmd in a new cmd window
% !runScript_ForSuite2p.cmd
% !runScript_ForSuite2p_pipeline.cmd
path_my2pScripts = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts'; %#ok<*UNRCH>
% cd(path_my2pScripts)

if if_runScript == 1
    %if if_NoRMCorre == 1
        delete(gcp); %#ok<*UNRCH>
    %end
    
    ! C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\forSuite2p\runScript_ForSuite2p_pipeline.cmd &
    %! C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\forSuite2p\runScript_ForSuite2p_pipeline.cmd
    
    m_runScriptState = memmapfile('m_runScriptState.dat', 'Writable', true, 'Format', 'double');
    m_runScriptState.Data = 0;
    while true
        pause(0.5);
        if m_runScriptState.Data == 1
            m_runScriptState.Data = 0;
            fprintf('runScript_ForSuite2p_pipeline is over!\n');
            break
        end
    end    
    gcp;
end

% modifyCellProbName = 'modifyCellProb';
modifyCellProbName = 'modifyCellProb_dff';
modifyCellProb_v = autoGetFunName(modifyCellProbName, [path_my2pScripts,'\','functions']);
fun_modifyCellProb = str2func(modifyCellProb_v);

%% To modify cellProb from suite2p result
path_plane = [output_path,'\plane0'];
% path_plane = ['C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230115A_Ding_Site06_sameFOV0112\Result20230201T191002\plane0'];
% path_plane = ['C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230122A_Ding_Site09B\Result20230202T101726\plane0'];
% path_plane = ['C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230122A_Ding_Site09B\Result20230202T204254\plane0'];
% path_plane = ['C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\Result20230202T231115\plane0']; %#ok<*NBRAK>

if if_runScript == 1
    % To conver Fall.mat to v7.3 mat version
    
    fileName_Fall = 'Fall.mat';
    fileName_Fall_new = 'Fall_new.mat';
    fullFileName_Fall = [path_plane,'\',fileName_Fall];

    fullFileName_Fall_new = [path_plane,'\',fileName_Fall_new];
    
    tempS = load(fullFileName_Fall);
    % Save mat by extracting every variable in a struct
    save(fullFileName_Fall_new,'-struct','tempS','-v7.3','-nocompression');
    recycle('off');
    delete(fullFileName_Fall);
    recycle('on');
    movefile(fullFileName_Fall_new,fullFileName_Fall);
    clear tempS
    
    fileName_iscell = 'iscell.npy';
    fileName_iscell_raw = 'iscell_raw.npy';
    fullFileName_iscell = [path_plane,'\',fileName_iscell];
    fullFileName_iscell_raw = [path_plane,'\',fileName_iscell_raw];
    
    copyfile(fullFileName_iscell,fullFileName_iscell_raw);
    
    if if_delete_dataBin == 1
        recycle('off');
        fileName_dataBin = 'data.bin';
        fullFileName_dataBin = [output_path_temporary,'\',fileName_dataBin];
        delete(fullFileName_dataBin);
        fileName_binnedData = 'binnedData.h5';
        fullFileName_binnedData = [output_path_temporary,'\',fileName_binnedData];
        delete(fullFileName_binnedData);      
        
        rmdir(output_path_temporary);
        
        recycle('on');
    end
end

calculate_dffName = 'calculate_dff';
calculate_dff_v = autoGetFunName(calculate_dffName, [path_my2pScripts,'\','functions']);
fun_calculate_dff = str2func(calculate_dff_v);
options_dff.fun = fun_calculate_dff;

%% To calculate F_dff and F0 from F
if if_dff == 1
    fprintf('Now calculating dff for current session.\n');
    fun_calculate_dff(path_plane,options_dff);
end

% if if_modifyCellProb == 1
if if_modifyCellProb_dff == 1
    fun_modifyCellProb(path_plane);
end

mappingTwoSessionFallName = 'mappingTwoSessionFall';
mappingTwoSessionFall_v = autoGetFunName(mappingTwoSessionFallName, [path_my2pScripts,'\','functions']);
fun_mappingTwoSessionFall = str2func(mappingTwoSessionFall_v);

categorySomaDendriteName = 'categorySomaDendrite';
categorySomaDendrite_v = autoGetFunName(categorySomaDendriteName, [path_my2pScripts,'\','functions']);
fun_categorySomaDendrite = str2func(categorySomaDendrite_v);
options_category.fun = fun_categorySomaDendrite;

%% To merge two sessions' results from suite2p
if if_mappingTwoSessionFall == 1 && if_useExistingTemplate == 1
    path_plane_sessionA = [Result_existingTemplateSession_path,'\plane0'];
    path_plane_sessionB = path_plane;
    %fun_mappingTwoSessionFall(path_plane_sessionA,path_plane_sessionB);
    fun_mappingTwoSessionFall(path_plane_sessionA,path_plane_sessionB,options_dff,options_category);
end



%% To category soma and dendrite
if if_categorySomaDendrite == 1
    fun_categorySomaDendrite(path_plane);
end


getMaxProj_Name = 'getMaxProj';
getMaxProj_v = autoGetFunName(getMaxProj_Name, [targetPATH,'\','functions']);
fun_getMaxProj = str2func(getMaxProj_v);

%% Get max projection
if if_getMaxProj== 1
    fun_getMaxProj(output_path);
end

%% Run suite2p gui
% !guiScript_ForSuite2p.cmd &
if if_guiScript == 1
    if if_mappingTwoSessionFall == 1 && if_useExistingTemplate == 1
        cd(path_plane_sessionA)
        ! C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\forSuite2p\guiScript_ForSuite2p.cmd &
    end
    
    cd(path_plane)
    ! C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\forSuite2p\guiScript_ForSuite2p.cmd &
end

% home
% fprintf('Time of the whole pipeline is %.1f secs.\n',toc(t0));
fprintf('Time of the current session pipeline is %.1f secs.\n',toc(t0));
pause(0.2);%1

%% END
% cd(targetPATH)
cd(path_my2pScripts)