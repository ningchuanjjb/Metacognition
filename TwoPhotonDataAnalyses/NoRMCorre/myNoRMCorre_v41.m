function [output_shortPath,output_shortPath_temporary,currentSession,existingTemplateSession,if_useExistingTemplate] = ...
    myNoRMCorre_v41(if_NoRMCorre,currentSession_external)
%% Initialization
% clear
close all
% home;% To scroll down in command window

if(~exist('if_NoRMCorre','var'))
    clear    
    if_NoRMCorre = 1;
end

t0 = tic;

% g = gpuDevice(1);
% reset(g);
% wait(g)

if if_NoRMCorre == 1
    gcp;% To open parallel pool
end
if false
    delete(gcp); %#ok<*UNRCH>
    gcp;
end

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)

% currentSession = '113Recording_20230122A_Ding_Site09B';
% currentSession = '113Recording_20230123A_Ding_Site09B_sameFOV0122';
% currentSession = '113Recording_20230423A_Ding_Site02';
% currentSession = '113Recording_20230424A_Ding_Site02_sameFOV0423';

% currentSession = '113Recording_20230426A_Ding_Site16';
% currentSession = '113Recording_20230427A_Ding_Site16_sameFOV0426';
% currentSession = '113Recording_20230502A_Ding_Site13';
% currentSession = '113Recording_20230503A_Ding_Site13_sameFOV0502';
% currentSession = '113Recording_20230504A_Ding_Site02';
% currentSession = '113Recording_20230508A_Ding_Site02_sameFOV0509';
% currentSession = '113Recording_20230509A_Ding_Site02';
% currentSession = '113Recording_20230510A_Ding_Site05_sameFOV0511';
% currentSession = '113Recording_20230510B_Ding_Site05_sameFOV0511';
% currentSession = '113Recording_20230511A_Ding_Site05';
% currentSession = '113Recording_20230512A_Ding_Site09';
% currentSession = '113Recording_20230513A_Ding_Site09_sameFOV0512';
% currentSession = '113Recording_20230515A_Ding_Site24_sameFOV0516';
% currentSession = '113Recording_20230516A_Ding_Site24';
% currentSession = '113Recording_20230517A_Ding_Site16B';
% currentSession = '113Recording_20230522A_Ding_Site05B';
% currentSession = '113Recording_20230525A_Ding_Site09B';
% currentSession = '113Recording_20230526A_Ding_Site09B_sameFOV0525';
% currentSession = '113Recording_20230527A_Ding_Site02B';
% currentSession = '113Recording_20230528A_Ding_Site02B_sameFOV0527';
% currentSession = '113Recording_20230530A_Ding_Site05C';
% currentSession = '113Recording_20230531A_Ding_Site05C_sameFOV0530';
% currentSession = '113Recording_20230601A_Ding_Site13B';
% currentSession = '113Recording_20230602A_Ding_Site13B_sameFOV0601';
% currentSession = '113Recording_20230604A_Ding_Site07';
% currentSession = '113Recording_20230605A_Ding_Site07_sameFOV0604';
% currentSession = 'test0608';
% currentSession = '113Recording_20230612A_Ding_Site14';
% currentSession = '113Recording_20230614A_Ding_Site05D';
% currentSession = '113Recording_20230615A_Ding_Site05D_sameFOV0614';
% currentSession = '113Recording_20230619A_Ding_Site02C';
% currentSession = '113Recording_20230620A_Ding_Site05E';

currentSession = currentSession_external;

% existingTemplateSession = '113Recording_20230122A_Ding_Site09B';
% existingTemplateSession = '113Recording_20230423A_Ding_Site02';
% existingTemplateSession = '113Recording_20230426A_Ding_Site16';
% existingTemplateSession = '113Recording_20230502A_Ding_Site13';
% existingTemplateSession = '113Recording_20230504A_Ding_Site02';
% existingTemplateSession = '113Recording_20230509A_Ding_Site02';
% existingTemplateSession = '113Recording_20230511A_Ding_Site05';
% existingTemplateSession = '113Recording_20230512A_Ding_Site09';
existingTemplateSession = '113Recording_20230516A_Ding_Site24';


% rawData_path = 'D:\twoPhotonRawData\ToBeProcessed';
rawData_path = 'R:\twoPhotonRawData\ToBeProcessed';
currentSession_path = [rawData_path '\' currentSession];
% output_shortPath = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected';
output_shortPath = 'D:\twoPhotonData_motionCorrected';
output_shortPath_temporary = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected_temporary';

datestr_now30 = datestr(now,30);
output_path = [output_shortPath '\' currentSession '\Result' datestr_now30];
output_path_temporary = [output_shortPath_temporary '\' currentSession '\Result' datestr_now30];
% if ~exist(output_path, 'dir')
if ~exist(output_path, 'dir') && if_NoRMCorre == 1
    mkdir(output_path);
end
if ~exist(output_path_temporary, 'dir') && if_NoRMCorre == 1
    mkdir(output_path_temporary);
end

existingTemplateSession_path = [output_shortPath '\' existingTemplateSession];

template_fileName = 'template';
raw_fileName = 'Image_001_001.raw';
% raw_fileName = 'Image_001_001_jjbNew.raw';
corrected_fileName = 'normCorrected.h5';
normName = 'normcorre_gpuBatch_jjb';
parellelTestName = 'parellelTest';
binMeanName = 'binMean';
first10000Name = 'first10000';

% if ~exist(output_path, 'dir')
%     mkdir(output_path);
% end

saveName = [output_path '\' corrected_fileName(1:end-3) currentSession(13:end) corrected_fileName(end-2:end)];
saveName_temporary = [output_path_temporary '\' corrected_fileName(1:end-3) currentSession(13:end) corrected_fileName(end-2:end)];
Y = [currentSession_path '\' raw_fileName];
d1 = 512;
d2 = 512;

if_useExistingTemplate = 0;
if_nonRigidExistingTemplate = 1;
if_nonRigidRefine = 1;
if_innerSessionNonRigid = 1;
if_runBinMean = 1;
if_runFirst10000 = 0;
if_profileOn = 0;
if_delete_partialFile = 1;
if_runParellelTest = 0;
if_plot = 0;
if_metrics = 0;
if_spmd = 1;
if_dftreg_bacth = 0;
if_spmd_singleGPU = 0;
if_imageAverage = 0;
imageAverageBin = 2;
mini_batch_num = 6;%10-->12
% numFrames = 20000;% -1 means all frames
% numFrames = 4011;% -1 means all frames
numFrames = -1;% -1 means all frames
workerNum = 10;%11-->10
spmdCut_local = 40;%5-->6-->7-->9
spmdCut_global = 40;%5-->6-->7-->9
cpuThread_taskCoeff = 0.17;%0.65-->0.54-->0.45-->0.85-->0.65-->0.44
spmdCut_4090 = 40;
init_iter = 2;
finalOutputType = 'suite2pbin'; %hdf5

if if_NoRMCorre == 0
    return
end

%% Set parameters (rigid motion correction)
options_rigid = NoRMCorreSetParms(...
    'd1',d1,...
    'd2',d2,...
    'grid_size',[512,512],...
    'overlap_pre',[0 0],...
    'overlap_post',[0 0],...
    'us_fac',50,...%50
    'max_dev',[8 8],...
    'max_shift',[512 512],...  %[64 64]
    'phase_flag',false,... %true --> false
    'upd_template',false,... %true
    'init_batch',10000,...%10000
    'bin_width',21,...%22-->19-->21
    'iter',1,... %2
    'boundary','zero',...
    ...%'mem_batch_size',25,...
    'output_type','hdf5',...
    'h5_filename',saveName,...
    'use_windowing',false,... %true
    'correct_bidir',false);
if_FTpad = 1;
options_rigid.numFrames = numFrames;
options_rigid.if_metrics = if_metrics;
options_rigid.if_spmd = if_spmd;
options_rigid.workerNum = workerNum;
options_rigid.h5_filename_raw = saveName;
options_rigid.h5_filename_temporary = saveName_temporary;
options_rigid.h5_filename_temporary_raw = saveName_temporary;
options_rigid.if_dftreg_bacth = if_dftreg_bacth;
options_rigid.mini_batch_num = mini_batch_num;
options_rigid.spmdCut_local = spmdCut_local;
options_rigid.spmdCut_global = spmdCut_global;
options_rigid.cpuThread_taskCoeff = cpuThread_taskCoeff;
options_rigid.if_FTpad = if_FTpad;
options_rigid.init_iter = init_iter;
options_rigid.mem_batch_size = options_rigid.bin_width;
options_rigid.if_spmd_singleGPU = if_spmd_singleGPU;
options_rigid.spmdCut_4090 = spmdCut_4090;
options_rigid.if_imageAverage = if_imageAverage;
options_rigid.imageAverageBin = imageAverageBin;
if if_imageAverage == 1
    options_rigid.bin_width = imageAverageBin*options_rigid.bin_width;
end
options_rigid.if_useExistingTemplate = if_useExistingTemplate;
options_rigid.currentSession_path = currentSession_path;
options_rigid.if_nonRigidExistingTemplate = if_nonRigidExistingTemplate;
options_rigid.if_nonRigidRefine = if_nonRigidRefine;
% if if_useExistingTemplate == 0
%     options_rigid.bin_width = options_rigid.bin_width + 2;
% end
options_rigid.datestr_now30 = datestr_now30;
options_rigid.if_innerSessionNonRigid = if_innerSessionNonRigid;

spmd(0, workerNum)
    spmdindex = labindex;
    for tempi=1:100
        a = gpuArray(magic(1000)) * spmdindex;
        a = gather(a.^2); %#ok<*NASGU>
    end
end


%% Perform motion correction
normName_v = autoGetFunName(normName, [targetPATH '\functions']);
% normName_v(end-2:end) = 'v19';
fun_normcorre = str2func(normName_v);
fprintf('Now runing %s with %s.\n',currentSession,normName_v);
if if_profileOn == 1
    profile on
end


if if_useExistingTemplate == 1    
    temp_if_max0_min1 = 0;
    Result_existingTemplateSession_path = autoGetFileName_general('Result', existingTemplateSession_path,temp_if_max0_min1);
    
    file_list=dir([Result_existingTemplateSession_path,'\template*']);
    existingTemplate_fileName = [Result_existingTemplateSession_path,'\',file_list(end).name];
    existingTemplate = single(read_file(existingTemplate_fileName));    
    
    [~,parts_existingTemplate_fileName,~] = fileparts(existingTemplate_fileName);
    fprintf('Now use existing template from %s.\n',parts_existingTemplate_fileName);
else
    existingTemplate = [];
end

t1 = tic; 
% [M1,shifts1,template1,options_rigid,col_shift,cY,cM1] = fun_normcorre(Y,options_rigid,existingTemplate); 
[M1,shifts1,template1,options_rigid,col_shift,cY,cM1] = fun_normcorre(Y,options_rigid,existingTemplate); 
% [M1,shifts1,template1,options_rigid,col_shift,cY,cM1] = normcorre_gpuBatch_jjb_v111(Y,options_rigid,existingTemplate); 

fprintf('\nNorm correction ending,   time = %.1f seconds. ',toc(t1));
fprintf('Efficiency = %.1f fps.', options_rigid.numFrames/toc(t1));

if if_spmd == 1 && if_spmd_singleGPU == 0
    time_spmd = options_rigid.time_spmd;
    sum((options_rigid.numFrames/workerNum)./time_spmd);
    totalEfficiency = sum(options_rigid.bin_width*options_rigid.i_worker./time_spmd);
    fprintf(' Total efficiency = %.1f fps.\n', totalEfficiency);
else
    fprintf('\n');
end
spmd, end % To release memory from spmd

if if_profileOn == 1
    profile viewer
end


%% Compute correlation threhold
if if_metrics == 1
    % outlierIndex_cM1 = isoutlier(cM1,'median');
    % outlierIndex_cM1 = isoutlier(cM1,'mean');
    outlierIndex_cM1 = isoutlier(cM1,'grubbs');
    sum(outlierIndex_cM1);
    if sum(outlierIndex_cM1)>0
        outlierThreshold_cM1 = max(cM1(outlierIndex_cM1));
    else
        outlierThreshold_cM1 = min(cM1);
    end
end
if if_plot == 1
    %% Plot template
    fig1 = figure(1);
    set(gcf,'Position',[5 35 350 350]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set (gca,'position',[0.01,0.01,0.98,0.98] )
    
    imagesc(template1.^(1/3));
    %imagesc(single(existingTemplate).^(1/3));
    %imagesc(template_in.^(1/3));
    my_gray = gray;
    colormap(my_gray);
    axis off
    axis equal
    
    %% Plot motion correction correlation
    fig2 = figure(2);
    set(gcf,'Position',[360 35 700 350]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
    
    plot(cY, '-', 'LineWidth', 1, 'Color', [0.40 0.40 0.40]);
    hold on
    plot(cM1, '-', 'LineWidth', 1, 'Color', [0.25 0.85 0.85]);
    hold on
    plot([0 length(cM1)], [outlierThreshold_cM1 outlierThreshold_cM1],...
        '-', 'LineWidth', 1, 'Color', [0.75 0.25 0.25]);
    
    hl = legend('Raw','Motion correction','Threshold',...
        'Location','southwest','fontsize',10);
    ylim([0 1]);
    set(gca, 'FontSize', 14)
    set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('Frame number, T=%d',length(cY)), 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Correlation coefficient', 'FontSize', 14, 'FontWeight', 'bold');
    title('Motion Correction Correlation', 'FontSize', 14);
    
end

[~,parts_fname,parts_ext] = fileparts(M1);
% M1_raw = M1;
% M1 = fullfile(output_path,[parts_fname,parts_ext]);



opts_tiff.overwrite = true;
opts_tiff.message = false;
saveastiff(uint16(template1),[output_path,'\',template_fileName,parts_fname(14:end),'.tif'],opts_tiff);


% template_2 = zeros(512,512,2);
% template_2(:,:,1) = template1;
% template_2(:,:,2) = template1;
% template_2_T = pagetranspose(template_2);
% template_fileFullName = [output_path,'\template_2.raw'];
% 
% fid = fopen(template_fileFullName,'w');
% fwrite(fid,template_2_T,'uint16',0,'l');
% fclose(fid);


% if_spmd=0;
if if_spmd == 1
    [parts_pathstr,parts_fname,h5_ext] = fileparts(M1);
    M1_Bin = fullfile(parts_pathstr,[parts_fname(1:13),'Bin',parts_fname(14:end),h5_ext]);
    file_list=dir([output_path_temporary,'\','partial_',parts_fname,'*']);
    file_list_bin=dir([output_path_temporary,'\','partial_',[parts_fname(1:13) 'Bin' parts_fname(14:end)],'*']);
     
    length_list = length(file_list);
    for tempi=1:length_list
        file_list(tempi).name_old = file_list(tempi).name;
        file_list(tempi).name = [output_path_temporary '\' file_list(tempi).name_old];
    end
        
    length_list_bin = length(file_list_bin);
    for tempi=1:length_list_bin
        file_list_bin(tempi).name_old = file_list_bin(tempi).name;
        file_list_bin(tempi).name = [output_path_temporary '\' file_list_bin(tempi).name_old];
    end
    concatenateFunName = 'concatenate_files_forParellel_jjb';
    concatenateFun_v = autoGetFunName(concatenateFunName, [targetPATH '\functions\PartA_PreProcess']);
    fun_concatenate = str2func(concatenateFun_v);

    %fprintf('Now is concatenating! (Bin, fast)...');
    fprintf('Concatenating Bin...');
    t0_cat1 = tic;
    [filename_bin_output,length_list_bin_output] = fun_concatenate(file_list_bin,M1_Bin,'hdf5');    
    %fprintf('Concatenate bin ending, time = %.1f seconds.\n',toc(t0_cat1));
    fprintf('Done, time = %.1f seconds.\n',toc(t0_cat1));
    if if_runBinMean == 1
        binMeanName_v = autoGetFunName(binMeanName, [targetPATH '\functions']);
        fprintf('Now runing %s. ', binMeanName_v);
        fun_binMean = str2func(binMeanName_v);
        temp_imageAverageBin = 2;
        temp_if_raw0_h51 = 1;
        binMeanFileName = fun_binMean(output_path_temporary,output_path,temp_imageAverageBin,temp_if_raw0_h51);        
    end    
    
    %fprintf('Now is concatenating! (Very time-consuming)...');
    fprintf('Concatenating all...');
    t0_cat2 = tic;
    %[filename_output,length_list_output,error_flag] = fun_concatenate(file_list,M1,'hdf5');
    [filename_output,length_list_output,error_flag] = fun_concatenate(file_list,M1,finalOutputType,if_delete_partialFile);
    %fprintf('Concatenate     ending, time = %.1f seconds.\n',toc(t0_cat2));
    fprintf('Done, time = %.1f seconds.\n',toc(t0_cat2));
else
    
    [parts_pathstr,parts_fname,h5_ext] = fileparts(M1); %#ok<*ASGLU>
    M1_Bin = fullfile(parts_pathstr,[parts_fname(1:13),'Bin',parts_fname(14:end),h5_ext]);
    
    file_list=dir([output_path_temporary,'\',parts_fname,'*']);    
    file_list_bin=dir([output_path_temporary,'\',[parts_fname(1:13) 'Bin' parts_fname(14:end)],'*']);
         
    length_list = length(file_list);
    for tempi=1:length_list
        file_list(tempi).name_old = file_list(tempi).name;
        file_list(tempi).name = [output_path_temporary '\' file_list(tempi).name_old];
    end
    
    length_list_bin = length(file_list_bin);
    for tempi=1:length_list_bin
        file_list_bin(tempi).name_old = file_list_bin(tempi).name;
        file_list_bin(tempi).name = [output_path_temporary '\' file_list_bin(tempi).name_old];
    end    
    
    concatenateFunName = 'concatenate_files_forParellel_jjb';
    concatenateFun_v = autoGetFunName(concatenateFunName, [targetPATH '\functions\PartA_PreProcess']);
    fun_concatenate = str2func(concatenateFun_v);
    
    if length(output_path_temporary) ~= length(output_path)
        [filename_bin_output,length_list_bin_output] = fun_concatenate(file_list_bin,M1_Bin,'hdf5');
        % Delete redundant h5 bin file in temporary dir
        recycle('off');
        delete(file_list_bin.name);
        recycle('on');        
    end
    
    if if_runBinMean == 1
        binMeanName_v = autoGetFunName(binMeanName, [targetPATH '\functions']);
        fprintf('Now runing %s for suite2p. ', binMeanName_v);
        fun_binMean = str2func(binMeanName_v);
        temp_imageAverageBin = 2;
        temp_if_raw0_h51 = 1;
        binMeanFileName = fun_binMean(output_path_temporary,output_path,temp_imageAverageBin,temp_if_raw0_h51);
    end
    
    [filename_output,length_list_output,error_flag] = fun_concatenate(file_list,M1,finalOutputType);
    if strcmp(finalOutputType,'suite2pbin') == 1
        % Delete redundant h5 file
        recycle('off');
        delete(file_list.name);
        recycle('on');
    end

end
% profile on
if if_runFirst10000 == 1
    first10000Name_v = autoGetFunName(first10000Name, [targetPATH '\functions']);
    fprintf('Now runing %s. ', first10000Name_v);
    fun_first10000 = str2func(first10000Name_v);
    fun_first10000(output_path_temporary,output_path,finalOutputType);
end
% profile viewer
if if_delete_partialFile == 1
    %recycle('on');
    recycle('off');
    delete([output_path_temporary '\partial*']);
    %delete([currentSession_path '\partial*']);
    recycle('on');
end
% if length(output_path_temporary) ~= length(output_path)
%     % Delete redundant temporary dir
%     recycle('off');
%     rmdir(output_path_temporary);
%     recycle('on');    
% end

if strcmp(finalOutputType,'suite2pbin') == 1
    % To rename corrected file as 'data.bin'
    movefile(filename_output,[parts_pathstr,'\data.bin']); 
    
    % To rename corrected binned file as 'binnedData.h5'
    %movefile(filename_bin_output,[parts_pathstr,'\binnedData.h5']);
    movefile(binMeanFileName,[parts_pathstr,'\binnedData.h5']);
    
    if if_runBinMean == 1
        fprintf('Now runing %s again for a raw file. ', binMeanName_v);
        temp_imageAverageBin = 4;
        temp_if_raw0_h51 = 0;
        temp_if_deleteBin = 1;
        fun_binMean(output_path_temporary,output_path,temp_imageAverageBin,temp_if_raw0_h51,temp_if_deleteBin);
    end     
end

if if_runParellelTest == 1
    parellelTestName_v = autoGetFunName(parellelTestName, [targetPATH '\functions']);
    fprintf('Now runing %s. ', parellelTestName_v);
    fun_parellelTest = str2func(parellelTestName_v);
    fun_parellelTest(output_path,currentSession,if_imageAverage,finalOutputType);
end

% % To create a fake h5 file to trick on suite2p
% temph5_fileFullName = [output_path,'/temp.h5'];
% if exist(temph5_fileFullName,'file') ~= 0
%     delete(temph5_fileFullName);
% end
% h5create(temph5_fileFullName,['/mov'],[512,512,Inf],'Chunksize',[512,512,1],'Datatype','uint16');
% h5write(temph5_fileFullName,['/mov'],...
%     uint16(rand(512,512,51)),...
%     [1,1,1],...
%     [512,512,51]);

% To create link files of cmd batch files
fileName_runScript_ForSuite2p = 'runScript_ForSuite2p.cmd.lnk';
fileName_guiScript_ForSuite2p = 'guiScript_ForSuite2p.cmd.lnk';
% linkPath = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\others';
linkPath = [output_shortPath,'\others'];
copyfile([linkPath,'\',fileName_runScript_ForSuite2p],[output_path,'\',fileName_runScript_ForSuite2p]);
copyfile([linkPath,'\',fileName_guiScript_ForSuite2p],[output_path,'\',fileName_guiScript_ForSuite2p]);

fprintf('---------------------------------------------------------------------------------------\n');

fprintf('myNoRMCorre ending,       time = %.1f seconds. ',toc(t0));

% delete(gcp);
%% End