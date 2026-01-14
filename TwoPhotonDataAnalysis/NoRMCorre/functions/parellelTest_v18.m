function [cY, cM] = parellelTest_v18(output_path,currentSession,if_imageAverage,output_type,if_plot,if_max0_min1)
%% Initialization
% clear
% close all

% gpuDevice(1);

gcp;% To open parallel pool
%delete(gcp)
if false
    delete(gcp); %#ok<*UNRCH>
    gcp;
end

output_type_flag_h50_bin1 = [];
if ~exist('output_type','var')
    output_type = 'suite2pbin';
end

if strcmp(output_type,'hdf5') == 1
    output_type_flag_h50_bin1 = 0;
    corrected_fileName = 'normCorrected_';
elseif strcmp(output_type,'suite2pbin') == 1
    output_type_flag_h50_bin1 = 1;
    corrected_fileName = 'data.bin';
end

if ~exist('currentSession','var')
    %currentSession = '113Recording_20221211A_Bubble_fixZtest';
    %currentSession = '113Recording_20230106A_Ding_Site09_sameFOV0105A_justTest';
    %currentSession = '113Recording_20230109A_Ding_Site09';
    currentSession = '113Recording_20230123A_Ding_Site09B_sameFOV0122';
end
if(~exist('output_path','var'))
    %output_path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\20230128T113714';
    output_path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\Result20230205T170425';    
end
if ~exist('if_imageAverage','var')
    if_imageAverage = 0;
end
if ~exist('if_plot','var')
    if_plot = 1;
end
if ~exist('if_max0_min1','var')
    if_max0_min1 = 0;
end
if_profileOn = 0;
% if_imageAverage = 1;
imageAverageBin = 2;

if if_profileOn == 1
    profile on
end

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)

% currentSession = '113Recording_Bubble_20221211A_fixZtest';

% rawData_path = 'D:\twoPhotonRawData\ToBeProcessed';
rawData_path = 'R:\twoPhotonRawData\ToBeProcessed';
currentSession_path = [rawData_path '\' currentSession];
% output_shortPath = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected';
% output_path = [output_shortPath '\' currentSession];


raw_fileName = 'Image_001_001.raw';

raw_fileFullName = [currentSession_path '\' raw_fileName];


t0 = tic;

if output_type_flag_h50_bin1 == 0
    corrected_fileFullName = autoGetFileName(corrected_fileName, output_path, if_max0_min1);
    [~,parts_fname,parts_ext] = fileparts(corrected_fileFullName);
    fprintf('Load %s.\n', [parts_fname,'.',parts_ext]);

    fileinfo = hdf5info(corrected_fileFullName); %#ok<*HDFI>
    sizY = fileinfo.GroupHierarchy.Datasets.Dims;
    frameNum = sizY(end);
    data_name = fileinfo.GroupHierarchy.Datasets.Name;
    
elseif output_type_flag_h50_bin1 == 1
    corrected_fileFullName = [output_path,'\',corrected_fileName];
    fprintf('Load %s.\n', corrected_fileName);

    fid_corrected = fopen(corrected_fileFullName);
    imsize = 512*512*2;
    fseek(fid_corrected, 0, 1);
    file_length = ftell(fid_corrected);
    frameNum = file_length/imsize;    
    fclose(fid_corrected);    
end

fprintf('Load %d frames.', frameNum);

        
% Y_corrected = read_file(corrected_fileFullName,1,frameNum);
% Y_raw = read_raw_file(raw_fileFullName,1,frameNum,[512 512],2);
% cM = motion_metrics(Y_corrected,10);
% cY = motion_metrics(Y_raw,10);

workerNum = 8;%9-->6
bin_width = 24;%25

T = frameNum;


% cY = zeros(1, T);
% cM = zeros(1, T);
% cY = zeros(1,T,'single','gpuArray');
% cM = zeros(1,T,'single','gpuArray');

if if_imageAverage == 1
    T = (T-1)*imageAverageBin + 1;
    
%     fid_raw = fopen(raw_fileFullName);
%     imsize = 512*512*2;
%     fseek(fid_raw, 0, 1);
%     file_length = ftell(fid_raw);
%     T_raw = round(file_length/imsize);        
    
    
    cY = zeros(1, ceil(T/imageAverageBin));
    cM = zeros(1, ceil(T/imageAverageBin));
else
    cY = zeros(1, T);
    cM = zeros(1, T);
end

bin_totalNum = length(1:bin_width:T);
N = bin_totalNum;

i_worker = floor(N/workerNum) * ones(1, workerNum);
i_worker(1:N - sum(i_worker)) = i_worker(1:N - sum(i_worker)) + 1;
i_worker_start = zeros(1, workerNum);
i_worker_end = zeros(1, workerNum);
for tempj=1:workerNum
    if tempj == 1
        i_worker_start(1) = 1;
    else
        i_worker_start(tempj) = sum(i_worker(1:tempj-1))+1;
    end
    i_worker_end(tempj) = i_worker_start(tempj)+i_worker(tempj)-1;
end

imsize = 512*512*2;

% for spmdindex=1:workerNum
spmd(0, workerNum)
    spmdindex = labindex;
        
    fid_raw = fopen(raw_fileFullName);                    
    if output_type_flag_h50_bin1 == 1
        fid_corrected = fopen(corrected_fileFullName);
    end
        
    tempi_worker = i_worker_start(spmdindex):i_worker_end(spmdindex);
    tempi_length = length(tempi_worker);
    for tempj=1:tempi_length
        n = tempi_worker(tempj);
        bin_count = n;
        t = (bin_count-1)*bin_width + 1;
        if if_imageAverage == 1
            t_valid = ceil(t/imageAverageBin);
        else
            t_valid = t;
        end
        
        window = min(t+bin_width-1,T)-t+1;
        feek_indicator = (t-1)*imsize;
        fseek(fid_raw,feek_indicator,'bof');
        Y_raw = gpuArray(fread(fid_raw,512*512*window,'uint16=>single',0,'l')');
        Y_raw = Y_raw/2;% Preprocess for suite2p bin (uint16-->int16, although not now)
        Y_raw = reshape(Y_raw,[512 512 window]);
        Y_raw = pagetranspose(Y_raw);
        size_3_Yraw = size(Y_raw,ndims(Y_raw));
        
        lY_raw = size_3_Yraw;
        
        %if t_valid == 1501
        %    ndims(Y_raw)
        %    a = 1;
        %end
        
        %if lY_raw == 1
        if ndims(Y_raw) == 2 %#ok<*ISMAT>
            cY(t_valid) = cY(t_valid-1);
            cM(t_valid) = cM(t_valid-1);
            continue
        else
            if if_imageAverage == 1
                if ceil(size_3_Yraw/imageAverageBin)*imageAverageBin == size_3_Yraw
                    Yraw_mean = squeeze(mean(reshape(Y_raw,[512,512,imageAverageBin,size_3_Yraw/imageAverageBin]), 3));
                else
                    Yraw_mean = squeeze(mean(reshape(...
                        Y_raw(:,:,1:floor(size_3_Yraw/imageAverageBin)*imageAverageBin),[512,512,imageAverageBin,floor(size_3_Yraw/imageAverageBin)]), 3));
                    tempi = (floor(size_3_Yraw/imageAverageBin)*imageAverageBin+1):size_3_Yraw;
                    Yraw_mean(:,:,(floor(size_3_Yraw/imageAverageBin)+1)) = mean(Y_raw(:,:,tempi), 3);
                end
                
                Y_raw = Yraw_mean;
                Yraw_mean = []; %#ok<*NASGU>
                size_3_Yraw = size(Y_raw,ndims(Y_raw));
            end
        end
        lY = size(Y_raw,3);
        if lY == 1
            continue            
        end
        
        if output_type_flag_h50_bin1 == 0
            Y_corrected = gpuArray(h5read(corrected_fileFullName,data_name,[1,1,t_valid],[512,512,lY]));
        elseif output_type_flag_h50_bin1 == 1
            fseek(fid_corrected,feek_indicator,'bof');
            Y_corrected = gpuArray(fread(fid_corrected,512*512*window,'int16=>single',0,'l')');
            Y_corrected = reshape(Y_corrected,[512 512 window]);
        end
      
        cY(t_valid:t_valid+lY-1) = motion_metrics_gpuBatch_jjb_v1(Y_raw,10);
        cM(t_valid:t_valid+lY-1) = motion_metrics_gpuBatch_jjb_v1(Y_corrected,10);                        
    end  
    fclose(fid_raw);
    if output_type_flag_h50_bin1 == 1
        fclose(fid_corrected);
    end    
end

if if_imageAverage == 1
    temp_T = ceil(T/imageAverageBin);
    temp_bin_width = ceil(bin_width/imageAverageBin);
else
    temp_T = T;
    temp_bin_width = bin_width;
end
    
    
cY_all = [];
cM_all = [];
for tempi=1:workerNum
    temp_cY = cY{tempi};
    temp_cM = cM{tempi};
    if tempi == 1
        temp_range = 1:i_worker_end(tempi)*temp_bin_width;
    else
        temp_range = (i_worker_end(tempi-1))*temp_bin_width+1:i_worker_end(tempi)*temp_bin_width;
    end
    temp_range(temp_range>temp_T) = [];
    cY_all = [cY_all temp_cY(temp_range)]; %#ok<*AGROW>
    cM_all = [cM_all temp_cM(temp_range)];
end
cY = cY_all;
cM = cM_all;


outlierIndex_cM1 = isoutlier(cM,'grubbs');
sum(outlierIndex_cM1);
if sum(outlierIndex_cM1)>0
    outlierThreshold_cM1 = max(cM(outlierIndex_cM1));
else
    outlierThreshold_cM1 = min(cM);
end

fprintf('time = %.2f seconds.\n', toc(t0));


%% Plot motion correction correlation
if if_plot == 1
    fig11 = figure(11);
    set(gcf,'Position',[360 35 700 350]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
    
    plot(cY, '-', 'LineWidth', 1, 'Color', [0.40 0.40 0.40]);
    hold on
    plot(cM, '-', 'LineWidth', 1, 'Color', [0.25 0.85 0.85]);
    hold on
    plot([0 length(cM)], [outlierThreshold_cM1 outlierThreshold_cM1],...
        '-', 'LineWidth', 1, 'Color', [0.75 0.25 0.25]);
    
    hl = legend('Raw','Motion correction','Threshold',...
        'Location','southwest','fontsize',10);
    ylim([0 1]);
    set(gca, 'FontSize', 14);
    set(gca,'YTick',[0:0.2:1]);%#ok<*NBRAK> %设置要显示坐标刻度的范围
    set(gca,'box','off');% 取消右、上边框
    xlabel(sprintf('Frame number, T=%d',length(cY)), 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Correlation coefficient', 'FontSize', 14, 'FontWeight', 'bold');
    title('Motion Correction Correlation', 'FontSize', 14);

end
% spmd(0, 16)
%     gpuDevice([]);    
% end

if if_profileOn == 1
    profile viewer
end
