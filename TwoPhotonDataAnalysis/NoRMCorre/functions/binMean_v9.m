function saveName = binMean_v9(intput_path,output_path,imageAverageBin,if_raw0_h51,if_deleteBin)
%% Initialization

a = 1;

gcp;% To open parallel pool
if false
    delete(gcp); %#ok<*UNRCH>
    gcp;
end

if(~exist('output_path','var'))
    %output_path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\20230128T113714';
    %output_path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230508A_Ding_Site02_sameFOV0504\Result20230509T132354';    
    output_path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230508A_Ding_Site02_sameFOV0504\Result20230509T125103';        
end
if(~exist('if_deleteBin','var'))
    if_deleteBin = 0;%1
end
if(~exist('if_raw0_h51','var'))
    if_raw0_h51 = 1;%0
end
if(~exist('imageAverageBin','var'))
    imageAverageBin = 2;%4-->2
end

if_profileOn = 0;
if_imageAverage = 1;

if if_profileOn == 1
    profile on
end

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)

corrected_fileName = 'normCorrectedBin_';


corrected_fileFullName = autoGetFileName(corrected_fileName, intput_path);
[h5_pathstr,h5_fname,h5_ext] = fileparts(corrected_fileFullName); %#ok<*ASGLU>

correctedBinMean_fileFullName = fullfile(output_path,[h5_fname(1:16),'Mean',h5_fname(17:end),h5_ext]);
correctedBinMean_raw_fileFullName = fullfile(output_path,[h5_fname(1:16),'Mean',h5_fname(17:end),'.raw']);

fprintf('Load %s.\n', h5_fname);
fileinfo = hdf5info(corrected_fileFullName); %#ok<*HDFI>
sizY = fileinfo.GroupHierarchy.Datasets.Dims;
frameNum = sizY(end);
fprintf('Load %d frames.', frameNum);

t0 = tic;
        

workerNum = 8;%9-->6
bin_width = 24;%25

data_name = fileinfo.GroupHierarchy.Datasets.Name;
T = frameNum;


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

if if_raw0_h51 == 1
    if exist(correctedBinMean_fileFullName,'file') ~= 0
        delete(correctedBinMean_fileFullName);
    end
    h5create(correctedBinMean_fileFullName,['/mov'],[512,512,Inf],'Chunksize',[512,512,1],'Datatype','uint16'); %#ok<*NBRAK>

elseif if_raw0_h51 == 0
    fid = fopen(correctedBinMean_raw_fileFullName,'w');
    imsize = 512*512*2;
end

for spmdindex=1:workerNum


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
        
        % Huge bug here! It cause binMean.m create wrong binMean image!!!
        %Y_corrected = gpuArray(h5read(corrected_fileFullName,data_name,[1,1,t_valid],[512,512,min(t+bin_width-1,T)-t+1]));
                
        Y_corrected = gpuArray(h5read(corrected_fileFullName,data_name,[1,1,t],[512,512,min(t+bin_width-1,T)-t+1]));        
        
        size_3_Y_corrected = size(Y_corrected,ndims(Y_corrected));
        
        lY_raw = Y_corrected;
        
        %if lY_raw == 1
        if ndims(Y_corrected) == 2 %#ok<*ISMAT>
            continue
        else
            if if_imageAverage == 1
                if ceil(size_3_Y_corrected/imageAverageBin)*imageAverageBin == size_3_Y_corrected
                    Y_corrected_mean = squeeze(mean(reshape(Y_corrected,[512,512,imageAverageBin,size_3_Y_corrected/imageAverageBin]), 3));
                else
                    Y_corrected_mean = squeeze(mean(reshape(...
                        Y_corrected(:,:,1:floor(size_3_Y_corrected/imageAverageBin)*imageAverageBin),[512,512,imageAverageBin,floor(size_3_Y_corrected/imageAverageBin)]), 3));
                    tempi = (floor(size_3_Y_corrected/imageAverageBin)*imageAverageBin+1):size_3_Y_corrected;
                    Y_corrected_mean(:,:,(floor(size_3_Y_corrected/imageAverageBin)+1)) = mean(Y_corrected(:,:,tempi), 3);
                end
                
                Y_corrected = Y_corrected_mean;
                Y_corrected_mean = []; %#ok<*NASGU>
                size_3_Y_corrected = size(Y_corrected,ndims(Y_corrected));
            end
        end
        if ndims(Y_corrected) == 2 %#ok<*ISMAT>
            continue
        end
        
        lY = size_3_Y_corrected;
        if lY == 1
            continue            
        end   
        
        a = 1;
        
        if if_imageAverage == 1            
            temp_bin_width = ceil(bin_width/imageAverageBin);
        else
            temp_bin_width = bin_width;
        end
        
        if if_raw0_h51 == 1
            h5write(correctedBinMean_fileFullName,['/mov'],...
                gather(uint16(Y_corrected)),...
                [ones(1,2),t_valid],...
                [512,512,lY]);
            
        elseif if_raw0_h51 == 0
            feek_indicator = (t_valid-1)*imsize;
            fseek(fid,feek_indicator,'bof');
            fwrite(fid,gather(uint16(Y_corrected)),'uint16',0,'l');
        end
    end        
end
if if_raw0_h51 == 0
    fclose(fid);
end
fprintf('Save %d frames, time = %.2f seconds\n',t_valid+lY,toc(t0));

if if_deleteBin == 1
    if exist(corrected_fileFullName,'file') ~= 0
        recycle('off');
        delete(corrected_fileFullName);
        recycle('on');
    end
end

if if_profileOn == 1
    profile viewer
end

saveName = correctedBinMean_fileFullName;
