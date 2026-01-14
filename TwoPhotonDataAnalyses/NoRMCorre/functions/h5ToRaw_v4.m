function h5ToRaw_v4(intput_path,output_path)
%% Initialization


if(~exist('intput_path','var'))
    intput_path = 'R:\twoPhotonRawData\ToBeProcessed\329Recording_20231227A_YRK';
end

if(~exist('output_path','var'))
    %output_path = 'D:\twoPhotonData_motionCorrected\329Recording_20231227A_YRK';
    output_path = 'R:\twoPhotonRawData\ToBeProcessed\329Recording_20231227A_YRK';    
end

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)



inputH5_fileName = '*.h5';


inputH5_fileFullName = autoGetFileName(inputH5_fileName, intput_path);
[parts_pathstr,parts_fname,~] = fileparts(inputH5_fileFullName); %#ok<*ASGLU>
t0 = tic;


fprintf('Load %s.\n', parts_fname);
fileinfo = hdf5info(inputH5_fileFullName); %#ok<*HDFI>
sizY = fileinfo.GroupHierarchy.Datasets.Dims;
frameNum = sizY(end);
fprintf('Load %d frames.', frameNum);


outputRaw_fileFullName = fullfile(output_path,[parts_fname,'.raw']); %#ok<*NASGU>


workerNum = 8;%9-->6
bin_width = 20;%25

data_name = fileinfo.GroupHierarchy.Datasets.Name;
T = frameNum;
% T = min(20000,frameNum);


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


fid = fopen(outputRaw_fileFullName,'w');
imsize = 512*512*2;



for spmdindex=1:workerNum

    tempi_worker = i_worker_start(spmdindex):i_worker_end(spmdindex);
    tempi_length = length(tempi_worker);
    for tempj=1:tempi_length
        n = tempi_worker(tempj);
        bin_count = n;
        t = (bin_count-1)*bin_width + 1;

        Y = h5read(inputH5_fileFullName,data_name,[1,1,t],[512,512,min(t+bin_width-1,T)-t+1]);
                        
        size_3_Y = size(Y,ndims(Y));
        
        
        feek_indicator = (t-1)*imsize;
        fseek(fid,feek_indicator,'bof');
        fwrite(fid,uint16(Y),'uint16',0,'l');

    end        
end
fclose(fid);

fprintf('Save %d frames in raw, time = %.2f seconds\n',T,toc(t0));


