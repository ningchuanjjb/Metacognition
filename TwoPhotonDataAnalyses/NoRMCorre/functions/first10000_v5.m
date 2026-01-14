function first10000_v5(intput_path,output_path,output_type)
%% Initialization

a = 1;

gcp;% To open parallel pool
if false
    delete(gcp); %#ok<*UNRCH>
    gcp;
end

output_type_flag_h50_bin1 = [];
if strcmp(output_type,'hdf5') == 1
    output_type_flag_h50_bin1 = 0;
elseif strcmp(output_type,'suite2pbin') == 1
    output_type_flag_h50_bin1 = 1;
end

if(~exist('output_path','var'))
    output_path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\20230128T113714';
end

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)


corrected_fileName = 'normCorrected_';


corrected_fileFullName = autoGetFileName(corrected_fileName, intput_path);
[parts_pathstr,parts_fname,~] = fileparts(corrected_fileFullName); %#ok<*ASGLU>
fprintf('Load %s.\n', parts_fname);
t0 = tic;

if output_type_flag_h50_bin1 == 0
    fileinfo = hdf5info(corrected_fileFullName); %#ok<*HDFI>
    sizY = fileinfo.GroupHierarchy.Datasets.Dims;
    frameNum = sizY(end);
    data_name = fileinfo.GroupHierarchy.Datasets.Name;
    T = min(10000,frameNum);
    fprintf('Load %d frames (from %d).', T, frameNum);
    Y_corrected_first10000 = h5read(corrected_fileFullName,data_name,[1,1,1],[512,512,T]);
    
elseif output_type_flag_h50_bin1 == 1
    fid = fopen(corrected_fileFullName);
    imsize = 512*512*2;
    current_seek = ftell(fid);
    fseek(fid, 0, 1);
    file_length = ftell(fid);
    fseek(fid, current_seek, -1);
    frameNum = floor(file_length/imsize);
    T = min(10000,frameNum);
    fprintf('Load %d frames (from %d).', T, frameNum);
    Y_corrected_first10000 = fread(fid,512*512*T,'int16=>uint16',0,'l')';
    Y_corrected_first10000 = reshape(Y_corrected_first10000,[512 512 T]);
    fclose(fid);    
end



correctedFirst10000_raw_fileFullName = fullfile(output_path,...
    [parts_fname(1:length(corrected_fileName)-1),'First',num2str(T),parts_fname(length(corrected_fileName):end),'.raw']);

fid = fopen(correctedFirst10000_raw_fileFullName,'w');
fwrite(fid,Y_corrected_first10000,'uint16',0,'l');
fclose(fid);

fprintf('Save first %d frames in raw, time = %.2f seconds\n',T,toc(t0));


