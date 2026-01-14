% function saveName = binMean_v8(output_path,imageAverageBin,if_raw0_h51,if_deleteBin)
%% Initialization
% clear
% close all


% output_path = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230508A_Ding_Site02_sameFOV0504\Result20230509T125103';
output_path = 'R:\twoPhotonRawData\ToBeProcessed\113Recording_20230510A_Ding_Site05';


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)


raw_fileName = 'Image_001_001_jjb';
newRaw_fileName = 'Image_001_001_jjbNew';

raw_fileFullName = autoGetFileName(raw_fileName, output_path);
newRaw_fileFullName = [output_path,'\',newRaw_fileName,'.raw'];

fid = fopen(raw_fileFullName);
imsize = 512*512*2;
current_seek = ftell(fid);
fseek(fid, 0, 1);
file_length = ftell(fid);
fseek(fid, current_seek, -1);
frame_length = file_length/imsize;


[~,fname,~] = fileparts(raw_fileFullName);
fprintf('Load %s.\n', fname);
fprintf('Load %d frames.', frame_length);


fid_new = fopen(newRaw_fileFullName,'w');


t0 = tic;


T = frame_length;
% T = 1000;

bin_width = 24;
bin_totalNum = length(1:bin_width:T);

savedFramesNum = 0;

for bin_count = 1:bin_totalNum

    t = (bin_count-1)*bin_width + 1;
    
    window = min(t+bin_width-1,T)-t+1;
    feek_indicator = (t-1)*imsize;
    fseek(fid,feek_indicator,'bof');
    Y_raw = fread(fid,512*512*window,'uint16=>uint16',0,'l')';
    Y_raw = reshape(Y_raw,[512 512 window]);
    
    size_3_Y_raw = size(Y_raw,ndims(Y_raw));
    lY = size_3_Y_raw;
    
    a = 1;
    temp_range = 1:2:lY;
    Y_newRaw = Y_raw(:,:,temp_range);
    lY_newRaw = size(Y_newRaw,3);
    savedFramesNum = savedFramesNum + lY_newRaw;

    fwrite(fid_new,Y_newRaw,'uint16',0,'l');    
end
fclose(fid);
fclose(fid_new);
fprintf('Save %d frames, time = %.2f seconds\n',savedFramesNum,toc(t0));

