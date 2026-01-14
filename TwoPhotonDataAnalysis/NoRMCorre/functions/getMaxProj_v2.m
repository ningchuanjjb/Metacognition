function getMaxProj_v2(output_path)
%% Initialization

a = 1; %#ok<*NASGU>

if(~exist('output_path','var'))
    output_path = 'D:\twoPhotonData_motionCorrected\113Recording_20230528A_Ding_Site02B_sameFOV0527\Result20230722T043326';
end

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)


binMean_fileName = 'normCorrectedBinMean_';


binMean_fileFullName = autoGetFileName(binMean_fileName, output_path);
[parts_pathstr,parts_fname,~] = fileparts(binMean_fileFullName); %#ok<*ASGLU>
fprintf('Load %s.\n', parts_fname);
t0 = tic;

fid = fopen(binMean_fileFullName);
imsize = 512*512*2;
current_seek = ftell(fid);
fseek(fid, 0, 1);
file_length = ftell(fid);
fseek(fid, current_seek, -1);
frameNum = floor(file_length/imsize);
T = frameNum;
fprintf('Load %d frames (from %d).', T, frameNum);
Y_binMean = fread(fid,512*512*T,'uint16=>uint16',0,'l')';
Y_binMean = reshape(Y_binMean,[512 512 T]);
fclose(fid);

Y_maxProj = max(Y_binMean,[],3);
Y_maxProj = Y_maxProj';

maxProj_tif_fileFullName = fullfile(output_path,'maxProjection.tif');

opts_tiff.overwrite = true;
opts_tiff.message = false;
saveastiff(uint16(Y_maxProj),maxProj_tif_fileFullName,opts_tiff);
fprintf('Save max projection image.\n');


