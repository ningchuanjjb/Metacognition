clear all
close all


% fileName1 = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\test1\plane0\data.bin';
fileName1 = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\test1\normCorrected_20230123A_Ding_Site09B_sameFOV0122_20230126T114039.h5';
fileName2 = 'mydata.bin';

% fid = fopen(fileName1);
% imsize = 512*512*2;
% current_seek = ftell(fid);
% fseek(fid, 0, 1);
% file_length = ftell(fid);
% fseek(fid, current_seek, -1);
% frame_length = file_length/imsize;
% Y_raw1 = fread(fid,512*512*frame_length,'int16=>uint16',0,'l')';
% Y_raw1 = single(2*Y_raw1);
% Y_raw1 = reshape(Y_raw1,[512 512 frame_length]);
% % Y_raw1 = pagetranspose(Y_raw1);
% fclose(fid);

info = h5info(fileName1);
dims = info.Datasets.Dataspace.Size;
sframe = 1;
num2read = dims(end)-sframe+1;
Y_h5 = h5read(fileName1,'/mov',[1,1,1],[512,512,num2read]);


fid2 = fopen(fileName2,'w');
fwrite(fid2,int16(uint16(Y_h5)/2),'int16',0,'l');
fclose(fid2);




