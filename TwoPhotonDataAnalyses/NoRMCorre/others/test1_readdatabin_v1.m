clear all
close all

% fid = fopen('data.bin','r');
% [data, file_length] = fread(fid,Inf,'*uint16');

fileName1 = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\test1\plane0\data.bin';
fileName2 = 'C:\ASDROOT\STUDY\twoPhotonData_motionCorrected\113Recording_20230123A_Ding_Site09B_sameFOV0122\test1\normCorrected_20230123A_Ding_Site09B_sameFOV0122_20230126T114039.h5';
fileName3 = 'mydata.bin';

fid = fopen(fileName1);
imsize = 512*512*2;
current_seek = ftell(fid);
fseek(fid, 0, 1);
file_length = ftell(fid);
fseek(fid, current_seek, -1);
frame_length = file_length/imsize;
Y_raw1 = fread(fid,512*512*frame_length,'int16=>uint16',0,'l')';
Y_raw1 = single(2*Y_raw1);
Y_raw1 = reshape(Y_raw1,[512 512 frame_length]);
% Y_raw1 = pagetranspose(Y_raw1);
fclose(fid);

info = h5info(fileName2);
dims = info.Datasets.Dataspace.Size;
sframe = 1;
num2read = dims(end)-sframe+1;
Y_h5 = single(h5read(fileName2,'/mov',[1,1,1],[512,512,num2read]));



fid2 = fopen(fileName3);
imsize = 512*512*2;
current_seek = ftell(fid2);
fseek(fid2, 0, 1);
file_length = ftell(fid2);
fseek(fid2, current_seek, -1);
frame_length = file_length/imsize;
Y_raw2 = fread(fid2,512*512*frame_length,'int16=>uint16',0,'l')';
Y_raw2 = single(2*Y_raw2);
Y_raw2 = reshape(Y_raw2,[512 512 frame_length]);
% Y_raw2 = pagetranspose(Y_raw2);
fclose(fid2);

% a = Y_raw2-Y_h5;
% b = sum(a,'all');
a = Y_raw1-Y_raw2;
b = sum(a,'all');

fig1 = figure('Name','Fig1, Y_h5','NumberTitle','off');
set(gcf,'Position',[5 35 350 350]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(Y_h5(:,:,1).^(1/3));
colormap(gray); axis off equal;

fig2 = figure('Name','Fig2, Y_raw1','NumberTitle','off');
set(gcf,'Position',[5+355 35 350 350]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(Y_raw1(:,:,1).^(1/3));
colormap(gray); axis off equal;

fig3 = figure('Name','Fig3, Y_raw2','NumberTitle','off');
set(gcf,'Position',[5+355*2 35 350 350]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
set (gca,'position',[0.01,0.01,0.98,0.98] )
imagesc(Y_raw2(:,:,1).^(1/3));
colormap(gray); axis off equal;



