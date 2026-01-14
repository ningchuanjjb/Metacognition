close all
clear all

targetPATH1 = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\EDF';
cd(targetPATH1)

dataPath = [targetPATH1 '\edfData'];
datasfile= dir([dataPath '\*.EDF']);
exeFileName = ['targetPATH1' 'edf2asc.exe -ntime_chcek'];
datasNum = length(datasfile);
% for j =3:datasNum
%     fileNamePath = [dataPath '\' datasfile(j).name ];
%     file = dir(fileNamePath);
%     filesNum = length(file);
%     for i =3:filesNum
for i = 1: datasNum
%         fileJIPath = [fileNamePath file(i).name '\'];
%         fileji = dir([fileJIPath '*.edf']);
       oldFileName = [dataPath '\' datasfile(i).name] ;
       cmd = [exeFileName 32 oldFileName];
       system(cmd);    
    end
