close all
clear all %#ok<CLALL>

targetPATH1 = 'C:\ASDROOT\STUDY\MonkeyDataAnalysis\EDF';
cd(targetPATH1)

dataPath = [targetPATH1 '\edfData'];
datasfile= dir([dataPath '\*.EDF']);
exeFileName = [targetPATH1 '\edf2asc\edf2asc.exe -ntime_chcek'];
datasNum = length(datasfile);
for i = 1: datasNum
    oldFileName = [dataPath '\' datasfile(i).name] ;
    cmd = [exeFileName 32 oldFileName];
    system(cmd);
end
