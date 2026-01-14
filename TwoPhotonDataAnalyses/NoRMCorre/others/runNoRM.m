clear
close all

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)

NoRMCorreName = 'myNoRMCorre';

NoRMCorreName_v = autoGetFunName(NoRMCorreName, [targetPATH]);
fprintf('Now runing %s.\n', NoRMCorreName_v);
pause(1);
fun_NoRMCorre = str2func(NoRMCorreName_v);
fun_NoRMCorre();
