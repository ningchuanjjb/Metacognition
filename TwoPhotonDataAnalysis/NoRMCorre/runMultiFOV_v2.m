clear
close all


targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre';
cd(targetPATH)


currentSession_multi = string;
% currentSession_multi = cell(1,1);

currentSession_multi = [currentSession_multi; '113Recording_20230426A_Ding_Site16'];
currentSession_multi = [currentSession_multi; '113Recording_20230427A_Ding_Site16_sameFOV0426'];
currentSession_multi = [currentSession_multi; '113Recording_20230502A_Ding_Site13'];
currentSession_multi = [currentSession_multi; '113Recording_20230503A_Ding_Site13_sameFOV0502'];
currentSession_multi = [currentSession_multi; '113Recording_20230504A_Ding_Site02'];
currentSession_multi = [currentSession_multi; '113Recording_20230508A_Ding_Site02_sameFOV0509'];
currentSession_multi = [currentSession_multi; '113Recording_20230509A_Ding_Site02'];% 660000 frames, easy to crash

currentSession_multi = [currentSession_multi; '113Recording_20230510A_Ding_Site05_sameFOV0511'];
currentSession_multi = [currentSession_multi; '113Recording_20230510B_Ding_Site05_sameFOV0511'];
currentSession_multi = [currentSession_multi; '113Recording_20230511A_Ding_Site05'];
currentSession_multi = [currentSession_multi; '113Recording_20230512A_Ding_Site09'];
currentSession_multi = [currentSession_multi; '113Recording_20230513A_Ding_Site09_sameFOV0512'];

currentSession_multi = [currentSession_multi; '113Recording_20230515A_Ding_Site24_sameFOV0516'];
currentSession_multi = [currentSession_multi; '113Recording_20230516A_Ding_Site24'];
currentSession_multi = [currentSession_multi; '113Recording_20230517A_Ding_Site16B'];
currentSession_multi = [currentSession_multi; '113Recording_20230522A_Ding_Site05B'];
currentSession_multi = [currentSession_multi; '113Recording_20230525A_Ding_Site09B'];
currentSession_multi = [currentSession_multi; '113Recording_20230526A_Ding_Site09B_sameFOV0525'];

currentSession_multi = [currentSession_multi; '113Recording_20230527A_Ding_Site02B'];
currentSession_multi = [currentSession_multi; '113Recording_20230528A_Ding_Site02B_sameFOV0527'];
currentSession_multi = [currentSession_multi; '113Recording_20230530A_Ding_Site05C'];
currentSession_multi = [currentSession_multi; '113Recording_20230531A_Ding_Site05C_sameFOV0530'];
currentSession_multi = [currentSession_multi; '113Recording_20230601A_Ding_Site13B'];
currentSession_multi = [currentSession_multi; '113Recording_20230602A_Ding_Site13B_sameFOV0601'];

currentSession_multi = [currentSession_multi; '113Recording_20230604A_Ding_Site07'];
currentSession_multi = [currentSession_multi; '113Recording_20230605A_Ding_Site07_sameFOV0604'];
currentSession_multi = [currentSession_multi; '113Recording_20230612A_Ding_Site14'];
currentSession_multi = [currentSession_multi; '113Recording_20230614A_Ding_Site05D'];
currentSession_multi = [currentSession_multi; '113Recording_20230615A_Ding_Site05D_sameFOV0614'];
currentSession_multi = [currentSession_multi; '113Recording_20230619A_Ding_Site02C'];
currentSession_multi = [currentSession_multi; '113Recording_20230620A_Ding_Site05E'];

currentSession_multi(1) = [];


runFullPipeline_Name_v = autoGetFunName('runFullPipeline', targetPATH);
fprintf('Now runing %s.\n', runFullPipeline_Name_v);
fun_runFullPipeline = str2func(runFullPipeline_Name_v);

num_FOV = length(currentSession_multi);

t0 = tic;
for tempi=1:num_FOV
    fun_runFullPipeline(currentSession_multi{tempi});
end
fprintf('Time of the whole pipeline is %.1f secs.\n',toc(t0));
fprintf('num_FOV=%d.\n',num_FOV);

