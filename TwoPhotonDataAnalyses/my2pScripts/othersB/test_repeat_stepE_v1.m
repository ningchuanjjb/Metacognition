clear 
close all

targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)

stepE_Name_v = autoGetFunName_myScripts('stepE', targetPATH);
fun_stepE = str2func(stepE_Name_v);


if_glm_12B0_24B1 = 1; %#ok<*NASGU>
fun_stepE();

if_glm_12B0_24B1 = 0;
fun_stepE();