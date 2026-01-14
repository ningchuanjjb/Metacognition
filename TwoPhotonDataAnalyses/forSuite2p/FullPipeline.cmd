::To hide command in cmd window

@echo off 
cd C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\NoRMCorre
@echo on

"C:\Program Files\Polyspace\R2021a\bin\matlab.exe" -batch "run('runNoRM.m')"


@echo off 
cd C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\ForSuite2p
set condaRoot=C:\Users\JJB\Anaconda3
call %condaRoot%\Scripts\activate.bat
call conda activate suite2p
@echo on


python script_jjb.py
suite2p
