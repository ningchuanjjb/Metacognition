::To hide command in cmd window
@echo off 
set condaRoot=C:\Users\JJB\Anaconda3
call %condaRoot%\Scripts\activate.bat
call conda activate suite2p
@echo on
::python hello.py
::python temp1.py
python C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\forSuite2p\script_jjb.py
::python script_jjb.py
::@echo off
pause