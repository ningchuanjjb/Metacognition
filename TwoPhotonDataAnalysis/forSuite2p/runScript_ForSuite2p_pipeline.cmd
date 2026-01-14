::To hide command in cmd window
@echo off 
set condaRoot=C:\Users\JJB\Anaconda3
call %condaRoot%\Scripts\activate.bat
call conda activate suite2p
::@echo on

::python hello.py  > log.txt 2>&1
::python C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\forSuite2p\hello.py
::python temp1.py
::python C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\forSuite2p\script_jjb_pipeline.py 2>&1 >log.txt 
::@echo off

::powershell "dir | tee test.txt"
::powershell "python C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\forSuite2p\hello.py | tee log.txt"
powershell "python C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\forSuite2p\script_jjb_pipeline.py | tee log.txt -append"

::pause
exit