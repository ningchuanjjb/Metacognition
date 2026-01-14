::To hide command in cmd window
@echo off 
set condaRoot=C:\Users\JJB\Anaconda3
call %condaRoot%\Scripts\activate.bat
call conda activate suite2p
@echo on
suite2p
::@echo off
::pause