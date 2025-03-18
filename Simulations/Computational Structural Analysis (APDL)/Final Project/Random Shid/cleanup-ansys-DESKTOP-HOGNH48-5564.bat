@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 632)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 19764)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 16036)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 5564)

del /F cleanup-ansys-DESKTOP-HOGNH48-5564.bat
