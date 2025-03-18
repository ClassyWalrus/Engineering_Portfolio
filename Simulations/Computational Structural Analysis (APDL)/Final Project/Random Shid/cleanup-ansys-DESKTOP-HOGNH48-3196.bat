@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 2856)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 18360)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 18140)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 3196)

del /F cleanup-ansys-DESKTOP-HOGNH48-3196.bat
