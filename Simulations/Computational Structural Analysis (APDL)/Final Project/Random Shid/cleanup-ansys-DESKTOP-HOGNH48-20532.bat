@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 11896)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 3224)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 11096)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 20532)

del /F cleanup-ansys-DESKTOP-HOGNH48-20532.bat
