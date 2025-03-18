@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 14152)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 968)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 3120)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 1052)

del /F cleanup-ansys-DESKTOP-HOGNH48-1052.bat
